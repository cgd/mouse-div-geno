\author{Hyuna Yang}
\name{mouseDivGenotype}
\alias{mouseDivGenotype}
\title{Genotype the Mouse Diversity Array}
\description{
  This function reads the .CEL files, normalizes the intensities, genotypes and vinotypes the Mouse Diversity Array. It also provides the input files used for the pennCNV (Wang, et al, 2007, a software to detect the copy number variance. Refer to CNV.R for detail information for detecting CNV.).  

  This function is highly tailored to Mouse Diversity Array. All annotation data required to use this function can be obtained from \url{http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. Coefficients to normalize the array (\code{snpProbeInfo$correction} and \code{referenceDistribution}) or hint file used to genotype the array (\code{snpInfo$snpHetHint}) were calculated based on 350 training .CEL files, and they can be obtained from the \url{http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}.  

}
\usage{
mouseDivGenotype(
    snpProbeInfo, snpInfo, referenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans"),
    celFiles = expandCelFiles(getwd()),
    chromosomes = c(1:19, "X", "Y", "M"), cacheDir = tempdir(),
    retainCache = FALSE, verbose = FALSE, cluster = NULL,
    probesetChunkSize = 1000, processResultsFunction = NULL)
}

\arguments{
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        \itemize{
            \item{\code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for B allele probes}
            \item{\code{snpId}: the ID of the SNP that this probe belongs to. These IDs should correspond to the same component in the \code{snpInfo} argument}
            \item{\code{probeIndex}: the index for the probe (used to get intensity data)}
            \item{\code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales}
        }
    }
    
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        \itemize{
            \item{\code{snpId}: the ID of this SNP}
            \item{\code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}}
            \item{\code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygotes}
            \item{\code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is in the pseudoautosomal region (PAR). This should only be set to \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}}
        }
    }
    
  \item{referenceDistribution}{Reference distribution used for the quantile normalization. If there is an evidence that the density distribution of samples follows more than one distribution or different from that of classical laboratory strain, this normalization should not be used and leave as a default i.e. \code{referenceDistribution = NULL}.}
  \item{transformMethod}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
  
  \item{celFiles}{this parameter should be one of the following:
        \itemize{
            \item{a vector of CEL files to process. In this case sample gender will be inferred}
            \item{a data-frame with a character "fileName" column and a boolean "isMale" column.}
        }
  }
  \item{chromosomes}{character vector of chromosome IDs that you want to genotype. These should correspond to the values in \code{snpInfo$chrId}. Eg: \code{chromosomes=c(1 : 19, "X", "Y", "M")}.}
  \item{cacheDir}{directory used to store intermediate results. By default a temporary directory is used. You can reuse a cache directory to save time in the normalization steps (this function will only perform CEL file normalization if it does not find the data that it needs in the cache)}
  \item{retainCache}{defaults to \code{FALSE}. If \code{FALSE} the data cache is deleted before returning. Setting this value to \code{TRUE} could save processing time on subsequent calls}
  \item{verbose}{if set to \code{TRUE} status information will be printed while the function executes}
  \item{cluster}{a cluster object as created by the \pkg{SNOW} package. If this value is non-NULL then this function will run genotyping computations in parallel on the given cluster}
  \item{probesetChunkSize}{determines how many SNPs will be processed at a time (in a single chunk). Keeping this number low reduces memory utilization and makes paralellism more fine-grained at the cost of managing more chunk objects. The default value provides a reasonable tradeoff}  
  \item{processResultsFunction}{an optional function argument which can process results chunk-by-chunk. Typically you should use the \code{createAppendResultsToCSVFunction} to create this function (see example below). Normally this parameter will just provide a convenient way to save results to a CSV file, but it can also be an important way to save memory if you are processing many CEL files in a single run (because results are processed chunk-by-chunk there is no longer a need to accumulate results for a single return value. Instead this function just returns NULL when \code{processResultsFunction} is non-NULL)}
}

\value{
    If \code{processResultsFunction} is non-NULL then this function returns a NULL value because it is assumed that the \code{processResultsFunction} is saving the results chunk-by-chunk allowing us to avoid the memory requirements imposed by accumulating all of the results into a return list. Assuming \code{processResultsFunction} is NULL this function returns an list object indexed by chromosome ID, containing another list indexed by "geno", "vino" or "conf", containing a matrix of values. So, for instance, \code{retVal[["X"]][["geno"]][3, 5]} will get you the genotype value for the third SNP on the X chromosome of the fifth CEL file (sample). It is probably simplest to use the \code{processResultsFunction} parameter so save your results to file unless you plan on further processing the genotype/vinotype data in R. The following describes the "geno", "vino" and "conf" list items
    \item{geno}{a matrix with a row per SNP and a column per CEL file. -1 = No call, 1 = AA, 2 = AB, and 3 = BB}
    \item{vino}{a matrix with a row per SNP and a column per CEL file. 1 indicates VINO, any other value indicates no-VINO}
    \item{conf}{a matrix with a row per SNP and a column per CEL file. These values are confidence scores based on Mohalanobis distance with chi-squre distribution approximation. Smaller confidence score implies less reliable data.}
}

\examples{
\dontrun{
baseDir <- "."
celfiledir <- file.path(baseDir, "celFiles/")

# create an output cache so that we only normalize CEL files once
outcache <- file.path(baseDir, "outcache")
dir.create(outcache)

################################################################################
# This command loads the snpProbeInfo, snpInfo and reference objects
# you should download this from the CGD website in order to be able to run
# through this example
################################################################################
load("MouseDivData.RData")

# genotype in serial and save the results to "serialGenoResults.csv"
cat("Running mouseDivGenotype in serial and saving to serialGenoResults.csv\n")
mouseDivGenotype(
    snpProbeInfo = snpProbeInfo, snpInfo = snpInfo, referenceDistribution=reference,
    transformMethod="CCStrans", chromosomes=c("19", "X", "Y"), celFiles=expandCelFiles(celfiledir),
    cacheDir = outcache, verbose=T, cluster=NULL,
    processResultsFunction = createAppendResultsToCSVFunction("serialGenoResults.csv"))

# Run in parallel on 2 cores (increase cores if you have more on your machine)
# and save the results to "parallelGenoResults.csv"
library("snow")
setDefaultClusterOptions(outfile = "/dev/tty")
numCores <- 2
cl <- makeSOCKcluster(rep("localhost", numCores))
clusterEvalQ(cl, library("MouseDivGeno"))

cat("Running MouseDivGenotype in parallel and saving to parallelGenoResults.csv\n")
mouseDivGenotype(
    snpProbeInfo = snpProbeInfo, snpInfo = snpInfo, referenceDistribution=reference,
    transformMethod="CCStrans", chromosomes=c("19", "X", "Y"), celFiles=expandCelFiles(celfiledir),
    cacheDir = outcache, verbose=T, cluster=cl,
    processResultsFunction = createAppendResultsToCSVFunction("parallelGenoResults.csv"))

# Run again in parallel but this time save the resulting R object to file rather
# than using the processResultsFunction argument
cat("Running MouseDivGenotype in parallel and saving to genoVinoResult-TIME_STAMP.RData.csv\n")
genoVinoResult <- mouseDivGenotype(
    snpProbeInfo = snpProbeInfo, snpInfo = snpInfo, referenceDistribution=reference,
    transformMethod="CCStrans", chromosomes=c("19", "X", "Y"), celFiles=expandCelFiles(celfiledir),
    cacheDir = outcache, verbose=T, cluster=cl,
    processResultsFunction = NULL)
timeStampFileName <- file.path(
    baseDir,
    format(Sys.time(), "genoVinoResult-%Y-%m-%d_%H.%M.%OS5.RData"))
save(genoVinoResult, file = timeStampFileName)

stopCluster(cl)
  }
}

\keyword{genotype}
