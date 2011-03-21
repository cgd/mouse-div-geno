\author{Hyuna Yang}
\name{mouseDivGenotype}
\alias{mouseDivGenotype}
\alias{mouseDivGenotypeTab}
\title{Genotype the Mouse Diversity Array}
\description{
  This function reads the .CEL files, normalizes the intensities, genotypes and vinotypes the Mouse Diversity Array. It also provides the input files used for the pennCNV (Wang, et al, 2007, a software to detect the copy number variance. Refer to CNV.R for detail information for detecting CNV.).  

  This function is highly tailored to Mouse Diversity Array. All annotation data required to use this function can be obtained from \url{http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. Coefficients to normalize the array (\code{snpProbeInfo$correction} and \code{referenceDistribution}) or hint file used to genotype the array (\code{snpInfo$snpHetHint}) were calculated based on 350 training .CEL files, and they can be obtained from the \url{http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}.  

}
\usage{
mouseDivGenotype(
        snpProbeInfo, snpInfo, referenceDistribution = NULL,
        transformMethod = c("CCStrans", "MAtrans"),
        celFiles = expandCelFiles(getwd()), isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        verbose = FALSE, cluster = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_")
mouseDivGenotypeTab(
        snpIntensityFile, snpInfo,
        transformMethod = c("CCStrans", "MAtrans"),
        isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        verbose = FALSE, cluster = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_")
}

\arguments{
    \item{snpIntensityFile}{a tab delimited file with columns for "sampleID", "snpId", "A", "B" where the A and B columns contain log-normalized summarized per-SNP intensities}
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales
    }
    
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{snpId}: the ID of this SNP
        
        \code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}
        
        \code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygous SNPs
        
        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is in the pseudoautosomal region (PAR). This should only be set to \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
    }
    
  \item{referenceDistribution}{Reference distribution used for the quantile normalization. If there is an evidence that the density distribution of samples follows more than one distribution or different from that of classical laboratory strain, this normalization should not be used and leave as a default i.e. \code{referenceDistribution = NULL}.}
  \item{transformMethod}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
  
  \item{celFiles}{this parameter should be a vector of CEL files to process}
  \item{isMale}{a logical vector the same length as celFiles or NULL if gender should be inferred from the CEL file intensities}
  \item{confScoreThreshold}{if confidence scores fall below this threshold we assign a genotype value of -1}
  \item{chromosomes}{character vector of chromosome IDs that you want to genotype. These should correspond to the values in \code{snpInfo$chrId}. Eg: \code{chromosomes=c(1 : 19, "X", "Y", "M")}.}
  \item{verbose}{if set to \code{TRUE} status information will be printed while the function executes}
  \item{cluster}{a cluster object as created by the \pkg{SNOW} package. If this value is non-NULL then this function will run genotyping computations in parallel on the given cluster}
  \item{probesetChunkSize}{determines how many SNPs will be processed at a time (in a single chunk). Keeping this number low reduces memory utilization and makes paralellism more fine-grained at the cost of managing more chunk objects. The default value provides a reasonable tradeoff}
  \item{outputDir}{an optional directory name argument which will cause genotype, vinotype and confidence results to be saved to a comma separated file rather than being returned as an R object. This parameter will provide a convenient way to save results to a CSV file, but it can also be an important way to save memory if you are processing many CEL files in a single run (because results are processed chunk-by-chunk there is no longer a need to accumulate results for a single return value.)}
  \item{outputFilePrefix}{this parameter is only used where outputDir is non-NULL. Three files will be saved in outputDir (for genotypes, vinotypes and confidence scores). This parameter determines the prefix used in those three file names.}
}

\value{
    If (and only if) the \code{outputDir} parameter is non-NULL then this function returns a
    list containing only \code{isMale} value because
    \code{outputDir} will contain data for \code{geno}, \code{vino} and
    \code{conf}. If \code{outputDir} is NULL however the return list
    will contain all of the following components:
    
    \item{geno}{a matrix with a row per SNP and a column per CEL file. -1 = No call, 1 = AA, 2 = AB, and 3 = BB}
    \item{vino}{a matrix with a row per SNP and a column per CEL file. 1 indicates VINO, 0 indicates no-VINO}
    \item{conf}{a matrix with a row per SNP and a column per CEL file. These values are confidence scores based on Mohalanobis distance with chi-squre distribution approximation. Smaller confidence score implies less reliable data.}
    \item{isMale}{A per-CEL-file vector of logicals indicating which were genotyped as male. If gender was supplied via the CEL files parameter then this value will match that otherwise it will be the inferred gender values.}
}

\examples{
\dontrun{
baseDir <- "."

# the directory that contains all of our CEL files
celfiledir <- file.path(baseDir, "celFiles/")

################################################################################
# This command loads the snpProbeInfo, snpInfo and reference objects
# you should download this from the CGD website in order to be able to run
# through this example
################################################################################
load("MouseDivData.RData")

# genotype in serial and save the results to "." (the current directory)
cat("Running mouseDivGenotype in serial and saving to serialGenoVinoResults.csv\n")
mouseDivGenotype(
    snpProbeInfo = snpProbeInfo, snpInfo = snpInfo, referenceDistribution=reference,
    chromosomes=c("19", "X", "Y"), celFiles=expandCelFiles(celfiledir),
    verbose = T, resultsOutputFile = ".")

# Run in parallel and this time save the resulting raw R object to file rather
# than using the processResultsFunction argument to create a CSV file
library("snow")
numCores <- 2
cl <- makeSOCKcluster(rep("localhost", numCores))
clusterEvalQ(cl, library("MouseDivGeno"))

cat("Running MouseDivGenotype in parallel and saving to parallelGenoVinoResult.RData.csv\n")
genoVinoResult <- mouseDivGenotype(
    snpProbeInfo = snpProbeInfo, snpInfo = snpInfo, referenceDistribution=reference,
    chromosomes=c("19", "X", "Y"), celFiles=expandCelFiles(celfiledir),
    verbose=T, cluster=cl)
save(genoVinoResult, file = "parallelGenoVinoResult.RData")

stopCluster(cl)
  }
}

\keyword{genotype}
