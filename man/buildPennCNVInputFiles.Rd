\author{Hyuna Yang}
\name{buildPennCNVInputFiles}
\alias{buildPennCNVInputFiles}
\title{Build PennCNV input files for the Mouse Diversity Array}
\description{
    This function creates the LRR/BAF and PFB inputs that are needed to run PennCNV
}
\usage{
buildPennCNVInputFiles(
    outdir = getwd(), allowOverwrite = FALSE,
    genotypes, snpProbeInfo, snpInfo, snpReferenceDistribution = NULL,
    invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans"), celFiles = expandCelFiles(getwd()),
    chromosomes = c(1:19, "X", "Y", "M"),
    chromosomeRenameMap = list(X = 20, Y = 21, M = 22),
    cacheDir = tempdir(),
    retainCache = FALSE, verbose = FALSE,
    probesetChunkSize = 1000)
}

\arguments{
    \item{outdir}{the directory to write PennCNV input files to}
    
    \item{allowOverwrite}{
        determines if this function is allowed to overwrite exising PennCNV
        input
    }
    
    \item{genotypes}{
        genotype data as returned by the mouseDivGenotype function
        genotypes must match the order of the SNPs in snpProbeInfo
    }
    
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are performing CNV on.
        The data-frame should contain the following components:
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales
    }
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are performing CNV on.
        The data-frame should contain the following components:

        \code{snpId}: the ID of this SNP

        \code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}

        \code{positionBp}: the chromosome position in base pairs

        \code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygotes

        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is in the pseudoautosomal region (PAR). This should only be set to \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
    }
    
    \item{snpReferenceDistribution}{
        A numeric vector which is the reference distribution used for the
        quantile normalization. If there is an evidence that the density
        distribution of samples follows more than one distribution or different
        from that of classical laboratory strain, this normalization should not
        be used and leave as a default i.e. \code{referenceDistribution = NULL}.
    }
    
    \item{invariantProbeInfo}{
        A data-frame containing probe-level data for the invariants we are
        performing CNV on, or a list of data frames where each data frame in the
        list represents an invariant group which should be normalized seperately
        from the others. The data-frame should contain the following components:
        
        \code{probesetId}: the ID of the probeset that this probe belongs to. These IDs should correspond to the same component in the \code{invariantProbesetInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales
    }
    \item{invariantProbesetInfo}{
        A data-frame containing probeset-level data for the invariants we are
        performing CNV on, or a list of data frames where each data frame in the
        list represents an invariant group which should be normalized seperately
        from the others.
        The data-frame should contain the following components:
        
        \code{probesetId}: the ID of this probeset
        
        \code{chrId}: the chromosome that this probeset belongs to. Like \dQuote{15} or \dQuote{X}
        
        \code{positionBp}: the chromosome position in base pairs
    }
    \item{invariantReferenceDistribution}{
        A numeric vector which is the reference distribution used for the
        quantile normalization, or a list of vectors where each data frame in
        the list represents an invariant group which should be normalized
        seperately from the others.
    }
  \item{transformMethod}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
  \item{celFiles}{the vector of CEL files to process}
  \item{chromosomes}{character vector of chromosome IDs that you want to CNV. These should correspond to the values in \code{snpInfo$chrId}. Eg: \code{chromosomes=c(1 : 19, "X", "Y", "M")}.}
  \item{chromosomeRenameMap}{A list that allows you to rename chromosomes for the PennCNV files. This is needed because PennCNV uses all numeric chromosome values where this packages uses "X", "Y" and "M"}
  \item{cacheDir}{directory used to store intermediate results. By default a temporary directory is used. You can reuse a cache directory to save time in the normalization steps (this function will only perform CEL file normalization if it does not find the data that it needs in the cache)}
  \item{retainCache}{defaults to \code{FALSE}. If \code{FALSE} the data cache is deleted before returning. Setting this value to \code{TRUE} could save processing time on subsequent calls}
  \item{verbose}{if set to \code{TRUE} status information will be printed while the function executes}
  \item{probesetChunkSize}{determines how many SNPs will be processed at a time (in a single chunk). Keeping this number low reduces memory utilization and makes paralellism more fine-grained at the cost of managing more chunk objects. The default value provides a reasonable tradeoff}  
}

\value{
    No return value. Instead PennCNV input files are written to \code{outdir}
}

\examples{
\dontrun{
baseDir <- "."
celfiledir <- file.path(baseDir, "celFiles/")
celfiles <- expandCelFiles(celfiledir)

# create an output cache so that we only normalize CEL files once
outcache <- file.path(baseDir, "outcache")
if(!file.exists(outcache))
{
    dir.create(outcache)
}

outdir <- file.path(baseDir, "cnv-out")
if(file.exists(outdir))
{
    unlink(outdir, recursive = TRUE)
}
dir.create(outdir)

# load the genotypes that we saved from a previous call to mouseDivGeno
load("genoVinoResult-2010-08-10_17.30.33.59530.RData")
genos <- genoVinoResult$geno
rm(genoVinoResult)


################################################################################
# This command loads the snpProbeInfo, snpInfo and reference objects
# you should download this from the CGD website in order to be able to run
# through this example
################################################################################
load("MouseDivData.RData")

buildPennCNVInputFiles(
    outdir = outdir,
    allowOverwrite = FALSE,
    genotypes = genos,
    snpProbeInfo = snpProbeInfo,
    snpInfo = snpInfo,
    snpReferenceDistribution = snpReferenceDistribution,
    invariantProbeInfo = invariantProbeInfo,
    invariantProbesetInfo = invariantProbesetInfo,
    invariantReferenceDistribution = invariantReferenceDistribution,
    transformMethod = "CCStrans",
    celFiles = celfiles,
    chromosomes = c("19", "X", "Y"),
    cacheDir = outcache,
    retainCache = retainCache,
    verbose = TRUE)
  }
}

\keyword{CNV}
