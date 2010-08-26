\author{Hyuna Yang}
\name{simpleCNV}
\alias{simpleCNV}
\title{Simplified CNV Analysis for the Mouse Diversity Array}
\description{
    This function performs a simplified CNV analysis which does not require
    the use of PennCNV (compare to buildPennCNVInputFiles)
}
\usage{
simpleCNV(
    snpProbeInfo, snpInfo, snpReferenceDistribution = NULL,
    invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution = NULL,
    celFiles = expandCelFiles(getwd()),
    referenceCelFile,
    chromosomes = c(1:19, "X", "Y", "M"),
    verbose = FALSE,
    cluster = NULL)
}

\arguments{
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are performing CNV on.
        The data-frame should contain the following components:
        \itemize{
            \item{\code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for B allele probes}
            \item{\code{snpId}: the ID of the SNP that this probe belongs to. These IDs should correspond to the same component in the \code{snpInfo} argument}
            \item{\code{probeIndex}: the index for the probe (used to get intensity data)}
            \item{\code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales}
        }
    }
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are performing CNV on.
        The data-frame should contain the following components:
        \itemize{
            \item{\code{snpId}: the ID of this SNP}
            \item{\code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}}
            \item{\code{positionBp}: the chromosome position in base pairs}
            \item{\code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygotes}
            \item{\code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is in the pseudoautosomal region (PAR). This should only be set to \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}}
        }
    }
  \item{snpReferenceDistribution}{Reference distribution used for the quantile normalization. If there is an evidence that the density distribution of samples follows more than one distribution or different from that of classical laboratory strain, this normalization should not be used and leave as a default i.e. \code{referenceDistribution = NULL}.}
  
    \item{invariantProbeInfo}{
        A data-frame containing probe-level data for the invariants we are performing CNV on.
        The data-frame should contain the following components:
        \itemize{
            \item{\code{probesetId}: the ID of the probeset that this probe belongs to. These IDs should correspond to the same component in the \code{invariantProbesetInfo} argument}
            \item{\code{probeIndex}: the index for the probe (used to get intensity data)}
            \item{\code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales}
        }
    }
    \item{invariantProbesetInfo}{
        A data-frame containing probeset-level data for the invariants we are performing CNV on.
        The data-frame should contain the following components:
        \itemize{
            \item{\code{probesetId}: the ID of this probeset}
            \item{\code{chrId}: the chromosome that this probeset belongs to. Like \dQuote{15} or \dQuote{X}}
            \item{\code{positionBp}: the chromosome position in base pairs}
        }
    }
  \item{invariantReferenceDistribution}{Reference distribution used for the quantile normalization.}
  
  \item{celFiles}{the vector of CEL files to process}
  \item{referenceCelFile}{the CEL file used as a reference when testing for CNV gains and losses}
  \item{chromosomes}{character vector of chromosome IDs that you want to CNV. These should correspond to the values in \code{snpInfo$chrId}. Eg: \code{chromosomes=c(1 : 19, "X", "Y", "M")}.}
  \item{verbose}{if set to \code{TRUE} status information will be printed while the function executes}
  \item{cluster}{a cluster object as created by the \pkg{SNOW} package. If this value is non-NULL then this function will run CNV computations in parallel on the given cluster}
}

\value{
    a list containing per-chromosome CNV results. Each list item will be a matrix
    where the columns correspond to samples (CEL files) and the rows
}

\examples{
\dontrun{
baseDir <- "."
celfiledir <- file.path(baseDir, "celFiles/")
celfiles <- expandCelFiles(celfiledir)


################################################################################
# This command loads the snpProbeInfo, snpInfo and reference objects
# you should download this from the CGD website in order to be able to run
# through this example
################################################################################
load("MouseDivData.RData")

# Do simple CNV analysis using a "local" SNOW cluster
library("snow")
setDefaultClusterOptions(outfile = "/dev/tty")
numCores <- 2
cl <- makeSOCKcluster(rep("localhost", numCores))
clusterEvalQ(cl, library("MouseDivGeno"))

library("time")

startTime <- getTime()
simpleCNVResultClustered <- simpleCNV(
    snpProbeInfo = snpProbeInfo,
    snpInfo = snpInfo,
    snpReferenceDistribution = snpReferenceDistribution,
    invariantProbeInfo = invariantProbeInfo,
    invariantProbesetInfo = invariantProbesetInfo,
    invariantReferenceDistribution = invariantReferenceDistribution,
    celFiles = celfiles,
    referenceCelFile = file.path(celfiledir, "SNP_mDIV_A1-1_081308.CEL"),
    chromosomes = c("19", "X", "Y"),
    verbose = TRUE,
    cluster = cl)
cat("time it took us to do clustered CNV\n")
timeReport(startTime)

timeStampFileName <- file.path(
    baseDir,
    format(Sys.time(), "simpleCNVResultClustered-%Y-%m-%d_%H.%M.%OS5.RData"))
save(simpleCNVResultClustered, file = timeStampFileName)
rm(simpleCNVResultClustered)

stopCluster(cl)

# Do simple CNV analysis without a cluster
startTime <- getTime()
simpleCNVResultNotClustered <- simpleCNV(
    snpProbeInfo = snpProbeInfo,
    snpInfo = snpInfo,
    snpReferenceDistribution = snpReferenceDistribution,
    invariantProbeInfo = invariantProbeInfo,
    invariantProbesetInfo = invariantProbesetInfo,
    invariantReferenceDistribution = invariantReferenceDistribution,
    celFiles = celfiles,
    referenceCelFile = file.path(celfiledir, "SNP_mDIV_A1-1_081308.CEL"),
    chromosomes = c("19", "X", "Y"),
    verbose = TRUE,
    cluster = NULL)
cat("time it took us to do non-clustered CNV\n")
timeReport(startTime)

timeStampFileName <- file.path(
    baseDir,
    format(Sys.time(), "simpleCNVResultNotClustered-%Y-%m-%d_%H.%M.%OS5.RData"))
save(simpleCNVResultNotClustered, file = timeStampFileName)
rm(simpleCNVResultNotClustered)
  }
}

\keyword{genotype}
