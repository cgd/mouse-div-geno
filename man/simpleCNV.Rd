\author{Hyuna Yang and Keith Sheppard}
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
        celFiles = getwd(),
        referenceCelFile,
        chromosomes = c(1:19, "X", "Y", "M"),
        verbose = FALSE,
        cluster = NULL,
        summaryOutputFile = NULL)
}

\arguments{
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are performing CNV on.
        The data-frame should contain the following components:
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for
        B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs
        should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied
        by adding it to probe mean intesity vales
    }
    
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are performing
        CNV on. The data-frame should contain the following components:
        
        \code{snpId}: the ID of this SNP
        
        \code{chrId}: the chromosome that this SNP belongs to.
        Like \dQuote{15} or \dQuote{X}
        
        \code{positionBp}: the chromosome position in base pairs
        
        \code{snpHetHint} (optional): provides a hint for the normalized mean
        intensity of heterozygous SNPs
        
        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is
        in the pseudoautosomal region (PAR). This should only be set to
        \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
    }
    
    \item{snpReferenceDistribution}{
        Reference distribution used for the quantile normalization. If there is
        an evidence that the density distribution of samples follows more than
        one distribution or different from that of classical laboratory strain,
        this normalization should not be used and leave as a default i.e.
        \code{referenceDistribution = NULL}.
    }
    
    \item{invariantProbeInfo}{
        A data-frame containing probe-level data for the invariants we are
        performing CNV on, or a list of data frames where each data frame in the
        list represents an invariant group which should be normalized separately
        from the others. The data-frame should contain the following components:
        
        \code{probesetId}: the ID of the probeset that this probe belongs to.
        These IDs should correspond to the same component in the
        \code{invariantProbesetInfo} argument

        \code{probeIndex}: the index for the probe (used to get intensity data)

        \code{correction} (optional): if present this correction will be applied
        by adding it to probe mean intensity vales
    }
    
    \item{invariantProbesetInfo}{
        A data-frame containing probeset-level data for the invariants we are
        performing CNV on, or a list of data frames where each data frame in the
        list represents an invariant group which should be normalized seperately
        from the others.
        The data-frame should contain the following components:
        
        \code{probesetId}: the ID of this probeset

        \code{chrId}: the chromosome that this probeset belongs to. Like
        \dQuote{15} or \dQuote{X}

        \code{positionBp}: the chromosome position in base pairs
    }
    
    \item{invariantReferenceDistribution}{
        A numeric vector which is the reference distribution used for the
        quantile normalization, or a list of vectors where each data frame in
        the list represents an invariant group which should be normalized
        separately from the others.
    }
  
    \item{celFiles}{
        either a vector of CEL file name(s) or directory name(s). Any directories
        will be expanded to the listing of CEL files found in the directory
    }
  
    \item{referenceCelFile}{
        the CEL file used as a reference when testing for CNV gains and losses
    }
    
    \item{chromosomes}{
        character vector of chromosome IDs that you want to CNV. These should
        correspond to the values in \code{snpInfo$chrId}.
        Eg: \code{chromosomes=c(1 : 19, "X", "Y", "M")}.
    }
    \item{verbose}{
        if set to \code{TRUE} status information will be printed while the
        function executes
    }
    \item{cluster}{
        a cluster object as created by the \pkg{SNOW} package. If this value is
        non-NULL then this function will run CNV computations in parallel on the
        given cluster
    }
    \item{summaryOutputFile}{
        the file that you would like to write a summary-level
        CNV table to. NULL (the default) indicates that no summary file should be
        generated. If a valid filename is given the output table will have the
        following columns:
        Name, Status, Chromosome, StartPosition, EndPosition,
        Number of Probe sets, Size, Number of SNP probe sets,
        Number of exon1 sets, Number of exon2 sets,
        Mean intensity of reference sample, Mean intensity of sample
    }
}

\value{
    a list containing per-chromosome CNV results. Each list item will be a matrix
    where the columns correspond to samples (CEL files) and the rows correspond
    to invariant probesets or SNPs. In the returned matrix list a value of
    2 indicates no copy change, a value of 1 indicates a copy loss with respect
    to the reference and a value of 3 indicates a copy gain with respect to the
    reference
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

# so that we can print timing information
library("time")

# First so simple CNV analysis without any clustering
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
    verbose = TRUE)
cat("time it took us to do non-clustered CNV\n")
timeReport(startTime)

timeStampFileName <- file.path(
    baseDir,
    format(Sys.time(), "simpleCNVResultNotClustered-%Y-%m-%d_%H.%M.%OS5.RData"))
save(simpleCNVResultNotClustered, file = timeStampFileName)

# Do simple CNV analysis using a "local" SNOW cluster
library("snow")
setDefaultClusterOptions(outfile = "/dev/tty")
numCores <- 2
cl <- makeSOCKcluster(rep("localhost", numCores))
clusterEvalQ(cl, library("MouseDivGeno"))

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
  }
}

\keyword{CNV}
