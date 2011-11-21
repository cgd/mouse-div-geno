\author{Hyuna Yang and Keith Sheppard}
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
        transformFunction = ccsTransform, celFiles = getwd(),
        isMale = NULL,
        chromosomes = c(1:19, "X", "Y", "M"),
        chromosomeRenameMap = list(X = "20", Y = "21", M = "22"),
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
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for
        B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs
        should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied
        by adding it to probe mean intesity vales
    }
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are performing CNV on.
        The data-frame should contain the following components:

        \code{snpId}: the ID of this SNP

        \code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15}
        or \dQuote{X}

        \code{positionBp}: the chromosome position in base pairs

        \code{snpHetHint} (optional): provides a hint for the normalized mean
        intensity of heterozygotes

        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is
        in the pseudoautosomal region (PAR). This should only be set to
        \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
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
        list represents an invariant group which should be normalized separately
        from the others. The data-frame should contain the following components:
        
        \code{probesetId}:
        the ID of the probeset that this probe belongs to. These IDs should
        correspond to the same component in the \code{invariantProbesetInfo}
        argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional):
        if present this correction will be applied by adding it to probe mean
        intensity vales
    }
    \item{invariantProbesetInfo}{
        A data-frame containing probeset-level data for the invariants we are
        performing CNV on, or a list of data frames where each data frame in the
        list represents an invariant group which should be normalized separately
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
        seperately from the others.
    }
    \item{transformFunction}{
        function used to transform SNP A/B intensities into intensity averages
        and contrasts. This function should either be ccsTransform, maTransform
        or your own function built using the CCS and MA code as an example
    }
    \item{celFiles}{
        either a vector of CEL file name(s) or directory name(s). Any directories
        will be expanded to the listing of CEL files found in the directory
    }
    \item{isMale}{
        a logical vector the same length as celFiles that you can use if gender
        is known
    }
    \item{chromosomes}{
        character vector of chromosome IDs that you want to CNV. These should
        correspond to the values in \code{snpInfo$chrId}.
        Eg: \code{chromosomes=c(1 : 19, "X", "Y", "M")}.
    }
    \item{chromosomeRenameMap}{
        A list that allows you to rename chromosomes for the PennCNV files. This
        is needed because PennCNV uses all numeric chromosome values where this
        packages uses "X", "Y" and "M"
    }
    \item{cacheDir}{
        directory used to store intermediate results. By default a temporary
        directory is used. You can reuse a cache directory to save time in the
        normalization steps (this function will only perform CEL file normalization
        if it does not find the data that it needs in the cache)
    }
    \item{retainCache}{
        defaults to \code{FALSE}. If \code{FALSE} the data cache is deleted before
        returning. Setting this value to \code{TRUE} could save processing time
        on subsequent calls
    }
    \item{verbose}{
        if set to \code{TRUE} status information will be printed while the
        function executes
    }
    \item{probesetChunkSize}{
        determines how many SNPs will be processed at a time (in a single chunk).
        Keeping this number low reduces memory utilization and makes paralellism
        more fine-grained at the cost of managing more chunk objects. The
        default value provides a reasonable tradeoff
    }
}

\value{
    No return value. Instead PennCNV input files are written to \code{outdir}
}

\examples{
\dontrun{
# to see this example in context please read the vignette

load("../../../MouseDivData.RData")
dir.create("../../../lrr-baf-output")
genoVinoResult <- mouseDivGenotype(
    snpProbeInfo            = snpProbeInfo,
    snpInfo                 = snpInfo,
    referenceDistribution   = snpReferenceDistribution,
    chromosomes             = c(19, "X", "Y", "M"),
    celFiles                = expandCelFiles("../../../celFiles"),
    confScoreThreshold      = 0)
buildPennCNVInputFiles(
    outdir                          = "../../../lrr-baf-output",
    allowOverwrite                  = TRUE,
    genotypes                       = genoVinoResult$geno,
    snpProbeInfo                    = snpProbeInfo,
    snpInfo                         = snpInfo,
    snpReferenceDistribution        = snpReferenceDistribution,
    invariantProbeInfo              = invariantProbeInfo,
    invariantProbesetInfo           = invariantProbesetInfo,
    invariantReferenceDistribution  = invariantReferenceDistribution,
    celFiles                        = expandCelFiles("../../../celFiles"),
    isMale                          = genoVinoResult$isMale,
    chromosomes                     = c("19", "X", "Y", "M"))
  }
}

\keyword{CNV}
