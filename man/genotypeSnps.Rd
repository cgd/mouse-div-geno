\author{Hyuna Yang and Keith Sheppard}
\name{genotypeSnps}
\alias{genotypeSnps}
\title{Genotype SNP Data Matrices}
\description{
    This is a simple general-purpose function for genotyping SNP intensity matrices
    If you are working with CEL or tab delimited files it's more simple and efficient
    to use \code{mouseDivGenotypeCEL} or \code{mouseDivGenotypeTab}.
}

\usage{
genotypeSnps(
        snpContAvg,
        snpInfo,
        isMale = NULL,
        samplesToKeep = NULL,
        chromosomes = c(1:19, "X", "Y", "M"),
        confScoreThreshold = 1e-05,
        numCores = NULL,
        probesetChunkSize = 3000,
        logFile = NULL)
}

\arguments{
    \item{snpContAvg}{
        an object of class \code{snpContAvg}. See documentation for
        \code{makeSnpContAvg(...)} for details about this class
    }
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{snpId}: the ID of this SNP
        
        \code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}
        
        \code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygous SNPs
        
        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is
        in the pseudoautosomal region (PAR). This should only be set to
        \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
    }
    \item{isMale}{
        a logical vector indicating which samples are male (\code{TRUE}) and
        which are female (\code{FALSE}). If \code{is.null(isMale)} then
        gender will be inferred from \code{snpContAvg}. Otherwise
        \code{length(isMale) == ncol(snpContAvg)} must be true.
    }
    \item{samplesToKeep}{
        this argument determines which samples we will generate output data for.
        The default value of NULL indicates that we should generate output
        for all available samples. The reason you may want to use this parameter
        is that the genotyping algorithm works best when it is run with many samples.
        So, you may want to genotype with many samples but you may only be interested
        in the genotype values of some subset of these samples
        This value can either be an integer vector
        which (will identify samples by index) or a character vector (which
        will identify samples by sample ID).
    }
    \item{chromosomes}{
        character vector of chromosome IDs that you want to genotype. These
        should correspond to the values in \code{snpInfo$chrId}. Eg:
        \code{chromosomes=c(1:19, "X", "Y", "M")}.
    }
    \item{confScoreThreshold}{
        if confidence scores fall below this threshold we assign a genotype
        value of -1
    }
    \item{numCores}{
        the number of cores to use. NULL indicates that we should use all
        available cores. This argument will have no affect if the \pkg{multicore}
        package is not installed.
    }
    \item{probesetChunkSize}{
        determines how many SNPs will be processed at a time (in a single chunk).
        Keeping this number low reduces memory utilization and makes parallelism
        more fine-grained at the cost of managing more chunk objects. The
        default value provides a reasonable tradeoff
    }
    \item{logFile}{
        the file name or connection that log messages should be written to.
        Use stdout() to print to standard-out
    }
}

\value{
    A list containing all of the following components:
    
    \item{geno}{
        a matrix with a row per SNP and a column per sample.
        -1 = No call, 1 = AA, 2 = AB, and 3 = BB}
    \item{vino}{
        a matrix with a row per SNP and a column per sample.
        1 indicates VINO, 0 indicates no-VINO}
    \item{conf}{
        a matrix with a row per SNP and a column per sample.
        These values are confidence scores based on Mohalanobis distance with
        chi-squre distribution approximation. Smaller confidence score implies
        less reliable data.}
    \item{isMale}{
        A per-CEL-file vector of logicals indicating which were genotyped as
        male. If gender was supplied via the CEL files parameter then this value
        will match that otherwise it will be the inferred gender values.}
}

\examples{
\dontrun{
# Using genotypeSnps on CEL files (if you have many CEL files it is much better
# to use mouseDivGenotypeCEL instead of genotypeSnps)
celData <- readCELFiles("./celFiles/", snpProbeInfo, snpReferenceDistribution, logFile=stdout())
transformedData <- ccsTransform(celData)
genoResults <- genotypeSnps(transformedData, snpInfo, chromosomes=c("19", "X", "Y", "M"), logFile=stdout())
}

#
}

\keyword{genotype}
