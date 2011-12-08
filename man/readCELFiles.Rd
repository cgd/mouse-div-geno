\author{Hyuna Yang and Keith Sheppard}
\name{readCELFiles}
\alias{readCELFiles}
\title{Read CEL Files}

\description{
    Reads SNP intensities from the given CEL files
}

\usage{
readCELFiles(
        celFiles,
        snpProbeInfo,
        referenceDistribution = NULL,
        logFile = NULL,
        numCores = NULL)
}

\arguments{
    \item{celFiles}{
        The CEL files to read. Must be a vector of CEL file name(s) or
        directory name(s). Any directories will be expanded to the listing of
        CEL files found in the directory
    }
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for
        B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs
        should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied
        by adding it to probe mean intesity vales
    }
    \item{referenceDistribution}{
        Reference distribution used for the quantile normalization. If there is
        an evidence that the density distribution of samples follows more than
        one distribution or different from that of classical laboratory strain,
        this normalization should not be used and leave as a default
        i.e. \code{referenceDistribution = NULL}.
    }
    \item{logFile}{
        if non-NULL this must be a connection type to which log messages will
        be written
    }
    \item{numCores}{
        the number of cores to use. NULL indicates that we should use all
        available cores. Using this argument requires that the \pkg{multicore}
        package is installed.
    }
}
\value{
    This function returns an object of class \code{snpIntensities}. See the
    documentation for the \code{makeSnpIntensities(...)} function for a detailed
    description of this class.
}

%\keyword{}
