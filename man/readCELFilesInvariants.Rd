\author{Hyuna Yang and Keith Sheppard}
\name{readCELFilesInvariants}
\alias{readCELFilesInvariants}
\title{Read Invariant Intensities From CEL Files}

\description{
    Reads invariant intensities from the given CEL files
}

\usage{
readCELFilesInvariants(
    celFiles,
    invariantProbeInfo,
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
        
        You can optionally use a list of data frames conforming to the
        description above rather than a single data frame. If you do so
        \code{referenceDistribution} must be either NULL or a list of numeric
        vectors. The length of the \code{referenceDistribution} list should also
        match the length of the \code{invariantProbeInfo} list
    }
    \item{referenceDistribution}{
        Reference distribution used for the quantile normalization. If there is
        an evidence that the density distribution of samples follows more than
        one distribution or different from that of classical laboratory strain,
        this normalization should not be used and leave as a default
        i.e. \code{referenceDistribution = NULL}.
        
        This parameter should either be a numeric vector or a list of numeric
        values depending on whether \code{invariantProbeInfo} is given as a simple
        data frame or a list of data frames.
    }
    \item{logFile}{
        the file name or connection that log messages should be written to. Use
        \code{stdout()} to print to standard-out. The default value of NULL
        indicates that no logging should be done
    }
    \item{numCores}{
        the number of cores to use. NULL indicates that we should use all
        available cores. Using this argument requires that the \pkg{multicore}
        package is installed.
    }
}
\value{
    This function returns a numeric matrix of intensities or a list of numeric
    matrices depending on whether invariantProbeInfo is given as a simple data
    frame or a list of data frames. The column names of the matrix(es) will
    correspond to sample IDs and the rows will correspond to probeset IDs
}

%\keyword{}
