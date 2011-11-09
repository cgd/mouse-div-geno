\author{Keith Sheppard}
\name{makeSnpContAvg}
\alias{makeSnpContAvg}
\title{Create a matrix-like object for storing SNP contrasts, SNP averages, SNP IDs and sample names}
\description{
    The object returned by this function has S3 classes \code{snpContAvg} and
    \code{probesetSampleMatrixes}. As such the documentation for
    \code{makeProbesetSampleMatrixes} applies to this object also (where the
    \code{mats} is set to \code{list(snpContrasts=snpContrasts, snpAverages=snpAverages)} and
    \code{probesetIds} is set to \code{snpIds})
}
\usage{makeSnpContAvg(snpContrasts, snpAverages, snpIds, sampleNames, transformMethod=NULL)}

\arguments{
    \item{snpContrasts}{
        the SNP allele intensity contrast matrix. The row count of this matrix should match
        \code{length(snpIds)} and the column count should match
        \code{length(sampleNames)}.
    }

    \item{snpAverages}{
        the SNP A,B intensity averages. The row count of this matrix should match
        \code{length(snpIds)} and the column count should match
        \code{length(sampleNames)}.
    }

    \item{snpIds}{
        this parameter must be a character vector. The given strings will serve
        as row labels for the matrix-like object created
    }

    \item{sampleNames}{
        this parameter must be a character vector. The given strings will serve
        as row labels for the matrix-like object created
    }
    
    \item{transformMethod}{
        this argument can be used to indicate which transformation was used
        to convert from A and B SNP intensities to average intensity and
        contrast. "MA" and "CCS" are the currently supported methods.
    }
}

\seealso{
\code{\link{makeProbesetSampleMatrixes}}
\code{\link{ccsTransform}}
\code{\link{maTransform}}
}

\value{
    a matrix-like S3 object with the class attribute set to
    \code{c("snpContAvg", "probesetSampleMatrixes", "list")}.
}

\keyword{class}
