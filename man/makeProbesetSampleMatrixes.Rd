\author{Hyuna Yang and Keith Sheppard}
\name{makeProbesetSampleMatrixes}
\alias{makeProbesetSampleMatrixes}
\title{Create a matrix-like object from data matrices, probeset IDs and sample names}
\description{
    The object returned by this function has S3 class \code{probesetSampleMatrixes}.
    The main idea of this class is that it holds a list of matrices all of which
    must have the same dimensions. The row count must match the length of
    \code{probesetIds} and the column count must match the length of
    \code{sampleNames} thus the matrix rows will correspond to probesets and
    columns will correspond to samples. There is no restriction on the number
    of matrices or the type of data held by these matrices.
    
    This class is a generalization of the \code{snpContrasts} and
    \code{snpIntensities} classes, both of which use \code{probesetSampleMatrixes}
    as a base class and specify further restrictions on what kind of matrices are
    stored.
}
\usage{makeProbesetSampleMatrixes(mats, probesetIds, sampleNames)}

\arguments{
    \item{mats}{
        the list of data matrices to include. For all matrices \code{m} in this list the
        following must hold true:
        \code{nrow(m) == length(probesetIds) && ncol(m) == length(sampleNames)}
    }

    \item{probesetIds}{
        this parameter must be a character vector. The given strings will serve
        as row labels for all of the given matrices
    }

    \item{sampleNames}{
        this parameter must be a character vector. The given strings will serve
        as row labels for all of the given matrices
    }
}

\seealso{
    \code{\link{makeSnpContAvg}}
    \code{\link{makeSnpIntensities}}
}

\value{
    a matrix-like S3 object with the class attribute set to
    \code{c("probesetSampleMatrixes", class(mats))}.
}

\examples{
# in this section we'll create the probesetSampleMatrixes
avgIntensities <- matrix(
    c(
        5.5, 3.3, 12.0,
        4.7, 8.6, 3.1,
        3.3, 5.8, 9.0,
        1.1, 3.1, 1.8),
    ncol=3,
    byrow=TRUE)
flags <- matrix(
    c(
        TRUE, FALSE, TRUE,
        FALSE, FALSE, TRUE,
        FALSE, FALSE, TRUE,
        TRUE, TRUE, TRUE),
    ncol=3,
    byrow=TRUE)
sampMats <- makeProbesetSampleMatrixes(
    mats=list(avgIntensities=avgIntensities, flags=flags),
    probesetIds=c("probe1", "probe2", "probe3", "probe4"),
    sampleNames=c("samp1", "samp2", "samp3"))

# prove that out sampMats object works as expected.
stopifnot(all(colnames(sampMats) == c("samp1", "samp2", "samp3")))
stopifnot(all(rownames(sampMats) == c("probe1", "probe2", "probe3", "probe4")))
stopifnot(nrow(sampMats) == 4 && ncol(sampMats) == 3)

# We can use matrix-like functions and indexing on sampMats
# Note we can index by numbers, logical values or labels just like with matrices.
# Here we just use numeric indexing.
sampMatsSubset <- sampMats[1:2, 1:2]
stopifnot(nrow(sampMatsSubset) == 2 && ncol(sampMatsSubset) == 2)
stopifnot(all(colnames(sampMatsSubset) == c("samp1", "samp2")))
stopifnot(all(rownames(sampMatsSubset) == c("probe1", "probe2")))
stopifnot(all(sampMatsSubset$avgIntensities == matrix(c(5.5, 3.3, 4.7, 8.6), ncol=2, byrow=TRUE)))
stopifnot(all(sampMatsSubset$flags == matrix(c(TRUE, FALSE, FALSE, FALSE), ncol=2, byrow=TRUE)))
}

\keyword{class}
