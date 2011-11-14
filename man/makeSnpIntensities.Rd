\author{Hyuna Yang and Keith Sheppard}
\name{makeSnpIntensities}
\alias{makeSnpIntensities}
\title{Create a matrix-like object for storing SNP A and B intensities, SNP IDs and sample names}
\description{
    The object returned by this function has S3 classes \code{snpIntensities} and
    \code{probesetSampleMatrixes}. As such the documentation for
    \code{makeProbesetSampleMatrixes} applies to this object also (where the
    \code{mats} is set to \code{list(A=aIntensities, B=bIntensities)} and
    \code{probesetIds} is set to \code{snpIds})
}
\usage{makeSnpIntensities(aIntensities, bIntensities, snpIds, sampleNames)}

\arguments{
    \item{aIntensities}{
        the A intensities matrix. This matrix must have the same dimensions as
        \code{bIntensities}. Additionally the row count should match
        \code{length(snpIds)} and the column count should match
        \code{length(sampleNames)}.
    }

    \item{bIntensities}{
        the B intensities matrix. This matrix must have the same dimensions as
        \code{aIntensities}. Additionally the row count should match
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
}

\seealso{\code{\link{makeProbesetSampleMatrixes}}}

\value{
    a matrix-like S3 object with the class attribute set to
    \code{c("snpIntensities", "probesetSampleMatrixes", "list")}.
}

\examples{
# create the SNP intensities object from two matrixes
aInten <- matrix(
    c(
        4.7, 8.6, 3.1,
        3.3, 5.8, 9.0),
    ncol=3,
    byrow=TRUE)
bInten <- matrix(
    c(
        5.5, 3.3, 12.0,
        1.1, 3.1, 1.8),
    ncol=3,
    byrow=TRUE)
snpIntens <- makeSnpIntensities(
    aIntensities=aInten,
    bIntensities=bInten,
    snpIds=c("snp1", "snp2"),
    sampleNames=c("samp1", "samp2", "samp3"))

# prove that out snpIntens object works as expected.
# matrix-like functions and indexing on snpIntens
stopifnot(all(colnames(snpIntens) == c("samp1", "samp2", "samp3")))
stopifnot(all(rownames(snpIntens) == c("snp1", "snp2")))
stopifnot(nrow(snpIntens) == 2 && ncol(snpIntens) == 3)

# We can use matrix-like functions and indexing on snpIntens
# Note we can index by numbers, logical values or labels just like with matrices.
# Here we use label indexing to subset and reorder samples.
snpIntensSubset <- snpIntens[, c("samp3", "samp2")]
stopifnot(all(colnames(snpIntensSubset) == c("samp3", "samp2")))
stopifnot(all(rownames(snpIntensSubset) == c("snp1", "snp2")))
stopifnot(all(snpIntensSubset$A == matrix(c(3.1, 8.6, 9.0, 5.8), ncol=2, byrow=TRUE)))
stopifnot(all(snpIntensSubset$B == matrix(c(12.0, 3.3, 1.8, 3.1), ncol=2, byrow=TRUE)))
}

\keyword{class}
