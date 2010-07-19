\author{Hyuna Yang}
\name{genotypeAnyChrChunk}
\alias{genotypeAnyChrChunk}

\title{Genotype the Mouse Diversity Array.}

\description{This function is to genotype and vinotype a "chunk" of data from an array. If you are working with the Mouse Diversity Array you should most likely use \code{\link[MouseDivGeno]{MouseDivGenotype}} to perform your genotyping, however this function may provide useful entry point for you to extend this package's genotyping and vinotyping capabilities for other array types.}

\usage{
genotypeAnyChrChunk <- function(chr, ms, ss, hint, parIndices, trans, isMale)
}

\arguments{
  \item{chr}{Chromosome index from which MM and SS obtain (1-19, X, Y and M (mitochondria)). This is used to save the file and distinguish Autosomes from sex chromosomes. When MM or SS are obtained from many chromosomes just give any fake chromosome index from 1 to 19. However, do not mix sex chromosome SNPs with Autosomes SNPs to make MM/SS matrix.}
  \item{mm}{log2 normalized contrast intensities matrix. Each row corresponds to each SNP and each column corresponds to each strain, and thus the size of MM should be number of SNP by number of arrays. This can be log2(A)-log2(B) (MA transformation) or asinh(K*(A-B)/(A+B))/asinh(K) (CCS transformation). }
  \item{ss}{log2 normalized average intensities ( (log2(A) + log2(B))/2 ) matrix. The size should be same to the \code{MM} matrix. }
  \item{hint}{A vector of mean of contrast intensities in AB group. Length of hint vector should be number of SNPs.}
  \item{parIndices}{a (possibly empty) vector of indices representing all mm and ss rows that are in the pseudoautosomal region. If this vector is non-NULL then chr should be "X"}
  \item{trans}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
  \item{isMale}{TRUE (if a sample is male) or FALSE (if female) vector indicating gender status. The length of isMale vector should be corresponds to the sample size, and the order must follow column order of MM and SS matrixes. }
}
\value{
    This function returns a list with three items named:
    \item{geno}{a matrix with a row per SNP and a column per CEL file. -1 = No call, 1 = AA, 2 = AB, and 3 = BB}
    \item{vino}{a matrix with a row per SNP and a column per CEL file. 1 indicates VINO, any other value indicates no-VINO}
    \item{conf}{a matrix with a row per SNP and a column per CEL file. These values are confidence scores based on Mohalanobis distance with chi-squre distribution approximation. Smaller confidence score implies less reliable data.}
}

\keyword{genotype}
