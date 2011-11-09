\author{Hyuna Yang and Keith Sheppard}
\name{convertToContrastAndAverage}
\alias{maTransform}
\alias{ccsTransform}
\title{Convert A/B Intensities to Contrast and Average}

\description{
	Given corresponding A/B allele intensities these functions will calculate
	intensity dQuote{average} and dQuote{contrast} using one of two transforms. 
}

\usage{
maTransform(ab)
ccsTransform(ab, k = 4)
}

\arguments{
    \item{ab}{intensities for the A and B alleles. This should have an S3 class of snpIntensities which you can create using makeSnpIntensities}
    \item{k}{this is the dQuote{k} hyperparameter used for CCS transform}
}
\value{
    CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm
    obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is
    intensity of A and B allele and K is a hyperparameter. MA transformation
    obtains contrast by log2(A) - log2(B). In both cases summation is defined
    by (log2(A) + log2(B))/2

    The return object will be of S3 class snpContAvg. To learn more about the
    snpContAvg class see documentation for the makeSnpContAvg function.
}
