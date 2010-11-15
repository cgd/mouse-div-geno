\author{Hyuna Yang and Keith Sheppard}
\name{convertToContrastAndAverage}
\alias{convertToContrastAndAverage}
\title{Convert A/B Intensities to Contrast and Average}

\description{
	Given corresponding A/B allele intensities this function will calculate
	intensity dQuote{average} and dQuote{contrast} using one of two transforms. 
}

\usage{
convertToContrastAndAverage(
        aIntensities,
        bIntensities,
        transformMethod = c("CCStrans", "MAtrans"),
        intensitiesAreLog2 = FALSE,
        k = 4)
}

\arguments{
    \item{aIntensities}{intensities for the A alleles. This argument can be a vector or a matrix but its shape must match up with bIntensities}
    \item{bIntensities}{intensities for the B alleles. This argument can be a vector or a matrix but its shape must match up with aIntensities}
    \item{transformMethod}{This package genotypes SNPs based on intensity contrasts and averages, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
    \item{intensitiesAreLog2}{should be set to TRUE if the given intensities are log2 normalized}
    \item{k}{this is the dQuote{k} parameter used for CCS transform}
}
\value{
    A list with dQuote{intensityConts} and dQuote{intensityAvgs} components.
    These components will have the same shape as the A and B intensities which
    were passed in.
}

\examples{
\dontrun{
# to see this example in context please read the vignette

# loads aIntensities and bIntensities objects which contain log2 normalized
# A/B SNP allele intensities. Matrix rows are SNPs and matrix columns are
# are samples
load("../../../abIntensities.RData")

# for this example I will only genotype SNPs 15 and 1281
intenCont <- convertToContrastAndAverage(
        aIntensities[c(15, 1281), ],
        bIntensities[c(15, 1281), ],
        intensitiesAreLog2 = TRUE)
}}

%\keyword{}
