\author{Hyuna Yang and Keith Sheppard}
\name{convertToContrastAndAverag}
\alias{convertToContrastAndAverag}
\title{Convert A/B Intensities to Contrast and Average}

\description{
	Given corresponding A/B allele intensities this function will calculate
	intensity dQuote{average} and dQuote{contrast} using one of two transforms. 
}

\usage{
convertToContrastAndAverag(
        aIntensities,
        bIntensities,
        transformMethod = c("CCStrans", "MAtrans"),
        intensitiesAreLog2 = FALSE,
        k = 4)
}

\arguments{
    \item{aIntensities}{intensities for the A alleles. This argument can be a vector or a matrix but its shape must match up with bIntensities}
    \item{bIntensities}{intensities for the B alleles. This argument can be a vector or a matrix but its shape must match up with aIntensities}
    \item{transformMethod}{determines which transformation method should be used}
    \item{intensitiesAreLog2}{should be set to TRUE if the given intensities are log2 normalized}
    \item{k}{this is the dQuote{k} parameter used for CCS transform}
}
\value{
    A list with "geno", "vino" and "conf" components
    \item{geno}{a matrix with a row per SNP and a column per sample. -1 = No call, 1 = AA, 2 = AB, and 3 = BB}
    \item{vino}{a matrix with a row per SNP and a column per sample. 1 indicates VINO, any other value indicates no-VINO}
    \item{conf}{a matrix with a row per SNP and a column per sample. These values are confidence scores based on Mohalanobis distance with chi-squre distribution approximation. Smaller confidence score implies less reliable data.}
}

\examples{
\dontrun{
#TODO add usage example
}}

%\keyword{}
