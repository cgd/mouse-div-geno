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
    \item{transformMethod}{determines which transformation method should be used}
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
#TODO add usage example
}}

%\keyword{}
