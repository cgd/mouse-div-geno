\author{Hyuna Yang}
\name{plotSNP}
\alias{plotSNP}

\title{plotSNP}

\description{
  Plot genotype/vinotype of a SNP.
}

\usage{
plotSNP(nm, ns, geno, vino, main = "", xlab = "Contrast", ylab = "Average Intensity")
}

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of average intensities of a SNP.}
  \item{geno}{A vector of genotype of a SNP.}
  \item{vino}{A vector of vinotype of a SNP.}
  \item{main}{the plot title}
  \item{xlab}{the x axis label}
  \item{ylab}{the y axis label}
}

\value{
 This function returns a scatter plot of one probe set. It uses contrast (x axes) and average intensities (y axes), and each color indicates red = AA, green = AB, blue = BB, purple = VINO, and grey = no call.
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
genoVinoConfForChr19 <- genotypeAnyChrChunk(
        "19",
        intenCont$intensityConts,
        intenCont$intensityAvgs)
plotSNP(intenCont$intensityConts[1, ],
        intenCont$intensityAvgs[1, ],
        genoVinoConfForChr19$geno[1, ],
        genoVinoConfForChr19$vino[1, ])
  }
}

\keyword{plot}
