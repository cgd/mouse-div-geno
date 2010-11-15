\author{Hyuna Yang and Keith Sheppard}
\name{genotypeAnyChrChunk}
\alias{genotypeAnyChrChunk}
\title{Platform Agnostic Genotyping}

\description{ 
  Unlike the mouseDivGenotype this function can be used to genotype SNPs data
  that did not originate from the mouse diversity array. If you are starting
  with A/B allele intensities you must first convert those intensities to
  contrast and average values. You can do this using a function like
  \code{convertToContrastAndAverag} or if you are starting
  with mouse diversity CEL files you can use \code{normalizeCelFileByChr}.
}

\usage{
genotypeAnyChrChunk(
        chr,
        intensityConts, intensityAvgs,
        hint = NULL,
        parIndices = NULL,
        transformMethod = c("CCStrans", "MAtrans"),
        isMale = NULL,
        confScoreThreshold = 1e-05)
}

\arguments{
    \item{chr}{a character string describing the chromosome being genotyped (this function expects that all SNPs being genotyped are on the same chromosome)}
    \item{intensityConts}{a matrix of (normalized) intensity contrasts. There must be a row per-SNP and a column per-sample. SNP ordering must be consistent with the ordering used for intensityAvgs but otherwise does not matter.}
    \item{intensityAvgs}{a matrix of (normalized) intensity averages. There must be a row per-SNP and a column per-sample. SNP ordering must be consistent with the ordering used for intensityConts but otherwise does not matter.}
    \item{hint}{an optional vector of hints for the (normalized) mean intensity of heterozygous SNPs}
    \item{parIndices}{this is only used in the case that chr == "X": a vector of indices indicating which SNPs fall into the pseudoautosomal region}
    \item{transformMethod}{indicates which transformation was used to generate intensityConts and intensityAvgs. Should be one of "CCStrans" or "MAtrans"}
    \item{isMale}{a vector of logicals indicating which samples are male. This is only used in the case that chr is "X" or "Y"}
    \item{confScoreThreshold}{the confidence threshold to use for genotyping. Any confidence scores below this threshold will be given a -1 rather than a valid genotype code}
}
\value{
    A list with "geno", "vino" and "conf" components
    \item{geno}{a matrix with a row per SNP and a column per sample. -1 = No call, 1 = AA, 2 = AB, and 3 = BB}
    \item{vino}{a matrix with a row per SNP and a column per sample. 1 indicates VINO, any other value indicates no-VINO}
    \item{conf}{a matrix with a row per SNP and a column per sample. These values are confidence scores based on Mohalanobis distance with chi-squre distribution approximation. Smaller confidence score implies less reliable data.}
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
}}

%\keyword{}
