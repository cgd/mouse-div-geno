\author{Hyuna Yang and Keith Sheppard}
\name{genotypeAnyChrChunk}
\alias{genotypeAnyChrChunk}
\title{Platform Agnostic Genotyping}

\description{ 
  Unlike the mouseDivGenotype this function can be used to genotype SNPs data
  that did not originate from the mouse diversity array.
}

\usage{
genotypeAnyChrChunk(chr, ms, ss, hint = NULL, parIndices = NULL, trans, isMale = NULL, confScoreThreshold = 1e-05)
}

\arguments{
    \item{chr}{a character string describing the chromosome being genotyped (this function expects that all SNPs being genotyped are on the same chromosome)}
    \item{ms}{a matrix of (normalized) intensity contrasts. There must be a row per-SNP and a column per-sample. SNP ordering must be consistent with the ordering used for ss but otherwise does not matter.}
    \item{ss}{a matrix of (normalized) intensity averages. There must be a row per-SNP and a column per-sample. SNP ordering must be consistent with the ordering used for ms but otherwise does not matter.}
    \item{hint}{an optional vector of hints for the (normalized) mean intensity of heterozygous SNPs}
    \item{parIndices}{this is only used in the case that chr == "X": a vector of indices indicating which SNPs fall into the pseudoautosomal region}
    \item{trans}{indicates which transformation was used to generate ms, and ss. Should be one of "CCStrans" or "MAtrans"}
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
#TODO add usage example
}}

\keyword{}
