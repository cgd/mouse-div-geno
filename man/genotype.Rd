\author{Hyuna Yang}
\name{genotype}
\alias{genotype}
\alias{genotypeHomozygous}

\title{Genotype a SNP.}

\description{This function genotype the array. This function is used inside \code{MouseDivGeno} and \code{genotypethis}, and also can be used to genotype one SNP. genotypeSexchr is a function to genotype a SNP having AA and BB clusters only, such as SNPs from male Y and X chromosomes, and Mitochondria.  }

\usage{
  genotype(nm, ns, hint1, trans)
  genotypeHomozygous(nm, ns, trans) 
}

\arguments{
  \item{nm}{A vector of intensity contrasts for a single SNP.}
  \item{ns}{A vector of intensity averages for a single SNP.}
  \item{hint1}{A mean of AB group contrast intensities.}
  \item{trans}{CCStrans or MAtrans.}
}
\value{
  It returns genotype (1 = AA, 2 = AB, 3 = BB), VINOtype (1 = VINO, 0 = not), and confidence score
}
\examples{
\dontrun{
    # TODO add examples
}
}

\keyword{genotype}
