\author{Hyuna Yang}
\name{vinotype}
\alias{vinotype}

\title{Vinotype a SNP.}

\description{This function identifies SNPs having low intensity, indicating a failure of hybridization in both allele, often caused by a new mutation somewhere in probe sequence or deletion (when it occurs consecutively in the genome).  We call them Variable INtensity Oligonucleotide (VINO) and this function is to find those SNPs. This function is used inside  \code{MouseDivGeno}, \code{genotypethis} and \code{genotype}. It also can be used to vinotype one SNP, but since it needs genotype information it is hardly used by itself.}

\usage{  vinotype(nm, ns, geno) }

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of average intensities of a SNP.}
  \item{geno}{A vector of genotype information of the SNP (1=AA, 2=AB, 3=BB).}
}
\value{
  It returns vinotype (1 = VINO, 0 = not), and confidence score
}
\examples{
\dontrun{
    # TODO add example code
}
}

\keyword{vinotype}
