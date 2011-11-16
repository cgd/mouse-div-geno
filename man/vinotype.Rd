\author{Hyuna Yang and Keith Sheppard}
\name{vinotype}
\alias{vinotype}

\title{Vinotype a SNP.}

\description{This function will not typically be used directly by users. This
function identifies SNPs having low intensity, indicating a failure of
hybridization in both allele, often caused by a new mutation somewhere in probe
sequence or deletion (when it occurs consecutively in the genome).  We call them
Variable INtensity Oligonucleotide (VINO) and this function is to find those
SNPs. This function is called by \code{genotype(...)}.}

\usage{vinotype(nm, ns, geno, logfn=NULL)}

\arguments{
    \item{nm}{A vector of contrast intensities of a SNP.}
    \item{ns}{A vector of average intensities of a SNP.}
    \item{geno}{A vector of genotype information of the SNP (1=AA, 2=AB, 3=BB).}
    \item{logfn}{
        A function used to write a line to our debug log. If non-NULL, the
        arguments accepted by this function must work the same way as the
        \code{sprintf(...)} function
    }
}
\value{
  It returns vinotype (1 = VINO, 0 = not), and confidence score
}

\keyword{vinotype}
