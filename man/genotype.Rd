\author{Hyuna Yang and Keith Sheppard}
\name{genotype}
\alias{genotype}
\alias{genotypeHomozygous}

\title{Genotype a SNP.}

\description{
If you want to directly genotype intensity contrasts and averages it is
recommended that you use the \code{genotypeAnyChrChunk(...)} function
rather than either of these, since it will call these functions for you.
These functions genotype several arrays at a single SNP. genotypeHomozygous
can be used when you know ahead of time that your samples can only have
homozygous genotypes such as SNPs from male Y and X chromosomes, and Mitochondria
}

\usage{
  genotype(nm, ns, hint1, trans, logfn=NULL)
  genotypeHomozygous(nm, ns, trans, logfn=NULL)
}

\arguments{
    \item{nm}{A vector of intensity contrasts for a single SNP.}
    \item{ns}{A vector of intensity averages for a single SNP.}
    \item{hint1}{A mean of AB group contrast intensities.}
    \item{trans}{CCStrans or MAtrans.}
    \item{logfn}{
        A function used to write a line to our debug log. If non-NULL, the
        arguments accepted by this function must work the same way as the
        \code{sprintf(...)} function
    }
}
\value{
  It returns genotype (1 = AA, 2 = AB, 3 = BB), VINOtype (1 = VINO, 0 = not), and confidence score
}

\keyword{genotype}
