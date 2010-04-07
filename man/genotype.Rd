\author{Hyuna Yang}
\name{genotype}
\alias{genotype}
\alias{genotypeSexchr}

\title{Genotype a SNP.}

\description{This function genotype the array. This function is used inside \code{MouseDivGeno} and \code{genotypethis}, and also can be used to genotype one SNP. genotypeSexchr is a function to genotype a SNP having AA and BB clusters only, such as SNPs from Y chromosomes (male samples only), male sample X chromosome, and Mitochondria.  }

\usage{
  genotype(nm, ns, hit, trans)
  genotypeSexchr(nm,ns,trans) 
}

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of contrast intensities of a SNP.}
  \item{hint}{A mean of AB group contrast intensities.}
  \item{trans}{CCStrans or MAtrans.}
}
\value{
  It returns genotype (1 = AA, 2 = AB, 3 = BB), VINOtype (1 = VINO, 0 = not), and confidence score.
}
\examples{
\dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  isMale = ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection, reference, trans="CCStrans", subset = TRUE)
  chr = 19
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, SS
  load('hint') # this is list
  genom = genotype(MM[1,], SS[1,], hint[[chr]][1], 'CCStrans')
  geno = genom$geno; vino = genom$vino; conf = genom$conf
  chr = 'M'
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, SS
  load('hint') # this is list
  genom = genotypeSexchr(MM[1,], SS[1,], 'CCStrans')
  geno = genom$geno; vino = genom$vino; conf = genom$conf
}
}

\keyword{genotype}
