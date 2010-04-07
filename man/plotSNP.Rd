\author{Hyuna Yang}
\name{plotSNP}
\alias{plotSNP}

\title{plot SNP}

\description{
  plot genotype/vinotype of a SNP 
}

\usage{
plotSN(nm, ns, geno, vino, main='', xlab='Contrast', ylab = 'Average Intensity')
}

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of genotype of a SNP.}
}

\examples{
\dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  dir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  isMale = ReadCel(dir, outfiledir, allid, ABid, chrid, CGFLcorrection, reference, 
    trans="CCStrans")
  chr = 19
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, SS
  load('hint') # this is list
  genom = genotype(MM[1,], SS[1,], hint[[chr]][1], 'CCStrans')
  geno = genom$geno; vino = genom$vino; conf = genom$conf
  plotSNP(MM[1,], SS[1,], geno, vino)
}
}

\keyword{plot}
