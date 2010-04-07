\author{Hyuna Yang}
\name{vinotype}
\alias{vinotype}

\title{Vinotype a SNP.}

\description{This function identifies SNPs having low intensity, indicating a failure of hybridization in both allele, often cased by a new mutation somewhere in probe sequence or deletion (when it occurs consecutively in the genome).  We call them Variable INtensity Oligonucleotide (VINO) and this function is to find those SNPs. This function is used inside  \code{MouseDivGeno}, \code{genotypethis} and \code{genotype}. It also can be used to vinotype one SNP, but since it needs genotype information and \code{genotype} uses \code{vinotype} inside the code, it is hardly used by itself.}

\usage{  vinotype(nm, ns, geno) }

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of contrast intensities of a SNP.}
  \item{geno}{A vector of genotype information of the SNP (1=AA, 2=AB, 3=BB).}
}
\value{
  It returns vinotype (1 = VINO, 0 = not) and confidence score. 
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
  vinom= vinotype(MM[1,], SS[1,], geno) 
  vino1 = vinom$vino; conf1 = vinom$conf # vino1 and conf1 are exactly same as vino and conf
}
}

\keyword{vinotype}
