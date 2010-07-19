\author{Hyuna Yang}
\name{vinotype}
\alias{vinotype}

\title{Vinotype a SNP.}

\description{This function identifies SNPs having low intensity, indicating a failure of hybridization in both allele, often caused by a new mutation somewhere in probe sequence or deletion (when it occurs consecutively in the genome).  We call them Variable INtensity Oligonucleotide (VINO) and this function is to find those SNPs. This function is used inside  \code{MouseDivGeno}, \code{genotypethis} and \code{genotype}. It also can be used to vinotype one SNP, but since it needs genotype information it is hardly used by itself.}

\usage{  vinotype(nm, ns, geno, doCNV = FALSE) }

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of average intensities of a SNP.}
  \item{geno}{A vector of genotype information of the SNP (1=AA, 2=AB, 3=BB).}
  \item{doCNV}{True or false indicator to generate CNV analysis input files.}
}
\value{
  It returns vinotype (1 = VINO, 0 = not), confidence score, and two CNV analysis related outputs (llr = log2(average intensity/median intensity) and BAF(B allele frequency)) if doCNV is TRUE.
}
\examples{
\dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  clefiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  isMale = ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection=NULL, 
    reference=NULL, trans="CCStrans",  subset = FALSE)
  chr = 19
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 contrast intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 average intensity, the output of ReadCel, SS
  load('hint') # this is list
  genom = genotype(MM[1,], SS[1,], hint[[chr]][1], 'CCStrans')
  geno = genom$geno; vino = genom$vino; conf = genom$conf
  # or if you see what vinotype does
  vinom= vinotype(MM[1,], SS[1,], geno) 
  vino1 = vinom$vino; conf1 = vinom$conf 
  # {\tt vino1} and {\tt conf1} are exactly same as {\tt vino} and {\tt conf}.
}
}

\keyword{vinotype}
