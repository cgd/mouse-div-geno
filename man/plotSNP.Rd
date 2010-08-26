\author{Hyuna Yang}
\name{plotSNP}
\alias{plotSNP}

\title{plotSNP}

\description{
  Plot genotype/vinotype of a SNP.
}

\usage{
plotSNP(nm, ns, geno, vino, main = "", xlab = "Contrast", ylab = "Average Intensity")
}

\arguments{
  \item{nm}{A vector of contrast intensities of a SNP.}
  \item{ns}{A vector of average intensities of a SNP.}
  \item{geno}{A vector of genotype of a SNP.}
  \item{vino}{A vector of vinotype of a SNP.}
  \item{main}{the plot title}
  \item{xlab}{the x axis label}
  \item{ylab}{the y axis label}
}

\value{
 This function returns a scatter plot of one probe set. It uses contrast (x axes) and average intensities (y axes), and each color indicates red = AA, green = AB, blue = BB, purple = VINO, and grey = no call.
}

\examples{
\dontrun{
  #== there are two ways to get genotype and vinotype
  #=== option 1 : get genotype using genotype function
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  isMale = ReadCel(celfiledir=celfiledir, outfiledir=outfiledir, allid, ABid, chrid, 
     CGFLcorrection=CGFLcorrection, reference=reference, trans="CCStrans")
  chr = 19
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 contrast intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 average intensity, the output of ReadCel, SS
  load('hint') # this is list
  genom = genotype(MM[1,], SS[1,], hint[[chr]][1], 'CCStrans')
  geno = genom$geno; vino = genom$vino;
  plotSNP(MM[1,], SS[1,], geno, vino)

  #=== option 2 : get genotype using MouseDivGenotype function
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference'); load('hint')
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  celnamefile = '~/celname.txt'
  chr = 19
  MouseDivGenotype = function(allid, ABid, chrid, CGFLcorrection=CGFLcorrection, reference=reference,
    hint=hint, trans='CCStrans', celnamefile=celnamefile, mchr=chr, celfiledir=celfiledir, 
    outfiledir=outfiledir, subset=TRUE, confscore = 0.00001)
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 contrast intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 average intensity, the output of ReadCel, SS
  xname = paste(outfiledir, '/genodata' chr, sep='', collapse='')
  load(xname) # this will read chr19 genotype data
  xname = paste(outfiledir, '/vinodata' chr, sep='', collapse='')
  load(xname) # this will read chr19 vinotype data
  plotSNP(MM[1,], SS[1,], geno[1,], vino[1,])
  }
}

\keyword{plot}
