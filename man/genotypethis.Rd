\author{Hyuna Yang}
\name{genotypethis}
\alias{genotypethis}

\title{Genotype the Mouse Diversity Array.}

\description{This function is for genotyping and vinotype the array. To genotype Mouse Diversity Array, \code{\link[MouseDivGeno]{MouseDivGenotype}} function would be handy, and for other arrays, obtain the log2 normalized intensities using other software, and genotype/vinotype the array using this function. }

\usage{  genotypethis(savefiledir, MM, SS, hint, isMale, trans, chr, pseudoautosomalID= c(21845:21870)) }

\arguments{
  \item{savefiledir}{Path to the directory in which output files will be saved. Default is the current directory.}
  \item{MM}{log2 normalized contrast intensities matrix. Each row corresponds to each SNP and each column corresponds to each strain, thus the size of MM should be number of SNP by number of strains. This can be log2(A)-log2(B) (MA transformation) or asinh(K*(A-B)/(A+B))/asinh(K) (CCS transformation). }
  \item{SS}{log2 normalized average intensities ( (log2(A) + log2(B))/2 ) matrix. The size should be same to MM matrix. }
  \item{hint}{A vector of mean intensity of AB group. Length of hint vector should be number of SNPs. User can use pre-calculated hint matrix, and it can be found at \tt{http:\\cgd.jax.org}. }
  \item{isMale}{TRUE (if a sample is male) or FALSE (if female) vector indicating gender status. The length of isMale vector should be corresponds to the sample size, and the order must follow column order of MM and SS matrixes. }
  \item{trans}{The MouseDivGenotype package genotype based on contrast (A allele intensity - B allele intensity) and average intensity dimensions rather than A allele and B allele dimension. To transform the data, it allows MA transformation, where M = A allele intensity - B allele intensity and A = average intensity, and CCS transformation, a slightly modified polar transformation used by BRLMM-P software.}
  \item{chr}{Chromosome in which MM and SS obtain. This is used to save the file. If MM or SS are obtained from many chromosomes just give any fake chromosome name. However, do not mix sex choromosome and autosome SNPs to make MM/SS matrix.}
  \item{pseudoautosomalID}{An index of psedoautosomal SNPs. Default is psedoautosomalID in Mouse Diversity Array, and it only affects to X chromosome. If you genotype other kind of array, this should be changed or marked as NULL. }
}
\value{
  It saves genodata* (genotype, 1=AA, 2=AB, 3=BB, -1 no call) vinodata* (vinotype, 1 = VINO, 0 or 2 otherwise) and confdata* (confidence score). 
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
  genotypethis(savefiledir, MM, SS, hint[[chr]], isMale, trans='CCStrans', chr, pseudoautosomalID= c(21845:21870))
  }
}

\keyword{genotype}
