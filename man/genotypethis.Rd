\author{Hyuna Yang}
\name{genotypethis}
\alias{genotypethis}

\title{Genotype the Mouse Diversity Array.}

\description{This function is to genotype and vinotype the array. To genotype Mouse Diversity Array, \code{\link[MouseDivGeno]{MouseDivGenotype}} function would be handy, and for other types of arrays obtain the log2 normalized intensities using other software and do genotype/vinotype using this function. }

\usage{  genotypethis(savefiledir, MM, SS, hint, isMale, trans, chr, confscore, doCNV) }

\arguments{
  \item{savefiledir}{Path to the directory in which output files will be saved. Default is the current directory.}
  \item{MM}{log2 normalized contrast intensities matrix. Each row corresponds to each SNP and each column corresponds to each strain, and thus the size of MM should be number of SNP by number of arrays. This can be log2(A)-log2(B) (MA transformation) or asinh(K*(A-B)/(A+B))/asinh(K) (CCS transformation). }
  \item{SS}{log2 normalized average intensities ( (log2(A) + log2(B))/2 ) matrix. The size should be same to the {\tt MM} matrix. }
  \item{hint}{A vector of mean of contrast intensities in AB group. Length of hint vector should be number of SNPs.}
  \item{isMale}{TRUE (if a sample is male) or FALSE (if female) vector indicating gender status. The length of isMale vector should be corresponds to the sample size, and the order must follow column order of MM and SS matrixes. }
  \item{trans}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
  \item{chr}{Chromosome index from which MM and SS obtain (1-19, X, Y and M (mitochondria)). This is used to save the file and distinguish Autosomes from sex chromosomes. When MM or SS are obtained from many chromosomes just give any fake chromosome index from 1 to 19. However, do not mix sex chromosome SNPs with Autosomes SNPs to make MM/SS matrix.}
  \item{confscore}{Threshold to assign the genotype based on confidence score. No call is assigned if the confidence score is less than {\tt confscore}. If {\tt confscore = 0}, it assigns genotype to all probe sets, and as the value gets big the assignment gets conservative. }
  \item{doCNV}{If this is true, it provides input files for the pennCNV. Refer to \code{\link[MouseDivGeno]{CNV.R}} function for detail information. The default is FALSE; no CNV analysis.}
}
\value{
  It saves genodata* (genotype, 1=AA, 2=AB, 3=BB, -1 no call) vinodata* (vinotype, 1 = VINO) and confdata* (confidence score), where input {\tt chr} is used in *. If {\tt doCNV = TRUE}, it also saves cnvChr* (LRR (log2(sample intensity/median intensity of each SNP in each group) BAF*) and cnv.pbf (population B allele frequency).  
}
\examples{
  \dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  isMale = ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection, reference, 
    trans="CCStrans")
  chr = 19
  xname = paste(outfiledir, '/rawdataMMchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, MM 
  xname = paste(outfiledir, '/rawdataSSchr' chr, sep='', collapse='')
  load(xname) # this will read chr19 raw intensity, the output of ReadCel, SS
  load('hint') # this is list
  genotypethis(savefiledir = outfiledir, MM = MM, SS = SS, hint = hint[[chr]], isMale = isMale, trans='CCStrans', chr = chr, confscore = 0.0001)
  # or without hint file and use current directory to save output files
  genotypethis(MM = MM, SS = SS, hint = NULL, isMale = isMale, trans='CCStrans', chr = chr)
  }
}

\keyword{genotype}
