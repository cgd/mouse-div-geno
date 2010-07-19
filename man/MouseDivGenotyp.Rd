\author{Hyuna Yang}
\name{MouseDivGenotype}
\alias{MouseDivGenotype}
\title{Genotype the Mouse Diversity Array}
\description{
  This function reads the .CEL files, normalizes the intensities, genotypes and vinotypes the Mouse Diversity Array. It also provides the input files used for the pennCNV (Wang, et al, 2007, a software to detect the copy number variance. Refer to CNV.R for detail information for detecting CNV.).  

  This function is highly tailored to Mouse Diversity Array, and if you want to genotype other types of arrays, refer to \code{\link[MouseDivGeno]{genotypethis}} function. All annotation files (allid, ABid, chrid) to run this function can be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. Coeffitients to normalized the array ({\tt CGFLcorrection} and {\tt reference}) or hint file used to genotype the array ({\tt hint}) were calculated based on 350 training .CEL files, and they can be obtained from the {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}.  

}
\usage{
MouseDivGenotype(allid, ABid, chrid, CGFLcorrection=NULL, reference=NULL, hint=NULL, trans=c("CCStrans","MAtrans"), celnamefile=NULL, mchr=c(1:19, 'X', 'Y', 'M'), celfiledir, outfiledir, subset=FALSE, confscore = 0.00001, doCNV = FALSE, exon1info=NULL, exon2info =NULL, exonoutfiledir, cnvoutfiledir)
}

\arguments{
  \item{allid}{SNP probes position at the .CEL file. } 
  \item{ABid}{A and B allele name and position. }
  \item{chrid}{chromosome and position of each probe set. }
  \item{CGFLcorrection}{Spline regression coefficients to correct the intensity variation due to different C or G contents in probe sequence and restriction fragment length. Default is not correcting intensity variation.}
  \item{reference}{Reference distribution used for the quantile normalization. If there is an evidence that the density distribution of samples follows more than one distribution or different from that of classical laboratory strain, this normalization should not be used and leave as a default i.e. {\tt reference = NULL}.}
 \item{hint}{Pre-calculated the mean of AB group in each probe set. }
  \item{trans}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
 \item{celnamefile}{A tab delimited .txt file containing the name of .CEL files. .CEL file name should be listed at the first column and the header MUST be 'celfile'. Optionally gender information can be provided by user, and if so it should be listed at the second column with header 'gender'. 'male' and 'female' (CASE SENSITIVE) only recognized and anything else (such as 'FEMALE', 'Male', 'Unknown' or 'NA') should be considered as unknown gender information. If 'filenames' file is not specified, the function reads all .CEL files under the celfiledir.}
  \item{mchr}{Index of chromosomes that you want to genotype. Chromosome 1 through 19 is indicated as number, X and Y chromosomes are 'X' and 'Y', and Mitochondria is indicated as 'M'. }
  \item{celfiledir}{Path to the directory in which all .CEL files are saved. Default is the current directory.}
  \item{outfiledir}{Path to the directory in which output will be saved. Default is the current directory.}
  \item{subset}{If this is TRUE, it only genotypes good probe sets. The definition of 'good' probe sets might be vary across laboratories, and current 'good' probe sets were identified by .CEL files processed at The Jackson Laboratory. }
  \item{confscore}{Threshold to assign the genotype based on confidence score. No call is assigned if the confidence score is less than {\tt confscore}. If {\tt confscore = 0}, it assigns genotype to all probe sets, and as the value gets big the assignment gets conservative. }
  \item{doCNV}{If this is true, it provides input files for the pennCNV. MouseDivGeno offers another option to study CNV, and refer to \code{\link[MouseDivGeno]{CNV.R}} function for detail information. }
  \item{exon1info}{Exon1 and exon2 are conventional name that Mouse Diversity Array used to distinguish types of exon. Exon1 and exon2 are invariant 25mer on exon, and if restriction fragment length is less than 1Kb, it is conventiallyed called an exon1 probes otherwise exon2 probes. Typically exon2 is dimmer than exon1, and exon1info and exon2info file provide position of probe sets, normalization coefficients, etc. If user wants to include exon1 for the pennCNV analysis, exon1info should be identified here. }
  \item{exon2info}{If user wants to include exon1 for the pennCNV analysis, exon1info should be identified here. }
  \item{exonoutfiledir}{Path to the directory in which exon intensities will be saved. Default is the current directory. This would not be used if doCNV = FALSE}
  \item{cnvoutfiledir}{Path to the directory in which pennCNV input files would be saved. Default is the current directory.}
}

\value{
 This function saves following information under the {\tt  outfiledir} directory. The name is rawdataMMchr* where * is chromosome index used at {\tt mchr}. The size of those files is $m X n$ where $m$ is number of SNPs and $n$ is number of arrays. 
 \item{rawdataMMchr}{Normalized log2 contrast intensities.}
 \item{rawdataSSchr}{Normalized log2 average intensities.}
 \item{genodata}{Genotype data. -1 = No call, 1 = AA, 2 = AB, and 3 = BB. Refer to anno files under the {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles} to get more information about each probe set.}
 \item{vinodata}{VINO data. 1 = VINO }
 \item{confdata}{Confidence score based on Mohalanobis distance with chi-squre distribution approximation. Smaller confidence score implies less reliable data.}  
 If {\tt doCNV = TRUE}, it also saves cnvChr* (LRR (log2(sample intensity/median intensity of each SNP in each group) BAF*) and cnv.pbf (population B allele frequency).  
}

\examples{
\dontrun{
  # you need to download all CDF and annotation files first. 
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference'); load('hint')
  load('exon1info'); load('exon2info')
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  celnamefile = '~/celname.txt'
  cnvoutfiledir = '~/cnvout'
  MouseDivGenotype(allid=allid, ABid=ABid, chrid=chrid, CGFLcorrection=CGFLcorrection, 
   reference=reference,  trans="CCStrans", hint = hint, celnamefile = celnamefile,
   mchr=c(17), celfiledir=celfiledir, subset=FALSE, confscore = 0.00001, doCNV = TRUE, 
   exon1info=exon1info, exon2info=exon2info, cnvoutfiledir=cnvoutfiledir)
  }
}

\keyword{genotype}
