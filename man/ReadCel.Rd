\author{Hyuna Yang}
\name{ReadCelFile}
\alias{ReadCelFile}
\title{Normalized the Mouse Diversity Array}
\description{
  This function reads the .CEL files, normalizes the intensities, and saves normalized intensities.}
\usage{
ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection=NULL, reference=NULL, 
  trans=c("CCStrans","MAtrans"), celnamefile, subset = FALSE)
}
\arguments{
  \item{celfiledir}{Path to the directory in which all .CEL files are saved. Default is the current directory.}
  \item{outfiledir}{Path to the directory in which output will be saved. Default is the current directory.}
  \item{allid}{SNP probes position at the .CEL file. } 
  \item{ABid}{A and B allele name and position. }
  \item{chrid}{chromosome and position of each probe set. }
 \item{CGFLcorrection}{Spline regression coefficients to correct the intensity variation due to different C or G contents in probe sequence and restriction fragment length. Default is not correcting intensity variation.}
  \item{reference}{Reference distribution used for the quantile normalization. If there is an evidence that the density distribution of samples follows more than one distribution or different from that of classical laboratory strain, this normalization should not be used and leave as a default i.e. {\tt reference = NULL}.}
  \item{trans}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
 \item{celnamefile}{A tab delimited .txt file containing the name of .CEL files. .CEL file name should be listed at the first column and the header MUST be 'celfile'. Optionally gender information can be provided by user, and if so it should be listed at the second column with header 'gender'. 'male' and 'female' (CASE SENSITIVE) only recognized and anything else (such as 'FEMALE', 'Male', 'Unknown' or 'NA') should be considered as unknown gender information. If 'filenames' file is not specified, the function reads all .CEL files under the celfiledir.}
  \item{subset}{If this is TRUE, it only genotypes good probe sets. The definition of 'good' probe sets might be vary across laboratories, and current 'good' probe sets were identified by .CEL files processed at The Jackson Laboratory. }
}

\value{
 This function saved $m X n$ matrix rawdataMMchr* and rawdataSSchr*, where $m$ is number of probe sets, $n$ is number of array, and * is chromosome index (1-19, X, Y, and M (Mitochondria)). It also returns length $m$ male indicator, i.e. when it is TRUE (FALSE) it indicates the corresponding array is from male (female) sample. 
 \item{rawdataMMchr}{Normalized log2 contrast intensities.}
 \item{rawdataSSchr}{Normalized log2 average intensities.}
}

\examples{
  \dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  isMale = ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection, reference, 
    trans="CCStrans", subset = TRUE)
  }
}

\keyword{genotype}
