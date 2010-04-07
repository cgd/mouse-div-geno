\author{Hyuna Yang}
\name{MouseDivGenotype}
\alias{MouseDivGenotype}
\title{Genotype the Mouse Diversity Array}
\description{
  This function reads the .CEL files, normalizes the intensities, genotypes and vinotypes the Mouse Diversity Array. 

  At the normalization step, user can correct hybridization bias due to C or G contents in 25 probe sequences and restriction fragment length correction. Quantile normalization can be done based on one reference distribution. 

  The genotype is based on contrast and summation dimension, and the function offers two different transformations to obtain contrast; contrast centers stretch (CCS) transformation propsed by by BRLMM-P algorithm obtains contrast by asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele, K is hyperparameter, and MA transformation obtains contrast by log2(A) - log2(B).

  This function is highly tailored to Mouse Diversity Array, and if you want to genotype other types of arrays, refer to \code{\link[MouseDivGeno]{genotypethis}} function. User must obtain all annotation files (allid, ABid and chrid) from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. C or G contents and fragment length correction, and quantile normalization is optional, and if you want to do it, you also need to download {\tt CGFLcorrection} and {\tt reference} files from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. Hint file used in the genotyping procedure is also optional and if you want to use 'hint' pre-calculated based on 350 .CEL file, you also need to obtain it from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. 
}
\usage{
MouseDivGenotype(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection=NULL, reference=NULL, hint = NULL,
  trans=c("CCStrans","MAtrans"), celnamefile =NULL, mchr = c(1:19, 'X', 'Y', 'M'), subset = FALSE)
}

\arguments{
  \item{celfiledir}{Path to the directory in which all .CEL files are saved. Default is the current directory.}
  \item{outfiledir}{Path to the directory in which output will be saved. Default is the current directory.}
  \item{allid}{CDF file identifies SNP position. Must be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. }
  \item{ABid}{CDF file identifies A and B allele name and position. Must be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}.}
  \item{chrid}{CDF file specifies indexes for each chromosome. Must be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. }
  \item{CGFLcorrection}{Spline regression coefficients to correct the intensity bias due to different C or G contents in probe sequence and restriction fragment length. The coefficients were pre-calculated based on 350 .CEL files, and can be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. Default is not correcting intensity bias.}
  \item{reference}{Reference distribution used for the quantile normalization. It can be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. If there is an evidence that the .CEL file density distributions follow more than one population, this normalization should not be used and leave as a default (NULL).}
 \item{hint}{Pre-calculated the mean of heterozygosity group in each probe set. Hint file is optional and can be downloaded at {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. }
  \item{trans}{The MouseDivGenotype package genotype based on contrast (A allele intensity - B allele intensity) and average intensity dimensions rather than A allele and B allele dimension. To transform the data, it allows MA transformation, where M = A allele intensity - B allele intensity and A = average intensity, and CCS transformation, a slightly modified polar transformation proposed by BRLMM-P software.}
 \item{celnamefile}{A tab delimited file .txt file containing the name of .CEL files. .CEL file name should be listed at the first column and the header MUST be 'celfile'. Optionally gender information can be provided by user, and if so it should be listed at the second column with header name 'gender'. 'Male' and 'Female' (CASE SENSITIVE) only recognized and anything else (such as 'Female', 'Male', 'Unknown' or 'NA') should be considered as unknown gender information. If 'filenames' file is not specified, the function reads all .CEL files under the celfiledir.}
\item{mchr}{Index of chromosomes that you want to genotype. Chromosome 1 through 19 is indicated as number, X and Y chromosomes are 'X' and 'Y', and Mitochondria is indicated as 'M'. }
\item{subset}{If this is TRUE, it only genotypes good probe sets. The definition of 'Good' probe sets can be vary across laboratory, and current 'Good' probe sets were identified by .CEL files processed at the Jackson Laboratory. }
}

\examples{
\dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference'); load('hint')
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  MouseDivGenotype(dir, outfiledir, allid, ABid, chrid, CGFLcorrection = CGFLcorrection, reference = reference, hint=hint, 
   trans="CCStrans", mchr = c(19, 'X', 'Y', 'M'), subset = TRUE)
  }
}

\keyword{genotype}
