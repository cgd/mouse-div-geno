\author{Hyuna Yang}
\name{ReadCelFile}
\alias{ReadCelFile}
\title{Normalized the Mouse Diversity Array}
\description{
  This function reads the .CEL files, normalizes the intensities, and saves normalized intensities. Refer to refer to the help file of \code{\link[MouseDivGeno]{MouseDivGenotype}} for detail. 
}
\usage{
ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection=NULL, reference=NULL, 
  trans=c("CCStrans","MAtrans"), celnamefile, subset = FALSE )
}
\arguments{
  \item{celfiledir}{Path to the directory in which all .CEL files are saved. Default is the current directory.}
  \item{outfiledir}{Path to the directory in which output will be saved. Default is the current directory.}
  \item{allid}{CDF file identifies SNP position. Must be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. }
  \item{ABid}{CDF file identifies A and B allele name and position. Must be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}.}
  \item{chrid}{CDF file specifies indexes for each chromosome. Must be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. }
  \item{CGFLcorrection}{Spline regression coefficients to correct the intensity bias due to different C or G contents in probe sequence and restriction fragment length. The coefficients were pre-calculated based on 350 .CEL files, and can be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. Default is not correcting intensity bias.}
  \item{reference}{Reference distribution used for the quantile normalization. It can be obtained from {\tt http://genomedynamics.org/tools/MouseDivGeno/CDFfiles}. If there is an evidence that the .CEL file density distributions follow more than one population, this normalization should not be used and leave as a default (NULL).}
 \item{trans}{The MouseDivGenotype package genotype based on contrast (A allele intensity - B allele intensity) and average intensity dimensions rather than A allele and B allele dimension. To transform the data, it allows MA transformation, where M = A allele intensity - B allele intensity and A = average intensity, and CCS transformation, a slightly modified polar transformation proposed by BRLMM-P software.}
 \item{celnamefile}{A tab delimited file .txt file containing the name of .CEL files. .CEL file name should be listed at the first column and the header MUST be 'celfile'. Optionally gender information can be provided by user, and if so it should be listed at the second column with header name 'gender'. 'Male' and 'Female' (CASE SENSITIVE) only recognized and anything else (such as 'Female', 'Male', 'Unknown' or 'NA') should be considered as unknown gender information. If 'filenames' file is not specified, the function reads all .CEL files under the celfiledir.}
\item{subset}{If this is TRUE, it only genotypes good probe sets. The definition of 'Good' probe sets can be vary across laboratory, and current 'Good' probe sets were identified by .CEL files processed at the Jackson Laboratory. }
}

\value{
 \item{isMle}{Computed gender information.}
}

\examples{
  \dontrun{
  load('allid'); load('ABid'); load('chrid'); load('CGFLcorrection'); load('reference');
  celfiledir = '~/projects/cel'
  outfiledir = '~/projects/outfile'
  ReadCel(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection, reference, 
    trans="CCStrans", subset = TRUE)
  }
}

\keyword{genotype}
