\author{Hyuna Yang}
\name{vinotype}
\alias{vinotype}

\title{Identify copy number variance.}

\description{This function is to identify copy number variance using a reference sample. It uses \code{simpleCNVdata} to get the intensities from SNP and exon, and \code{simpleCNVsummary} which uses {\tt HiddenMarkov}, a Hidden Markov Model R package, to identify copy number gain or loss by comparing a sample and a reference sample.}

\usage{simpleCNP(allid, ABid, chrid, CGFLcorrection, reference, exon1info, exon2info, celfiledir, filenames, 
cnvoutfiledir, mchr, stname, refsample) }

\arguments{
  \item{allid}{SNP probes position at the .CEL file. } 
  \item{ABid}{A and B allele name and position. }
  \item{chrid}{chromosome and position of each probe set. }
  \item{CGFLcorrection}{Spline regression coefficients to correct the intensity variation due to different C or G contents in probe sequence and restriction fragment length.}
  \item{reference}{Reference distribution used for the quantile normalization.}
  \item{exon1info}{Exon1 and exon2 are conventional name that Mouse Diversity Array used to distinguish types of exon. Exon1 and exon2 are invariant 25mer on exon, and if restriction fragment length is less than 1Kb, it is conventiallyed called an exon1 probes otherwise exon2 probes. Typically exon2 is dimmer than exon1, and exon1info and exon2info file provide position of probe sets, normalization coefficients, etc. If user wants to include exon1 for the CNV analysis, exon1info should be identified here. }
  \item{exon2info}{If user wants to include exon2 for the CNV analysis, exon2info should be identified here. }
  \item{celfiledir}{Path to the directory in which all .CEL files are saved. Default is the current directory.}
  \item{filenames}{A vector of .CEL file names.}
  \item{cnvoutfiledir}{Path to the directory in which all intensities from SNP, exon1 and exon2 are saved. Default is the current directory.}
   \item{mchr}{Index of chromosomes that you want to genotype}
   \item{stname}{Strain names corresponding to the .CEL files listed in the {\tt filenames}. If it is not specified, .CEL file names are used.}
  \item{refsample}{.CEL file name (must be one of .CEL files listed in the {\tt filenames}) of reference strain to which every other sample is compared. }
}
\value{
  It saves intensities from SNP, exon1 and exon2 by chromosomes, HMM states (copy number gain = 1 and loss = 3 compared to the reference strain), and summary files only showing copy number variance region. The summary table shows the name of strain, status (gain or loss), start and end position of CNV, number of total, SNP, exon1 and exon2 probe sets in that region, mean intensity of the  reference and a sample. 
}
\examples{
\dontrun{
load('allid'); load('ABid'); load('chrid');load('CGFLcorrection'); load('reference');
load('exon1info');load('exon2info');
celfiledir = '~/celfiles/'; cnvoutfiledir = '~/stest2/'
filenames = c('ex1.CEL','ex2.CEL','ex3.CEL'); stname =c('B6','DBA/2','A/J')
mchr = 1:19
refsample = 'ex1.CEL'
simpleCNV(allid, ABid, chrid, CGFLcorrection, reference, exon1info, exon2info, celfiledir, filenames, 
cnvoutfiledir, mchr = c(1:19), stname, refsample)}
}

\keyword{simpleCNV}
