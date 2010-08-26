\author{Hyuna Yang}
\name{densityplot}
\alias{densityplot}
\title{plot density of probe set intensities.}

\description{ 
  This function plot the density of probe set intensities. 
}

\usage{
densityplot(filenames, allid, ABid, type = c("Average", "MatchedSet"), plotname = "densityplot.jpg")
}

\arguments{
  \item{filenames}{Name of .CEL files that you want to draw density plot.}
  \item{allid}{SNP probes position at the .CEL file. } 
  \item{ABid}{A and B allele name and position. }
  \item{type}{Type of SNP intensity that you want to draw a density plot. 'Average' option draws the density plot using the average intensities, and 'MatchedSet' option draws the density plot using higher (matched) intensities between A and B allele intensities.}
  \item{plotname}{Name of density file. The function returns .jpg file, thus the file name should be *.jpg. The default is densityplot.jpg. }
}

\examples{
\dontrun{
load('allid'); load('ABid')
filenames = c('~/projects/cel/SNPex1.CEL', '~/projects/cel/SNPex1.CEL')
densityplot(filenames, allid, ABid, plotname = 'mydensityplot.jpg')
}}

\keyword{plot}
