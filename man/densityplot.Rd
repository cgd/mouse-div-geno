\author{Hyuna Yang}
\name{densityplot}
\alias{densityplot}
\title{plot density of probe set intensities.}

\description{ 
  This function plot the density of probe set intensities. 
}

\usage{
densityplot(filenames, plotname)
}

\arguments{
  \item{filenames}{Name of .CEL files that you want to draw density plot.}
  \item{plotname}{Name of density file. The function returns .jpg file, thus the file name should be *.jpg. The default is densityplot.jpg. }
}

\examples{
\dontrun{
filenames = c('~/projects/cel/SNPex1.CEL', '~/projects/cel/SNPex1.CEL')
densityplot(filenames, plotname = 'mydensityplot.jpg')
}}

\keyword{plot}
