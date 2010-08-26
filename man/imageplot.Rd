\author{Hyuna Yang}

\name{imageplot}
\alias{imageplot}

\title{.CEL file image plot}

\description{
  This function generates a image plot of .CEL file.
}

\usage{
imageplot(filename, plotname = 'imageplot.jpg')
}

\arguments{
  \item{filename}{Name of a .CEL file that you want to draw image plot.}
  \item{plotname}{Name of density file. The function returns .jpg file, thus the file name should be *.jpg. The default is imageplot.jpg. }
}

\examples{
\dontrun{
filename = '~/projects/cel/SNPex1.CEL'
densityplot(filename, plotname = 'myimageplot.jpg')
}}


\keyword{plot}
