\author{Hyuna Yang}

\name{imageplot}
\alias{imageplot}

\title{.CEL file image plot}

\description{
  This function generates a image plot of .CEL file.
}

\usage{
imageplot(filename, plotname = 'imageplot.jpg', celfiledir, outfiledir)
}

\arguments{
  \item{filenames}{Name of a .CEL file that you want to draw image plot.}
  \item{plotname}{Name of density file. The function returns .jpg file, thus the file name should be *.jpg. The default is imageplot.jpg. }
  \item{celfiledir}{Path to the directory in which the .CEL files is saved. Default is the current directory.}
  \item{outfiledir}{Path to the directory in which plot will be saved. Default is the current directory.}
}

\examples{
\dontrun{
filename = '~/projects/cel/SNPex1.CEL'
celfiledir = '~/projects/cel'
outfiledir = '~/projects/outfile'
densityplot(filename, plotname = 'myimageplot.jpg', celfiledir = celfiledir, outfiledir = outfiledir)
}}


\keyword{plot}
