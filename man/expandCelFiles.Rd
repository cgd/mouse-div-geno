\author{Hyuna Yang}
\name{expandCelFiles}
\alias{expandCelFiles}
\title{Expand CEL Files}
\description{
    This function will expand a directory to a vector of the CEL files which are
    found in that directory. If given a vector of CEL files the vector will be
    returned unchanged.
}
\usage{
expandCelFiles(filename)
}

\arguments{
    \item{filename}{the directory or CEL file vector}
}

\value{
    If \code{filename} is a directory the vector of CEL file in the directory
    will be returned. If given a vector of CEL files this function will return
    the vector unchanged.
}
