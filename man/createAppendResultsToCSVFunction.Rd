\author{Hyuna Yang}
\name{createAppendResultsToCSVFunction}
\alias{createAppendResultsToCSVFunction}
\title{Create a Function for Incrementally Saving Genotyping/Vinotyping Data}
\description{
Create a function for incrementally saving genotyping/vinotyping data. See
mouseDivGenotype documentation for example usage.
}
\usage{
createAppendResultsToCSVFunction(fileName)
}

\arguments{
    \item{fileName}{the file that we should write results to}
}

\value{
    a function which will take incremental genotyping/vinotyping results
    and append them to \code{fileName}. This function can be passed into
    this mouseDivGenotype function via the processResultsFunction in order
    to incrementally save genotyping results to file.
}
