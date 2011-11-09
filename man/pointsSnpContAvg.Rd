\author{Hyuna Yang and Keith Sheppard}
\name{pointsSnpContAvg}
\alias{pointsSnpContAvg}

\title{pointsSnpContAvg}

\description{
    Scatter-plot the intensity contrasts and averages in the given snpContAvg object.
    These points are rendered in an existing plot as opposed to
    \code{plotSnpContAvg} which creates a new plot.
}

\usage{
pointsSnpContAvg(
        snpContAvg,
        renderMask=matrix(TRUE, nrow=nrow(snpContAvg), ncol=ncol(snpContAvg)),
        col=rgb(0, 0, 1),
        pch=20,
        ...)
}

\arguments{
    \item{snpContAvg}{
        the SNP contrast and average data. This must be of S3 class snpContAvg
        as defined by the \code{makeSnpContAvg(...)} function
    }
    \item{renderMask}{
        for each point in the given \code{snpContAvg} parameter the \code{renderMask}
        determines whether that point should be rendered or not. The default matrix
        of all TRUE will render every point
    }
    \item{col}{the color as defined by the \code{par} function}
    \item{pch}{the symbol to plot as defined by the \code{par} function}
    \item{...}{all remaining arguments which are passed on to the \code{points} function}
}

\value{
}

\examples{
\dontrun{
  }
}

\keyword{plot}
