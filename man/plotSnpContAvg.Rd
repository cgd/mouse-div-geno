\author{Hyuna Yang and Keith Sheppard}
\name{plotSnpContAvg}
\alias{plotSnpContAvg}

\title{plotSnpContAvg}

\description{
    Scatter-plot the intensity contrasts and averages in the given snpContAvg object.
    These points are rendered to a new plot as opposed to
    \code{plotSnpContAvg} which renders to an existing plot.
}

\usage{
plotSnpContAvg(
        snpContAvg,
        renderMask=matrix(TRUE, nrow=nrow(snpContAvg), ncol=ncol(snpContAvg)),
        xlab="Contrast",
        ylab="Average Intensity",
        main="",
        col=rgb(0, 0, 0, 0.25),
        pch=20,
        xlim=c(-max(snpContAvg$snpContrasts), max(snpContAvg$snpContrasts)),
        ylim=c(0, max(snpContAvg$snpAverages)),
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
    \item{xlab}{the X axis label as defined by the \code{plot} function}
    \item{ylab}{the Y axis label as defined by the \code{plot} function}
    \item{main}{the main title as defined by the \code{plot} function}
    \item{col}{the color as defined by the \code{par} function}
    \item{pch}{the symbol to plot as defined by the \code{par} function}
    \item{xlim}{a vector of two numbers defining the minimum and maximum extent of the X axis}
    \item{ylim}{a vector of two numbers defining the minimum and maximum extent of the Y axis}
    \item{...}{all remaining arguments which are passed on to the \code{plot} function}
}

\value{
}

\examples{
\dontrun{
  }
}

\keyword{plot}
