\author{Hyuna Yang}
\name{mouseDivDensityPlot}
\alias{mouseDivDensityPlot}
\title{plot density of probe set intensities.}

\description{ 
  This function plot the density of probe set intensities. 
}

\usage{
mouseDivDensityPlot(celFilenames, snpProbeInfo, type = c("Average", "MatchedSet"))
}

\arguments{
    \item{celFilenames}{Names of .CEL files that you want to draw density plot for.}
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)

        \code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales
    }
    \item{type}{Type of SNP intensity that you want to draw a density plot. 'Average' option draws the density plot using the average intensities, and 'MatchedSet' option draws the density plot using higher (matched) intensities between A and B allele intensities.}
}

\examples{
\dontrun{
load("MouseDivData.RData") # loads MouseDivData into scope
celFilenames <- c("~/projects/cel/SNPex1.CEL", "~/projects/cel/SNPex1.CEL")
mouseDivDensityPlot(celFilenames, MouseDivData)
}}

\keyword{plot}
