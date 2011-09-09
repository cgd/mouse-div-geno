\author{Hyuna Yang and Keith Sheppard}
\name{normalizeCelFileByChr}
\alias{normalizeCelFileByChr}
\title{Normalize a CEL File and Break Into Chromosome Chunks}

\description{
    Normalizes CEL file SNP data and calculates the average and contrast for
    each SNP. The result is returned as a per-chromosome list where each list
    element contains two matrices (one matrix for SNP average and the other
    for contrast)
}

\usage{
normalizeCelFileByChr(
    celFileName,
    snpProbeInfo,
    snpInfo,
    logFile = NULL,
    chromosomes = c(as.character(1 : 19), "X", "Y", "M"),
    referenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans", "None"),
    keepSNPIds = FALSE)
}

\arguments{
    \item{celFileName}{the file path to the CEL file that will be read and normalized}
    \item{snpProbeInfo}{
        A data-frame containing probe-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{isAAllele}: \code{TRUE} for A allele probes and \code{FALSE} for B allele probes
        
        \code{snpId}: the ID of the SNP that this probe belongs to. These IDs should correspond to the same component in the \code{snpInfo} argument
        
        \code{probeIndex}: the index for the probe (used to get intensity data)
        
        \code{correction} (optional): if present this correction will be applied by adding it to probe mean intesity vales
    }
    
    \item{snpInfo}{
        A data-frame containing SNP-level data for the SNPs we are genotyping.
        The data-frame should contain the following components:
        
        \code{snpId}: the ID of this SNP
        
        \code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}
        
        \code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygous SNPs
        
        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is in the pseudoautosomal region (PAR). This should only be set to \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
    }
    \item{logFile}{if non-NULL this must be a connection type to which log messages will be written}
    \item{chromosomes}{the vector of chromosomes that should be normalized. Data from other chromosomes will be discarded}
    \item{referenceDistribution}{Reference distribution used for the quantile normalization. If there is an evidence that the density distribution of samples follows more than one distribution or different from that of classical laboratory strain, this normalization should not be used and leave as a default i.e. \code{referenceDistribution = NULL}.}
    \item{transformMethod}{The MouseDivGenotype package genotype the array based on contrast and summation, and offers two transformations to obtain contrast; CCS and MA transformation. CCS (contrast centers stretch) transformation proposed by BRLMM-P algorithm obtains contrast via asinh(K*(A-B)/(A+B))/asinh(K), where A and B is intensity of A and B allele and K is a hyperparameter. MA transformation obtains contrast by log2(A) - log2(B). In both cases summation is defined by (log2(A) + log2(B))/2}
    \item{keepSNPIds}{determines if the return value should preserve the SNP IDs or not}
}
\value{
    A list object where chromosome names can be used to index. Each list item
    will contain two matrix components:
    \item{M}{a matrix of (normalized) intensity contrasts. There will be a row per-SNP and a column per-sample. SNP ordering will be consistent with the ordering used for S component.}
    \item{S}{a matrix of (normalized) intensity averages. There will be a row per-SNP and a column per-sample. SNP ordering will be consistent with the ordering used for M component.}
    The two matrix components along with the chromosome index can be used as
    inputs to the genotypeAnyChrChunk function in order to genotype and
    vinotype the SNPs
}

%\keyword{}
