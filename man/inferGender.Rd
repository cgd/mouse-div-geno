\author{Hyuna Yang and Keith Sheppard}
\name{inferGender}
\alias{inferGender}
\alias{inferGenderFromSnpContAvg}
\title{Infer Sample Gender}

\description{ 
  infer sample gender from some SNP probe intensity statistics
}

\usage{
inferGender(meanIntensityXPerArray, meanIntensityYPerArray, meanIntensityPerAutosome)
inferGenderFromSnpContAvg(snpContAvg, snpInfo)
}

\arguments{
    \item{meanIntensityXPerArray}{a vector containing the mean normalized X chromosome intensity for each sample}
    \item{meanIntensityYPerArray}{a vector containing the mean normalized Y chromosome intensity for each sample}
    \item{meanIntensityPerAutosome}{a vector containing the mean normalized autosome intensity for each sample}
    \item{snpContAvg}{
        an object of class \code{snpContAvg}. See documentation for
        \code{makeSnpContAvg(...)} for details about this class
    }
    \item{snpInfo}{
        A data-frame containing SNP-level metadata data corresponding to \code{snpContAvg}.
        The data-frame should contain the following components:
        
        \code{snpId}: the ID of this SNP
        
        \code{chrId}: the chromosome that this SNP belongs to. Like \dQuote{15} or \dQuote{X}
        
        \code{snpHetHint} (optional): provides a hint for the normalized mean intensity of heterozygous SNPs
        
        \code{isPAR} (optional): if \code{TRUE} this indicates that this SNP is
        in the pseudoautosomal region (PAR). This should only be set to
        \code{TRUE} for SNPs where \code{snpInfo$chrId == "X"}
    }
}
\value{
  a logical vector set to true for samples that were inferred as male and false for samples that were inferred as female
}

%\keyword{}
