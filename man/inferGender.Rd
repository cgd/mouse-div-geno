\author{Hyuna Yang and Keith Sheppard}
\name{inferGender}
\alias{inferGender}
\title{Infer Sample Gender}

\description{ 
  infer sample gender from some SNP probe intensity statistics
}

\usage{
inferGender(meanIntensityXPerArray, meanIntensityYPerArray, meanIntensityPerAutosome)}

\arguments{
    \item{meanIntensityXPerArray}{a vector containing the mean normalized X chromosome intensity for each sample}
    \item{meanIntensityYPerArray}{a vector containing the mean normalized Y chromosome intensity for each sample}
    \item{meanIntensityPerAutosome}{a vector containing the mean normalized autosome intensity for each sample}
}
\value{
  a logical vector set to true for samples that were inferred as male and false for samples that were inferred as female
}

%\keyword{}
