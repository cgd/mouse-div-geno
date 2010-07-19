\author{ Hyuna Yang }

\name{MouseDivGeno-internal}

\alias{ccstrans}
\alias{computegender}
\alias{test123}
\alias{test12}
\alias{test13}
\alias{vdist}
\alias{E.step3}
\alias{M.step3}
\alias{E.step2}
\alias{M.step2}
\alias{bicov}
\alias{bivar}
\alias{CNVinput}
\alias{igp}
\alias{pennCNVinput}
\alias{simpleCNVdata}
\alias{simpleCNVsummary}

\usage{
ccstrans(a,b,k=4)
computegender(intX,intY, autoint)
test123 (tscore2, tscore3, nm, tgeno2, tgeno3, rmid, thres, iig, nsize)
test12( tscore2, nm, tgeno, rmid, thres, iig, nsize)
test13( tscore3, nm, tgeno3, rmid, geno)
vdist(ta,tb,t2)
E.step3(theta,data)
M.step3(T,data)
E.step2(theta,data)
M.step2(T,data)
bicov(x,y)
bivar(x)
CNVinput(nm, ns, geno, mm, ms, ss, nsize, ngeno, iig)
igp(celfiledir, outfiledir, outfilename, id, chrid, mpos, ename, CGFLcorrection, reference, filenames, mchr)
pennCNVinput(chrid,mpos, exon1info, exon2info, celfiledir, filenames, snpsavefiledir, exonoutfiledir, cnvoutfiledir, mchr)
simpleCNVdata(allid, ABid, chrid, CGFLcorrection, reference, exon1info, exon2info, celfiledir, filenames, cnvoutfiledir, mchr)
simpleCNVsummary(cnvoutfiledir, filenames, refid, mchr, stname)
}

\title{Internal MouseDivGeno functions}

\description{
  Internal MouseDivGeno functions. These are generally not to be called by
  the user. 
}


\keyword{internal}
