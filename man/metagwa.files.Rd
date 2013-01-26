\name{metagwa.files}
\alias{metagwa.files}
\title{Meta-analysis of seriens of results GWA scans}
\description{
Performes meta-analysis of results of multiple GWA studies
stored in files.
}

\usage{
  metagwa.files(dir = ".", pops, extens, maf = 5, call = 0.95, phwe = 1e-08,
              precorrect=TRUE, correct.pooled = FALSE)
}
\arguments{
  \item{dir}{Path to directory containing files with GWA results}
  \item{pops}{A vector specifying study names, which make the first
        part of the file names}
  \item{extens}{Fixed extension to the file name}
  \item{maf}{Filter threshold for the absolute number/frequency of minor allele.
        If the parameter is < 1, this MAF threshold is used for filtering before
        meta-analysis; if this parameter is >1, n*MAF -- absolute number of
        the minor allele copies -- is used}
  \item{call}{Filter threshold for SNP call rate}
  \item{phwe}{Filter threshold for SNP P-value for HWE}
  \item{precorrect}{Should GC be applied to the original data before pooling}
  \item{correct.pooled}{Whether to apply Genomic Control correction to
    the study named "POOLED"}
}

\details{
  The function looks for files named "POPSextens" in the directory
  "dir", reads them and pool consecutively using
  \code{\link{metagwa.tables}} function.

  The source files should contain a number of variables, such as allelic
  coding, code for the effective allele, etc. Please refer to the
  formetascore function of GenABEL package for details.

  When the file is not available for the first population, or there are
  non-unique population names, the function stops with an error message.

  If there are no files corresponding to next populations, these are
  skipped in analysis with a warning message (see details in the
  description of returned values).
}

\value{
  Results are dumped to file named "POOLEDextens" in the directory
  "dir".  Additionally, a list is returned with elements
  \item{analysed.pops}{list of populations successfully analysed}
}

%\references{ ~put references to the literature/web site here ~ }
\author{Yurii Aulchenko}
%\note{ ~~further notes~~
%}
\seealso{\code{\link{metagwa.tables}}}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
