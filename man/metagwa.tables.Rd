\name{metagwa.tables}
\alias{metagwa.tables}
\title{Pairwise meta-analysis of results from GWA scans}
\description{
Performes meta-analysis of results of two individual GWA studies
}

\usage{
  metagwa.tables(data.x, data.y, name.x = "P1", name.y = "P2",
               precorrect=TRUE, correct.pooled = FALSE)
}
\arguments{
  \item{data.x}{First data frame with GWA data}
  \item{data.y}{Second data frame with GWA data}
  \item{name.x}{First study name}
  \item{name.y}{Second study name}
  \item{precorrect}{Should GC be applied to the original data before pooling}
  \item{correct.pooled}{Whether to apply Genomic Control correction to
    the study named "POOLED"}
}

\details{
  Original data frames should contain a number of variables,
  such as allelic coding, code for the effective allele,
  etc. Please refer to formetascore function of GenABEL package
  for details.
}

\value{
  A data frame containing summary statistics and essential details
  of the individual studies
}

%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{\code{\link{metagwa.files}}}
%\examples{}
\keyword{htest}
