\name{forestplot}
\alias{forestplot}
\title{Forest plot}
\description{
  Function to draw meta-analysis forest plots
}

\usage{
  forestplot(estimate, se, labels = paste("Study", c(1:length(estimate))),
           CI = 0.95, xexp = FALSE, ...)
}

\arguments{
  \item{estimate}{Vector of effect estimates}
  \item{se}{Vector of standard errors}
  \item{labels}{Individual study labels}
  \item{CI}{Confidence interval, e.g. 0.95}
  \item{xexp}{Whether the effect values are to be depicted on exponential scale}
  \item{\dots}{Arguments passed to plot function, e.g. main = "My plot"}
}

%\details{
%}
%\value{
%}
%\references{ ~put references to the literature/web site here ~ }
\author{Yurii Aulchenko}
%\note{ ~~further notes~~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
  beta <- c(0.16, 0.091, 0.072, -0.03)
  se   <- c(0.07, 0.042, 0.048, 0.12)
  forestplot(beta, se, main="Example plot")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
