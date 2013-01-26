\name{forestplot}
\alias{forestplot}
\title{Function to draw meta-analysis forest plots}
\usage{
forestplot(estimate, se, labels = paste("Study", c(1:length(estimate))),
  CI = 0.95, xexp = FALSE, xlab = expression(beta), ylab = "", ...)
}
\arguments{
  \item{estimate}{Vector of effect estimates}

  \item{se}{Vector of standard errors}

  \item{labels}{Vector of labels for the individual studies
  (default: Study 1, Study 2, etc.)}

  \item{CI}{Confidence interval (default: 0.95)}

  \item{xexp}{Whether the effect values are to be depicted
  on an exponential scale (default: \code{FALSE})}

  \item{xlab}{Label for the horizontal axis (default:
  \eqn{\beta})}

  \item{ylab}{Label for the horizontal axis (default:
  empty)}

  \item{...}{Arguments passed to the \code{plot} function,
  e.g. \code{main="My plot"}}
}
\description{
This function creates forest plots from meta-analysis data.
}
\examples{
beta <- c(0.16, 0.091, 0.072, -0.03)
se   <- c(0.07, 0.042, 0.048, 0.12)
forestplot(beta, se, main="Example plot")
}
\author{
Yurii Aulchenko, Lennart C. Karssen
}
\keyword{hplot}

