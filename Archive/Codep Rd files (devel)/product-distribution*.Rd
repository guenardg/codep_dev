\encoding{utf8}
\name{Product-distribution}
\alias{Product-distribution}
\alias{dphi}
\alias{pphi}
\alias{dtau}
\alias{ptau}
\title{Frequency distributions for MCA parametric testing}
\description{
  Density and distribution functions of the phi statistic, which is the
  product of two Fisher-Snedecor distributions with particular degrees
  of freedom.
}
\usage{
dphi(x, nu1, nu2, tol = .Machine$double.eps ^ 0.5)
pphi(q, nu1, nu2, lower.tail = TRUE, tol = .Machine$double.eps ^ 0.5)
dtau(x, nu, tol = .Machine$double.eps ^ 0.5)
ptau(q, nu, lower.tail = TRUE, tol = .Machine$double.eps ^ 0.5)
}
\arguments{
  \item{x, q}{ vector of quantile. }
  \item{nu1, nu2, nu}{ degrees of freedom (>0, may be
    non-integer). \code{Inf} is allowed. }
  \item{lower.tail}{ logical; if TRUE (default), probabilities are P[X
    <= x], otherwise, P[X > x]. }
  \item{tol}{ tolerance used for numerical estimation. }
}
\value{
  \code{dphi} and \code{dtau} return the density functions whereas
  \code{pphi} and \code{ptau} return the distribution functions.
}
\details{
  The density distribution of a variable \code{z} that is the product of
  two random variables \code{x} and \code{y} with density distributions
  f(x) and g(y), respectively, is the integral over the intersection of
  the domains of \code{x} and \code{y} of f(x) * g(z/x) / abs(x) dx.

  \code{dphi} estimates density values using numerical integration
  (\code{\link{integrate}}) of the Fisher-Scedecor \code{\link{FDist}} density
  distribution function. Following the algebra of Multiscale
  Codependence Analysis, f(x) has df1 = nu1 and df2 = nu1 * nu2 degrees
  of freedom; g(x) has 'df1 = 1' and 'df2 = nu2' degrees of
  freedom. Hence, the distribution of that product has two parameters.

  \code{pphi} integrates \code{dphi} in the interval
  [\code{0},\code{q}] when "lower.tail = TRUE" (the default) and on
  the interval [\code{q},\code{Inf}] when "lower.tail = FALSE"".

  \code{dtau} and \code{ptau} are similar to \code{dphi} integrates
  \code{pphi}, but with f(x) and f(y) being two Student's t distributions
  with \code{nu} degrees of freedom. \code{} It is called by functions
  \code{\link{test.cdp}} and \code{\link{permute.cdp}} to perform tests
  of hypotheses for single response variables, in which case
  one-tailed tests can be performed.
}
\references{
  Springer, M. D. 1979. The algebra of random variables. John Wiley and
  Sons Inc., Hoboken, NJ, USA.

  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964

  Guénard, G. Legendre, P. 2018. Bringing multivariate support to
  multiscale codependence analysis: Assessing the drivers of community
  structure across spatial scales. Meth. Ecol. Evol. 9: 292-304
}
\seealso{
  \link{test.cdp}
}
\examples{
#
### Display the phi probability distribution for five different pairs
### of degrees of freedom
#
x <- 10^seq(-4, 0.5, 0.05)
plot(y = dphi(x, 1, 10), x = x, type = "l", col = "black", las = 1, ylab = "pdf",
  ylim = c(0, 0.5))
lines(y = dphi(x, 3, 10), x = x, col = "purple")
lines(y = dphi(x, 5, 70), x = x, col = "blue")
lines(y = dphi(x, 12, 23), x = x, col = "green")
lines(y = dphi(x, 35, 140), x = x, col = "red")
#
### Display the density distribution function for 10 degrees of freedom
### and the cumulative probability above x = 1
#
x <- 10^seq(-4, 0.5, 0.05)
y <- dphi(x, 5, 70)
plot(y = y, x = x, type = "l", col = "black", las = 1, ylab = "Density",
  ylim = c(0, 0.5))
polygon(x = c(x[81L:91], x[length(x)], 1), y = c(y[81L:91], 0, 0),
  col = "grey")
text(round(pphi(1, 5, 70, lower.tail=FALSE), 3), x = 1.75, y = 0.05)
#
### Same for the tau distribution
#
x <- c(-(10^seq(0.5, -4, -0.05)), 10^seq(-4, 0.5, 0.05))
plot(y = dtau(x, 1), x = x, type = "l", col = "black", las = 1,
   ylab = "pdf", ylim = c(0, 0.5))
lines(y = dtau(x, 2), x = x, col = "purple")
lines(y = dtau(x, 5), x = x, col="blue")
lines(y = dtau(x, 10), x = x, col="green")
lines(y = dtau(x, 100), x = x, col="red")
#
y <- dtau(x, 10)
plot(y = y, x = x, type = "l", col = "black", las = 1, ylab = "Density",
  ylim = c(0, 0.5))
polygon(x = c(x[which(x==1):length(x)], x[length(x)],1),
  y = c(y[which(x==1):length(x)], 0, 0), col = "grey")
text(round(ptau(1, 10, lower.tail = FALSE), 3), x = 1.5, y = 0.03)
polygon(x = c(-1, x[1L], x[1L:which(x==-1)]),
  y = c(0, 0, y[1L:which(x==-1)]), col="grey")
text(round(ptau(-1, 10), 3), x = -1.5, y = 0.03)
#
}
\author{ Guillaume Guénard, Département de sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
