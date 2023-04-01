\encoding{utf8}
\name{MCA}
\alias{MCA}
\alias{test.cdp}
\alias{permute.cdp}
\alias{parPermute.cdp}
\title{ Multiple-descriptors, Multiscale Codependence Analysis }
\description{ Functions to perform Multiscale Codependence Analysis
  (MCA) for univariate or multivariate response data }
\usage{
MCA(Y, X, emobj)
test.cdp(object, alpha = 0.05, max.step, response.tests = TRUE)
permute.cdp(object, permute, alpha = 0.05, max.step,
            response.tests = TRUE)
parPermute.cdp(object, permute, alpha = 0.05, max.step,
               response.tests = TRUE, nnode, seeds, verbose = TRUE,
               ...)
}
\arguments{
  \item{Y}{ a numeric matrix or vector containing the response
    variable(s). }
  \item{X}{ a numeric matrix or vector containing the explanatory
    variable(s). }
  \item{emobj}{ a \code{\link{eigenmap-class}} object. }
  \item{object}{ a \code{\link{cdp-class}} object. }
  \item{alpha}{ type I (\eqn{\alpha}{alpha}) error threshold used by the
    testing procedure. }
  \item{max.step}{ maximum number of steps to perform when testing for
    statistical significance. }
  \item{response.tests}{ a boolean specifying whether to test individual
    response variables. }
  \item{permute}{ The number of random permutations used for
    testing. When omitted, the number of permutations is calculated
    using function \code{\link{minpermute}}.}
  \item{nnode}{ The number of parallel computation nodes. }
  \item{seeds}{ Random number generator seeds for parallel the
    computation nodes. }
  \item{verbose}{ Whether to return user notifications. }
  \item{...}{ Parameters to be passed to \code{parallel::makeCluster()}}
}
\details{
  Multiscale Codependence Analysis (MCA) allows the calculation of
  correlation-like (i.e.codependence) coefficients between two variables
  with respect to structuring variables (Moran's eigenvector maps). The
  purpose of this function is limited to parameter fitting. 
  
  Test procedures are handled through \code{test.cdp} (parametric testing)
  or \code{permute.cdp} (permutation testing). Additionaly, methods are
  provided for printing (\code{print.cdp}), displaying a summary of the tests 
  (\code{summary.cdp}), plotting results (\code{plot.cdp}), calculating fitted 
  (\code{fitted.cdp}) and residuals values (\code{redisuals.cdp}), and making
  predictions (\code{predict.cdp}).
  
  It is noteworthy that the test procedures used by
  \code{MCA} differ from the standard R workflow since intermediate
  testing functions (\code{test.cdp} and \code{permute.cdp}) need
  first to be called before testing be performed. 
  ### Éliminer cette dernière phrase qui ne me semble pas très utile: ### For \code{MCA},
  testing functionalities have been moved away from summary.cdp because
  testing is computationally intensive. ###
  
  Function \code{parPermute.cdp} (permutation test using parallel computing)
  allows the user to spread the number of permutation on many
  computation nodes. It relies on package \code{parallel}. Omitting
  parameter \code{nnode} lets function \code{parallel::detectCores()}
  specify the number of nodes. Similarly, omitting parameter \code{seeds}
  makes the function draw seeds uniformly at random between 1 and
  \code{.Machine$integer.max}. If needed, one may pass initialization
  parameters to \code{parallel::makeCluster()}.
}
\value{
  A \code{\link{cdp-class}} object.
}
\references{
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
  
  Guénard, G. and Legendre, P. 2018. Bringing multivariate support to
  multiscale codependence analysis: Assessing the drivers of community
  structure across spatial scales. Meth. Ecol. Evol. 9: 292-304
}
\author{ Guillaume Guénard, Département de sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
\seealso{
  \code{\link{eigenmap}}
}
\examples{
#
###### Begin {Salmon example}
#
data(Salmon)
#
## Converting the data from data frames to matrices;
## drop=FALSE prevents dimension reduction of the matrix "Abundance" to a vector
Abundance <- log1p(as.matrix(Salmon[,"Abundance",drop=FALSE]))
Environ <- as.matrix(Salmon[,3:5])
#
## Creating a spatial eigenvector map from the site positions along the river 
map1 <- eigenmap(x=Salmon[,"Position"],weighting=Wf.binary,boundaries=c(0,20))
summary(map1)
# The eigenvalues are in map1$lambda, the MEM eigenvectors in matrix map1$U
#
## MCA with a single response variable; parametric test
mca1 <- MCA(Y=Abundance,X=Environ[,"Substrate",drop=FALSE],emobj=map1)
summary(mca1)
mca1_partest <- test.cdp(mca1)      # Parametric test
summary(mca1_partest)
par(mar = c(6,4,2,4))
## In the following plot, the significant dbMEMs are shown with a red frame
plot(mca1_partest, las = 3)
#
## Same analysis with permutation test
mca1_pertest <- permute.cdp(mca1)   # Permutation test
\dontrun{     # Code not run during CRAN software tests
## or:
mca1_pertest <- parPermute.cdp(mca1,permute=999999)
}
summary(mca1_pertest)
## In the following plot, the significant dbMEMs are shown with a red frame
plot(mca1_pertest, las = 3)
#
# The array containing the codependence coefficients is found in:
mca1_pertest$UpYXcb$C 
#
## Analysis with all three explanatory variables
mca2 <- MCA(Y=log1p(as.matrix(Salmon[,"Abundance",drop=FALSE])),
            X=as.matrix(Salmon[,3L:5]),emobj=map1)
summary(mca2)
mca2_partest <- test.cdp(mca2)      # Parametric test
summary(mca2_partest)
par(mar = c(6,4,2,4))
## In the following plot, significant dbMEMs for explanatory variables have a red frame
plot(mca2_partest, las = 3)
#
## Same analysis with permutation test
mca2_pertest <- permute.cdp(mca2)   # Permutation test
\dontrun{     # Code not run during CRAN software tests
or:
mca2_pertest <- parPermute.cdp(mca2,permute=999999)
}
summary(mca2_pertest)
## In the following plot, significant dbMEMs for explanatory variables have a red frame
plot(mca2_pertest, las = 3)
#
mca2_pertest$UpYXcb$C      # Array containing the codependence coefficients
test.res.mat <- mca2_pertest$UpYXcb$C[,1L,]   # ...now turned into a matrix
#
###### End {Salmon example}
#
###### Begin {Doubs example}
#
# The Doubs.rda file contains three data files; see ?Doubs
data(Doubs)
#
## Creating a spatial eigenvector map from the DSF = "Distance from source" variable
map2 <- eigenmap(x=Doubs.geo[,"DFS"])
# The eigenvalues are in map2$lambda, the MEM eigenvectors in matrix map2$U
#
## MCA with multivariate response data
mca3 <- MCA(Y=log1p(Doubs.fish),X=Doubs.env,emobj=map2)
mca3
mca3_pertest <- permute.cdp(mca3)
\dontrun{     # Code not run during CRAN software tests
## or:
mca3_pertest <- parPermute.cdp(mca3,permute=999999)
}
mca3_pertest
summary(mca3_pertest)
par(mar = c(6,4,2,4))
plot(mca3_pertest, las = 2)
mca3_pertest$UpYXcb$C # Array containing the codependence coefficients
#
## Display the results along the transect
spmeans <- colMeans(log1p(Doubs.fish))
pca1 <- svd(log1p(Doubs.fish) - rep(spmeans,each=nrow(Doubs.fish)))
par(mar = c(5,5,2,5)+0.1)
plot(y = pca1$u[,1L], x = Doubs.geo[,"DFS"], pch = 21L, bg = "red",
     ylab = "PCA1 loadings", xlab = "Distance from river source (km)")
#
x <- seq(0,450,1)
newdists <- matrix(NA, length(x), nrow(Doubs.geo))
for(i in 1L:nrow(newdists))
  newdists[i,] <- abs(Doubs.geo[,"DFS"] - x[i])
#
## Calculating predictions for arbitrary sites under the same set of
## environmental conditions that the codependence model was built with.
prd1 <- predict(mca3_pertest,
                newdata=list(target = eigenmap.score(map2, newdists)))
#
## Projection of the predicted species abundance on pca1:
Uprd1 <- (prd1 - rep(spmeans, each = nrow(prd1))) \%*\% pca1$v \%*\% diag(pca1$d^-1)
lines(y = Uprd1[,1L], x = x, col=2, lty = 1)
#
## Projection of the predicted species abundance on pca2:
plot(y = pca1$u[,2L], x = Doubs.geo[,"DFS"], pch = 21L, bg = "red",
     ylab = "PCA2 loadings", xlab = "Distance from river source (km)")
lines(y = Uprd1[,2L], x = x, col=2, lty = 1)
#
## Displaying only the observed and predicted abundance for Brown Trout.
par(new=TRUE)
plot(y = log1p(Doubs.fish[,"TRU"]),Doubs.geo[,"DFS"],pch=21L,bg="green",
     ylab="",xlab="",new=FALSE,axes=FALSE)
axis(4)
lines(y = prd1[,"TRU"], x = x, col=3)
mtext(side=4, "log(Abundance+1)", line = 2.5)
#
###### End {Doubs example}
#
###### Begin {Oribatid mite example}
#
# The Mite.rda file contains three data files; see ?Mite
data(Mite)
#
## Create a spatial eigenvector map from the coordinates of the soil cores in mite.geo
map3 <- eigenmap(x = mite.geo)
# The eigenvalues are in map3$lambda, the MEM eigenvectors in matrix map3$U
#
## MCA with multivariate response data and explanatory environmental variables
mca4 <- MCA(Y = log1p(mite.species), X = mite.env, emobj = map3)
mca4_partest <- test.cdp(mca4, response.tests = FALSE)      # Parametric tests
summary(mca4_partest)
## Greyscale plot; significant dbMEMs for explanatory variables have a red frame
plot(mca4_partest, las = 2, lwd = 2)
## Same plot in colour; significant dbMEMs for explanatory variables have a white frame
plot(mca4_partest, col = rainbow(1200)[1L:1000], las = 3, lwd = 4,
     main = "Codependence diagram", col.signif = "white")
#
## *** Guillaume: explain the purpose of the code that follows. What does it do?
rng <- list(x = seq(min(mite.geo[,"x"]) - 0.1, max(mite.geo[,"x"]) + 0.1, 0.05),
            y = seq(min(mite.geo[,"y"]) - 0.1, max(mite.geo[,"y"]) + 0.1, 0.05))
grid <- cbind(x = rep(rng[["x"]], length(rng[["y"]])),
              y = rep(rng[["y"]], each = length(rng[["x"]])))
newdists <- matrix(NA, nrow(grid), nrow(mite.geo))
for(i in 1L:nrow(grid)) {
  newdists[i,] <- ((mite.geo[,"x"] - grid[i,"x"])^2 +
                   (mite.geo[,"y"] - grid[i,"y"])^2)^0.5
}
# Compute mean abundance across the sites for each species
spmeans <- colMeans(mite.species)
# PCA by SVD of ### Data incorrectly centred ###
pca2 <- svd(log1p(mite.species) - rep(spmeans, each = nrow(mite.species)))
#
prd2 <- predict(mca4_partest,
          newdata = list(target = eigenmap.score(map3, newdists)))
Uprd2 <- (prd2 - rep(spmeans, each = nrow(prd2))) \%*\% pca2$v \%*\% diag(pca2$d^-1)
#
### Printing the response variable     ### Printing what? ### 
prmat <- Uprd2[,1L]
dim(prmat) <- c(length(rng$x),length(rng$y))
zlim <- c(min(min(prmat),min(pca2$u[,1L])),max(max(prmat),max(pca2$u[,1L])))
image(z = prmat, x = rng$x, y = rng$y, asp = 1, zlim = zlim,
      col = rainbow(1200L)[1L:1000], ylab = "y", xlab = "x")
points(x = mite.geo[,"x"], y = mite.geo[,"y"], pch = 21,
  bg = rainbow(1200L)[round(1+(999*(pca2$u[,1L]-zlim[1L])/(zlim[2L]-zlim[1L])),0)])
#
###### End {Oribatid mite example}
#
}
