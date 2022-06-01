\encoding{utf8}
\name{eigenmap}
\alias{eigenmap}
\alias{Wf.sqrd}
\alias{Wf.RBF}
\alias{Wf.binary}
\alias{Wf.PCNM}
\alias{Wf.Drayf1}
\alias{Wf.Drayf2}
\alias{Wf.Drayf3}
\alias{eigenmap.score}
\title{Spatial eigenvector maps}
\description{
  Function to calculate spatial eigenvector maps of a set of locations
  in Euclidean space with an arbitrary number of dimensions. A variety of 
  transformations of the geographic distances through weighting functions are proposed.
}
\usage{
eigenmap(x,opt.coord=NA,weighting=Wf.sqrd,boundaries,wpar,select=.Machine$double.eps^0.5)
Wf.sqrd(D)
Wf.RBF(D,wpar=1)
Wf.binary(D,boundaries)
Wf.PCNM(D,boundaries)
Wf.Drayf1(D,boundaries)
Wf.Drayf2(D,boundaries,wpar=1)
Wf.Drayf3(D,boundaries,wpar=1)
eigenmap.score(object,target)
}
\arguments{
  \item{x}{ A set of coordinates defined in one dimension (numeric vector), or in 
    many dimensions (in a matrix of coordinates with any number of dimensions). 
    Or, alternatively, a geographic distance matrix computed by \code{stats} function 
    \code{\link{dist}}. The coordinates are treated as Cartesian coordinates and, if a 
    distance matrix is provided, the distances are assumed to be Euclidean. }
  \item{opt.coord}{ Coordinates to be used for plotting when a distance matrix is
    provided in x. }
  \item{weighting}{ Select one of the following transformation functions to compute the 
    edge weighting matrix. That function must have the raw distances as first
    argument; optionally, a second argument called \code{boundaries} is provided,
    giving the smallest (usually 0) and the largest distance within which locations are 
    considered to be neighbours; and a third argument, called \code{wpar}, containing 
    any other needed parameter for the weighting function. The available functions are:

    \code{Wf.sqrd} (default value) consists in taking \eqn{w_{i,j} =
      -0.5*d_{i,j}}{w_{i,j} = -0.5*d_{i,j}}; it does not involve any truncation.

    \code{Wf.RBF} consists in taking \eqn{w_{i,j} =
      exp(-wpar*d_{i,j}^2)}{w_{i,j} = exp(-wpar*d_{i,j}^2)}; it does not involve any 
      truncation.

    \code{Wf.binary}; the spatial weighting matrix is
     simply the connectivity matrix, i.e. a binary matrix indicating presence or 
     absence of a connexion between two sites;

    \code{Wf.PCNM} is \eqn{a_{i,j} = 1 - (d_{i,j} /
      (4*boundaries_2))^2}{a_{i,j} = 1 - (d_{i,j} / (4*boundaries_2))^2}
      where \code{boundaries_2} is the maximum threshold value in the 
      \code{boundaries} range of values. See paragraph \code{boundaries}.
      Distance-based MEM (dbMEM eigenfunctions) were first described by 
      Borcard and Legendre (2002) under the acronym PCNM.

    \code{Wf.Drayf1} is \eqn{a_{i,j} = 1 - (d_{i,j} / d_{max})}{a_{i,j}
      = 1 - (d_{i,j} / d_max)} where \eqn{d_max}{dmax} is the distance
      between the two most distant locations in the set,

    \code{Wf.Drayf2} is \eqn{a_{i,j} = 1 - (d_{i,j} /
      d_{max})^{wpar}}{a_{i,j} = 1 - (d_{i,j} / d_max)^wpar},

    \code{Wf.Drayf3} is \eqn{a_{i,j} = 1 / d_{i,j}^{wpar}}{a_{i,j} = 1 /
     d_{i,j}^wpar}.

    Functions \code{Wf.PCNM}, \code{Wf.Drayf1}, \code{Wf.Drayf2}, and \code{Wf.Drayf3}
    were proposed by Dray et al. (2006). They are also described in Legendre and 
    Legendre (2012, Section 14.2). 
    
    The \code{Wf.sqrd} weighting approach is equivalent to submitting
    a matrix of square-rooted dissimilarities to a principal
    coordinate analysis. That option is not often used in the
    ecological litterature, except for beta diversity analysis of non-Euclidean 
    dissimilarity matrices (Legendre and De C&aacute;ceres 2013), but it is actually 
    equivalent, for evenly spaced transect or surfaces (square or rectangular), to 
    using the basis functions of type II discrete cosine basis transforms. }

   \item{boundaries}{ (optional) A pair of threshold values (minimum and maximum, 
    see examples below) used to construct the connectivity matrix. Locations 
    whose distance to one another are between these values (limits included) are 
    considered to be neighbours (\eqn{b_{i,j}=1}{b_{i,j}=1}), whereas values located below
    the mininum or above the maximum are considered to be unconnected, 
    and given the value (\eqn{b_{i,j}=0}{b_{i,j}=0}). 
    Defaults are \eqn{0}{0} for the minimum value and \code{NA}
    for the maximum. \code{NA} values tells the function to take the
    minimum distance value that allows all locations to form a single cluster in 
    single linkage clustering. That value is obtained
    from \code{\link{hclust}}. Ignored when \code{weighting="Wf.sqrd"}
    or when \code{weighting="Wf.RBF"}. }

  \item{wpar}{ A weighting function argument. }
  \item{select}{ The smallest absolute eigenvalue required for an eigenfunction to
    be considered as a suitable predictive variable in MCA modelling. Default value
    depends on one's particular computer. It is set to the square-root
    of \code{.Machine$double.eps} }
  \item{D}{ A distance matrix. }
  \item{object}{ An \code{\link{eigenmap-class}} object. }
  \item{target}{ A distances matrix among the study sites, including all sites used for
    modelling and for prediction, passed on to \code{\link{eigenmap}}. See example 2. }
}
\details{
  Spatial eigenvector maps are sets of eigenfunctions obtained from the
    locations of the observations in a structuring framework, e.g.,
    space or time, or in a graph representing connected objects like atoms 
    in a molecule or a trophic network, or members of a social network 
    in sociology. It is obtained by eigenvalue
    decomposition of a spatial weighting matrix, computed as described
    in Dray et al. (2006) and Legendre & Legendre (2012, Section
    14.2). That square matrix is Gower-centred before
    eigen-decomposition. The spatial weighting matrix is the Hadamard
    product of a connectivity matrix \eqn{\mathbf{B}}{B} and an edge
    weighting matrix \eqn{\mathbf{A}}{A} as described in Legendre & Legendre 
    (2012, Fig. 14.10). 
    
    The function described herein can handle user-chosen transformation functions 
    containing a truncation argument (\code{boundaries}) to producee the binary 
    connectivity matrix \eqn{\mathbf{B}}{B}. The \code{eigenmap} function provides a 
    default method to estimate the boundary parameters should they be missing. 
    It also offers six different ways of computing the weighting matrix 
    \eqn{\mathbf{A}}{A} using arguments \code{boundaries} and \code{wpar}.

    In the present implementation, matrix
    \eqn{\mathbf{B}}{B} is not obtained using a minimum spanning tree as
    suggested by Dray et al. (2006). Instead, it uses a simpler approach where
    all distances within a user-defined truncation interval indicate neighbouring sites 
    in matrix \eqn{\mathbf{B}}{B}.

    Functions \code{\link{Wf.sqrd}}, \code{\link{Wf.RBF}},
    \code{\link{Wf.binary}}, \code{\link{Wf.PCNM}},
    \code{\link{Wf.Drayf1}}, \code{\link{Wf.Drayf2}}, and
    \code{\link{Wf.Drayf3}} are not intended to be called directly, but
    through \code{\link{eigenmap}} (and within
    \code{\link{eigenmap.score}}). Other functions could be written by users and added 
    to the R-code file \code{\link{eigenmap.R}}.
    
    For function \code{\link{eigenmap.score}}, the distance matrix among sampling
    locations, which includes the sites for prediction, must be of the same type as
    the matrix that has been passed to \code{\link{eigenmap}}. If Cartesian
    coordinates were passed to \code{\link{eigenmap}} for the active sites used for 
    modelling, and a pre-computed distance matrix is now passed to the 
    \code{\link{eigenmap.score}} function for all sites (active plus prediction), that 
    distance matrix must also be Euclidean and computed in the same geographic units as 
    the one passed to \code{\link{eigenmap.score}}.
}
\value{
  Function \code{\link{eigenmap}} returns an \code{\link{eigenmap-class}} object.
    Function \code{\link{eigenmap.score}} returns the eigenvector scores for each
    site in the \code{target} distance matrix file. 
}
\references{
  Borcard, D. and Legendre, P. 2002. All-scale spatial analysis of
  ecological data by means of principal coordinates of neighbour
  matrices. Ecol. Model. 153: 51-68
  
  Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
  comprehensive framework for principal coordinate analysis of neighbor
  matrices (PCNM). Ecol. Modelling 196: 483-493

  Legendre, P. & M. De Cáceres. 2013. Beta diversity as the variance of community 
  data: dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963. 

  Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English
  edition. Elsevier Science B.V., Amsterdam, The Neatherlands.
}
\author{ Guillaume Guénard, Departement de sciences biologiques,
  Universite de Montréal, Montréal, Quebec, Canada. }

\seealso{
  \code{\link{MCA}}
  \code{\link{eigenmap-class}}
}
\examples{
#
### Example 1: A linear transect.
#
data(Salmon)
# Variable "Position" in file Salmon: positions of the 76 study sites in the river. 
# See the documentation of data file Salmon in {codep}.
#
## A warning is issued when no boundaries are provided for a function that requires them 
## Example:
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.binary)
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
## In the following examples, boundaries are provided; they are needed by the functions 
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.binary,boundaries=c(0,20))
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.Drayf1,boundaries=c(0,20))
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.Drayf2,boundaries=c(0,20))
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.Drayf3,boundaries=c(0,20),wpar=2)
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.PCNM,boundaries=c(0,20))
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
## In the following examples, boundaries are not provided because they are not needed 
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.sqrd)
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.RBF,wpar=0.001)
tmp <- # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
### Example 2, with predictor scores
#
# smpl: 8 sites removed from the modelling phase, to be used for prediction
smpl <- c(4,7,10,14,34,56,61,64)  
# map: output of eigenmap for the 68 active sites   # see: summary(map)
map <- eigenmap(x=Salmon[-smpl,"Position"],weighting=Wf.sqrd)
# scr: the spatial eigenvector scores of all 76 study sites (file 76 x 67)
scr <- eigenmap.score(object=map,target=as.matrix(dist(Salmon[,"Position"]))[,-smpl])
# Check: the eigenvector scores of the 68 active sites are the same in 'scr' and 'map$U' 
all(round(scr[-smpl,] - map$U, 10) == 0) 
# Extract the spatial eigenvector scores of the 8 prediction sites (file 8 x 67)
tmp <- scr[smpl,]
#
# Plot the 8 prediction sites along particular eigenvectors, here eigenvector #5. 
# The notation 'L' produces a value with class integer instead of numeric
wh <- 5L   # You can also try with other eigenvectors
plot(map$U[,wh]~Salmon[-smpl,"Position"], ylab = expression(U[5]),
     xlab = "Position along transect")
points(y=scr[smpl,wh],x=Salmon[smpl,"Position"],pch=21,bg="black")
#
# Try with a different weighting function, Wf.binary
map <- eigenmap(x=Salmon[-smpl,"Position"],weighting=Wf.binary,boundaries=c(0,20))
scr <- eigenmap.score(object=map,target=as.matrix(dist(Salmon[,"Position"]))[smpl,-smpl])
#
# Plot the 8 prediction sites along particular eigenvectors, here eigenvector #1. 
wh <- 1L   # You can also try with other eigenvectors
plot(map$U[,wh]~Salmon[-smpl,"Position"], ylab = expression(U[1]),
     xlab = "Position along transect (m)")
points(y=scr[,wh],x=Salmon[smpl,"Position"],pch=21,bg="black")
#
# Try with a different weighting function, Wf.PCNM
map <- eigenmap(x=Salmon[-smpl,"Position"],weighting=Wf.PCNM,boundaries=c(0,100))
scr <- eigenmap.score(object=map,target=as.matrix(dist(Salmon[,"Position"]))[smpl,-smpl])
#
# Plot the 8 prediction sites along particular eigenvectors, here eigenvector #1. 
wh <- 1L   # You can also try with other eigenvectors 
plot(map$U[,wh]~Salmon[-smpl,"Position"], ylab = expression(U[1]),
     xlab = "Position along transect (m)")
points(y=scr[,wh],x=Salmon[smpl,"Position"],pch=21,bg="black")
#
### Example 3: An unevenly sampled map surface
#
data(Mite)
map <- eigenmap(x=as.matrix(mite.geo),weighting=Wf.sqrd)
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file
#
map <- eigenmap(x=as.matrix(mite.geo),weighting=Wf.RBF)
tmp <- map # summary(tmp) # plot(tmp)   # Warning: 'map' may be a large output data file 
#
}
