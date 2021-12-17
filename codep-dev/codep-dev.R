##
### Development script.
##
## rm(list=ls())
compile <- function() {
  try(dyn.unload("../codep/src/codep.so"),silent=TRUE)
  system("R CMD SHLIB -o ../codep/src/codep.so ../codep/src/*.c")
  dyn.load("../codep/src/codep.so")
  for(i in list.files("../codep/R","*.R"))
    source(file.path("../codep/R",i))
}
## compile()
library(codep)
##
### A set of reference points:
x <- cbind(c(1,4,5,2,8,4),c(3,6,7,1,3,2))
dimnames(x) <- list(LETTERS[1:6],c("x","y"))
##
### The pairwise Euclidean distances among the reference points: 
d1 <- Euclid(x)
d1
##
### That result is the same as that obtained from function dist:
d2 <- dist(x, method = "euclidean")
all(d1 == d2)
##
### A second set of points:
y <- cbind(c(3,5,7),c(3,6,8))
dimnames(y) <- list(LETTERS[7:9],c("x","y"))
##
### The distances between each point of y (rows) and each point of x (columns):
Euclid(x, y)
##
load("../codep/data/Doubs.rda")
load("../codep/data/salmon.rda")
load("../codep/data/mite.rda")




















##
### Eventual update to the handling of the weighting functions:
##
dist.transform <-
  function(
    type=c("sqrd","RBF","PCNM","binary","Drayf1","Drayf2","Drayf3"),
    boundaries,
    wpar
  ) {
  type <- match.arg(type)
  if(type %in% c("PCNM","binary","Drayf1","Drayf2","Drayf3")) {
    if(!missing(boundaries)) {
      binary <- function(d) {
        b <- d
        b[] <- 0
        b[d > boundaries[1L] & d <= boundaries[2L]] <- 1
        attr(b,"method") <- sprintf("binary(%s)",attr(d,"method"))
        attr(b,"call") <- c(attr(d,"call"),match.call())
        return(b)
      }
    } else stop("Missing argument 'boundaries' for type (",type,")")
  }
  if(missing(wpar))
    if(type %in% c("RBF","Drayf2","Drayf3")) {
      wpar <- 1
    } else if(type=="PCNM") wpar <- 4
  return(
    switch(
      type,
      sqrd = function(d) {
        d[] = -0.5 * d
        attr(d,"method") <- sprintf("sqrd(%s)",attr(d,"method"))
        attr(d,"call") <- c(attr(d,"call"),match.call())
        attr(d,"class") <- c("eigwts",attr(d,"class"))
        d
      },
      RBF = function(d) {
        w <- exp(-wpar * d**2)
        w[d==0] <- 0
        attr(w,"method") <- sprintf("RBF(%s)",attr(w,"method"))
        attr(w,"call") <- c(attr(d,"call"),match.call())
        attr(d,"class") <- c("eigwts",attr(d,"class"))
        w
      },
      PCNM = function(d) {
        d[!(d > boundaries[1L] & d <= boundaries[2L])] <- wpar * boundaries[2L]
        d <- -0.5 * d**2
        attr(d,"method") <- sprintf("PCNM(%s)",attr(d,"method"))
        attr(d,"call") <- c(attr(d,"call"),match.call())
        attr(d,"class") <- c("eigwts",attr(d,"class"))
        d
      },
      binary = binary,
      Drayf1 = function(d) {
        wh0 <- which(d==0)
        d[] <- binary(d) * (1 - d / max(d))
        d[wh0] <- 0
        attr(d,"method") <- sprintf("Drayf1(%s)",attr(d,"method"))
        attr(d,"call") <- c(attr(d,"call"),match.call())
        attr(d,"class") <- c("eigwts",attr(d,"class"))
        d
      },
      Drayf2 = function(d) {
        wh0 <- which(d==0)
        d[] <- binary(d) * (1 - (d / max(d))**wpar)
        d[wh0] <- 0
        attr(d,"method") <- sprintf("Drayf2(%s)",attr(d,"method"))
        attr(d,"call") <- c(attr(d,"call"),match.call())
        attr(d,"class") <- c("eigwts",attr(d,"class"))
        d
      },
      Drayf3 = function(d) {
        wh0 <- which(d==0)
        d[] <- binary(d) / (d**wpar)
        d[wh0] <- 0
        attr(d,"method") <- sprintf("Drayf3(%s)",attr(d,"method"))
        attr(d,"call") <- c(attr(d,"call"),match.call())
        attr(d,"class") <- c("eigwts",attr(d,"class"))
        d
      }
    )
  )
}
##
dtr <- dist.transform("s")
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
dtr <- dist.transform("R")
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
dtr <- dist.transform("P",boundaries=c(0,2))
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
dtr <- dist.transform("b",boundaries=c(0,2))
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
dtr <- dist.transform("Drayf1",boundaries=c(0,2))
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
dtr <- dist.transform("Drayf2",boundaries=c(0,2),wpar=0.5)
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
dtr <- dist.transform("Drayf3",boundaries=c(0,2),wpar=0.5)
dtr(dist(1:10))
attributes(dtr(dist(1:10)))
##
dtr <- dist.transform("Drayf3",c(0,max(hclust(d1,method="single")$height)),wpar=0.5)
## as.matrix(dtr(d1))
##
y <- dtr(d1)
attributes(y)



c(0,max(hclust(dist(1:10),method="single")$height))





