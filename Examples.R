##
### codep package examples
##
## rm(list=ls())
##
library(codep)
##
compile <- function() {
  try(dyn.unload("./src/Mahal.so"),silent=TRUE)
  system("R CMD SHLIB ./src/Mahal.c")
  dyn.load("./src/Mahal.so")
  source("./R/Mahal.R")
}
##
compile()
##
x <- matrix(rnorm(30,5,2),10,3)
cov <- matrix(c(1,2,3,2,1.5,3,3,2,0.5),3,3)
## solve(cov)
Mahal(x, cov = cov)
dist(x)
##
