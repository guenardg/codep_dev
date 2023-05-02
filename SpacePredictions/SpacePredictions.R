##
###
##
## rm(list=ls())
##
##
fun1 <- function(d, wpar = 2*max(d), ...)
  ifelse(
    d < wpar,
    1 - 3*(d/wpar)^2 + 2*(d/wpar)^3,
    0
  )
##
fun2 <- function(d, wpar = 1, ...)
  exp(-wpar * d^2)
##
fun3 <- function(d, wpar = NA, ...)
  -0.5*d
##
fun4 <- function(d, wpar = max(d), ...)
  ifelse(
    d < wpar,
    1 - d/wpar,
    0
  )
##
fun5 <- function(d, wpar = c(max(d),1), ...)
  ifelse(
    d < wpar[1L],
    1 - (d/wpar[1L])^wpar[2L],
    0
  )
##
fun6 <- function(d, wpar = 1, ...)
  1/(1 + d)^wpar
##
par(mfrow=c(3,2),mar=c(4.1,4.1,0.5,0.5))
d <- seq(0,5,0.001)
plot(x=d, y=fun1(d,1), type="l", ylim=c(0,1), las=1L, xlab="", ylab="")
lines(x=d, y=fun1(d,0.5), col="red")
lines(x=d, y=fun1(d,2), col="blue")
plot(x=d, y=fun2(d), type="l", ylim=c(0,1), las=1L, xlab="", ylab="")
lines(x=d, y=fun2(d,0.5), col="red")
lines(x=d, y=fun2(d,2), col="blue")
plot(x=d, y=fun3(d), type="l", las=1L, xlab="", ylab="")
plot(x=d, y=fun4(d,1), type="l", ylim=c(0,1), las=1L, xlab="", ylab="")
lines(x=d, y=fun4(d,0.5), col="red")
lines(x=d, y=fun4(d,2), col="blue")
plot(x=d, y=fun5(d,c(1,1)), type="l", ylim=c(0,1), las=1L, xlab="", ylab="")
lines(x=d, y=fun5(d,c(2,0.5)), col="red")
lines(x=d, y=fun5(d,c(3,0.5)), col="blue")
plot(x=d, y=fun6(d,1), type="l", ylim=c(0,1), las=1L, xlab="", ylab="")
lines(x=d, y=fun6(d,0.5), col="red")
lines(x=d, y=fun6(d,2), col="blue")
rm(d)
dev.copy2eps(file="Common weighting functions.eps")
dev.off()
##
draw <- function(x, fun, ext, by, ..., d0w0 = FALSE,
                 tol = .Machine$double.eps^0.5) {
  n <- NROW(x)
  d <- dist(x)
  w <- as.matrix(d)
  w[] <- fun(w[], ...)
  if(d0w0)
    diag(w) <- 0
  sw <- sum(w)
  g <- t(t(w - rowMeans(w)) - colMeans(w)) + mean(w)
  eigen(g) -> eig
  eig$vectors <- eig$vectors[,abs(eig$values) > tol, drop=FALSE]
  eig$values <- eig$values[abs(eig$values) > tol]
  I <- (eig$values - 1)*n/(sw - n)
  ##
  seq(min(x) - ext, max(x) + ext, by) -> ss
  t(sapply(ss, function(x, y) abs(x - y), y = x)) -> dd
  ww <- dd
  ww[] <- fun(ww[], ...)
  if(d0w0)
    ww[dd==0] <- 0
  gg <- t(t(ww - rowMeans(ww)) - colMeans(w)) + mean(w)
  gg %*% eig$vectors %*% diag(eig$values^(-1)) -> scr
  ##
  list(x = x, eig = eig, ss = ss, scr = scr, w = w, I = I)
}
##
n <- 3
## x <- round(c(0, cumsum(runif(n - 1L, 0.25, 1.75))),2)
## x <- 0:(n-1)
x <- c(0,0.5,1)
## max(diff(x))
draw(x, fun1, 0.1, 0.01, wpar = 2*max(diff(x))) -> drw
## draw(x, fun2, 0.1, 0.01, wpar=1/(2*max(diff(x)))) -> drw
## draw(x, fun3, 0.1, 0.01) -> drw
##
## draw(x, fun4, 0.1, 0.01, wpar = 2*max(diff(x))) -> drw
## draw(x, fun5, 0.1, 0.01, wpar = c(2*max(diff(x)),2)) -> drw
## draw(x, fun6, 0.1, 0.01, wpar=1) -> drw
## 
## drw$eig$values
##
## par(mfrow=c(1,1))
m <- c(1,2)#,3)#,4)
ylim <- max(abs(drw$eig$vectors[,m]),abs(drw$scr[,m])) * c(-1,1)
par(mfrow=c(4,1), mar=c(4.1,4.1,0.5,0.5))
for(i in m) {
  drw$scr[1,i] -> s
  plot(x = drw$ss, y = s*drw$scr[,i], type = "l", ylim = ylim,
       las = 1L, ylab="", xlab="")
  abline(h = 0, lty = 3L)
  points(x = drw$x, y = s*drw$eig$vectors[,i])
}
##
drw$w -> w
diag(w) <- 0
sum(w) -> sw0
g <- t(t(w - rowMeans(w)) - colMeans(w)) + mean(w)
eigen(g) -> eig
eig$vectors <- eig$vectors[,abs(eig$values) > .Machine$double.eps^0.5]
eig$values <- eig$values[abs(eig$values) > .Machine$double.eps^0.5]
##
cbind(drw$eig$values, drw$eig$values - 1, eig$values)
cbind(Moran = eig$values*n/sw0, recalculated = drw$I)
all((round(eig$vectors - drw$eig$vectors,8) == 0) |
      (round(eig$vectors + drw$eig$vectors,8) == 0))

nrow(w) -> n

colSums(w)
matrix(1/n,n,n) %*% w
w %*% matrix(1/n,n,n)






##
## dd[d==0,]

### Observations:
### Il semble que le fait que la fonction ait une dérivée de 0 (ou d'une petite
### valeur) à d = 0 soit une propriété importante



