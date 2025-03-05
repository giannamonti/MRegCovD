r6pack = function (x, h, full.h, scaled = TRUE, scalefn = Qn){
  initset <- function(data, scalefn, P, h) {
    stopifnot(length(d <- dim(data)) == 2, length(h) == 1, 
              h >= 1)
    n <- d[1]
    stopifnot(h <= n)
    lambda <- doScale(data %*% P, center = median, scale = scalefn)$scale
    sqrtcov <- P %*% (lambda * t(P))
    sqrtinvcov <- P %*% (t(P)/lambda)
    estloc <- colMedians(data %*% sqrtinvcov) %*% sqrtcov
    centeredx <- (data - rep(estloc, each = n)) %*% P
    sort.list(mahalanobisD(centeredx, FALSE, lambda))[1:h]
  }
  ogkscatter <- function(Y, scalefn, only.P = TRUE) {
    stopifnot(length(p <- ncol(Y)) == 1, p >= 1)
    U <- diag(p)
    for (i in seq_len(p)[-1L]) {
      sYi <- Y[, i]
      ii <- seq_len(i - 1L)
      for (j in ii) {
        sYj <- Y[, j]
        U[i, j] <- (scalefn(sYi + sYj)^2 - scalefn(sYi - 
                                                     sYj)^2)/4
      }
      U[ii, i] <- U[i, ii]
    }
    P <- eigen(U, symmetric = TRUE)$vectors
    if (only.P) 
      return(P)
    Z <- Y %*% t(P)
    sigz <- apply(Z, 2, scalefn)
    lambda <- diag(sigz^2)
    list(P = P, lambda = lambda)
  }
  stopifnot(length(dx <- dim(x)) == 2)
  n <- dx[1]
  p <- dx[2]
  #scalefn <- robScalefn(scalefn, n)
  if (!scaled) {
    x <- doScale(x, center = median, scale = scalefn)$x
  }
  nsets <- 6 #KB
  hsets <- matrix(integer(), h, nsets)
  y1 <- tanh(x)
  R1 <- cor(y1)
  P <- eigen(R1, symmetric = TRUE)$vectors
  hsets[, 1] <- initset(x, scalefn = scalefn, P = P, h = h)
  R2 <- cor(x, method = "spearman")
  P <- eigen(R2, symmetric = TRUE)$vectors
  hsets[, 2] <- initset(x, scalefn = scalefn, P = P, h = h)
  y3 <- qnorm((apply(x, 2L, rank) - 1/3)/(n + 1/3))
  R3 <- cor(y3, use = "complete.obs")
  P <- eigen(R3, symmetric = TRUE)$vectors
  hsets[, 3] <- initset(x, scalefn = scalefn, P = P, h = h)
  znorm <- sqrt(rowSums(x^2))
  ii <- znorm > .Machine$double.eps
  x.nrmd <- x
  x.nrmd[ii, ] <- x[ii, ]/znorm[ii]
  SCM <- crossprod(x.nrmd)
  P <- eigen(SCM, symmetric = TRUE)$vectors
  hsets[, 4] <- initset(x, scalefn = scalefn, P = P, h = h)
  ind5 <- order(znorm)
  half <- ceiling(n/2)
  Hinit <- ind5[1:half]
  covx <- cov(x[Hinit, , drop = FALSE])
  P <- eigen(covx, symmetric = TRUE)$vectors
  hsets[, 5] <- initset(x, scalefn = scalefn, P = P, h = h)
  hsets[, 5] <- initset(x, scalefn = scalefn, P = P, h = h)
  hsets[,5] <- hsets[,4]
  P <- ogkscatter(x, scalefn, only.P = TRUE)
  hsets[, 6] <- initset(x, scalefn = scalefn, P = P, h = h)
  if (full.h) 
    hsetsN <- matrix(integer(), n, nsets)
  for (k in 1:nsets) {
    xk <- x[hsets[, k], , drop = FALSE]
    svd <- classPC(xk, signflip = FALSE)
   # if (svd$rank < p) 
  #    stop("More than half of the observations lie on a hyperplane.")
    score <- (x - rep(svd$center, each = n)) %*% svd$loadings
    ord <- order(mahalanobisD(score, FALSE, sqrt(abs(svd$eigenvalues))))
    if (full.h) 
      hsetsN[, k] <- ord
    else hsets[, k] <- ord[1:h]
  }
  if (full.h) 
    hsetsN
  else hsets
}