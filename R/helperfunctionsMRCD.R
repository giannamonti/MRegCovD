############################################################################################################################################################
## Helper functions for MRCD
############################################################################################################################################################

scfactor <- function(alpha, ip)
  # compute the scaling factor
  # input:
  #   alpha: the proportion of the contamination (between 0.5 and 1)
  #   ip: the dimension
  # output:
  #   the scaling factor
{
  return( alpha / pchisq( qchisq(alpha, df=ip) ,df=ip+2) )
}

condnumber <- function(X)
  # compute the condition number(i.e., the ratio between the largest and smallest eigenvalues) of a matrix
  # input:
  #   X: a matrix
  # output:
  #   the condition number
{
  temp=eigen(X)$values
  condnr=max(temp)/min(temp)
}

############################################################################################################################################################
############################################################################################################################################################

TargetCov <- function(mX,target=1)
  # return target covariance matrix with given structure
  # input:
  #   mX: the p by T matrix of the data
  #   target: structure of the robust positive definite target matrix (default=1)
  #           0: target matrix is diagonal matrix with robustly estimated univariate scales on the diagonal  
  #           1: non-diagonal target matrix that incorporates an equicorrelation structure (see (17) in paper) 
  #           2: spearman or kendall rank correlation
  #           3: non-diagonal target matrix with equicorrelation structure using Qn estimator
  # output:
  #   the target matrix
{
  ip = dim(mX)[1]
  vD = apply(mX,1,scaleTau2)
  vD[vD==0] = min(vD[vD!=0])#0 is not allowed on diagonal, replace by small value (0.1?)
  #vD[vD<0.1]=1
  mD = diag(vD, nrow=ip)
  I = diag(1,ip)
  if (target==0){
    R=I
  } else if (target==1){
    cortmp = cor.fk(t(mX)) 
    cortmp=sin(1/2*pi*cortmp) 
    constcor = mean(cortmp[upper.tri(cortmp,diag=FALSE)])
    # KB: add bound to ensure positive definiteness; see paper page 7, below (17)
    if(constcor<= min(c(0,(-1/(ip-1)+0.01))) ){  constcor =  min(c(0,(-1/(ip-1)+0.01)))  }
    # KB: code below can be removed? Also seems not correct as it uses p instead of ip
    # if(mindet>0){
    #   detR2 = function(c){
    #     out = ((1-c)^(p-1))*(1+(p-1)*c)
    #     print(c("det in mindet", out))
    #     return(out-mindet)
    #   }
    #   maxc = uniroot(f=detR2,lower=0,upper=1)$root
    #   print(c("constcor",constcor,"maxc", maxc))
    #   constcor = min( c(constcor,maxc))
    # }
  } else if (target==2){
    Rq  <- Qncormat(t(mX))
    constcor  = median(Rq[upper.tri(Rq,diag=FALSE)])
    J = matrix(1,ip,ip)
    R = constcor*J+(1-constcor)*I
  } else if(target==3){
    #R=cor(t(mX), method = "spearmann") 
    #R = cor(t(mX), method = "kendall")#cor.fk(t(mX)) #
    maxcond=1000
    R=cor(t(mX), method = "spearman") 
    rrr=eigen(R,symmetric=TRUE)
    #R = vectors%*%diag(Lambda)%*%t(vectors)
    Lambda=pmax(rrr$values,rrr$values[1]/maxcond)
    R = rrr$vectors%*%diag(Lambda)%*%t(rrr$vectors)
    #R = cor(t(mX), method = "kendall")#cor.fk(t(mX)) #
  }
  equiCor.mat = mD%*%R%*%mD
  return(equiCor.mat)
}

TargetCorr <- function(mX,target=1,mindet=0)
  # return target correlation matrix with given structure
  # input:
  #   mX: the p by T matrix of the data
  #   target: structure of the robust positive definite target matrix (default=1)
  #           0: identity matrix  
  #           1: non-diagonal matrix with an equicorrelation structure (see (17) in paper) 
  # output:
  #   the target correlation matrix
{
  ip = dim(mX)[1]
  I = diag(1,ip)
  if (target==0){
    R=I
  } 
  else if (target==1){
    cortmp = cor.fk(t(mX)) 
    cortmp=sin(1/2*pi*cortmp) 
    constcor = mean(cortmp[upper.tri(cortmp,diag=FALSE)])
    # KB: add bound to ensure positive definiteness; see paper page 7, below (17)
    if(constcor<= min(c(0,(-1/(ip-1)+0.01))) ){  constcor =  min(c(0,(-1/(ip-1)+0.01)))  }
    # KB: code below can be removed? Also seems not correct as it uses p instead of ip
    # if(mindet>0){
    #   detR2 = function(c){
    #     out = ((1-c)^(p-1))*(1+(p-1)*c)
    #     print(c("det in mindet", out))
    #     return(out-mindet)
    #   }
    #   maxc = uniroot(f=detR2,lower=0,upper=1)$root
    #   print(c("constcor",constcor,"maxc", maxc))
    #   constcor = min( c(constcor,maxc))
    # }
    J = matrix(1,ip,ip)
    R = constcor*J+(1-constcor)*I
  } 
  else if (target==2){
    Rq  <- Qncormat(t(mX))
    constcor  = median(Rq[upper.tri(Rq,diag=FALSE)])
    # KB: add bound to ensure positive definiteness; see paper page 7, below (17)
    if(constcor<= min(c(0,(-1/(ip-1)+0.01))) ){  constcor =  min(c(0,(-1/(ip-1)+0.01)))  }
    J = matrix(1,ip,ip)
    R = constcor*J+(1-constcor)*I
  } 
  else if(target==3){
    maxcond=1000
    R=cor(t(mX), method = "spearman") 
    rrr=eigen(R,symmetric=TRUE)
    #R = vectors%*%diag(Lambda)%*%t(vectors)
    # Clean the smallest eigenvalues; cfr Molenberghs and Rousseeuw, but then with the objective of achieving max cond number
    Lambda=pmax(rrr$values,rrr$values[1]/maxcond)
    R = rrr$vectors%*%diag(Lambda)%*%t(rrr$vectors)
    #R = cor(t(mX), method = "kendall")#cor.fk(t(mX)) #
  }
  equiCor.mat = R
  return(equiCor.mat)
}

RCOV <- function(XX, vMu, rho=NULL, mT, scfac, bcd=NULL, target=1, invert=F)
  # compute the regularized covariance
  # input:
  #   XX, the p by T (or h) matrix of the data (not necessarily demeaned)
  #   vMu: the initial mean (p-vector)
  #   rho: the regularization parameter
  #   mT: the target matrix
  #   scfac: the scaling factor
  #   bcd: diagonal matrix used for rescaling (ratio of scales of target and mcd component)
  #   target: structure of the robust positive definite target matrix (default=1)
  #           0: identity matrix  
  #           1: non-diagonal matrix with an equicorrelation structure (see (17) in paper)
#   invert: if true, gives also inverted regularized covariance 
# output (a list):
#   rho: the regularization parameter
#   mT: the target matrix
#   cov: the covariance matrix based on subset (without regularization step)
#   rcov: the regularized covariance matrix
#   inv_rcov: the inverse of the regularized covariance matrix (if invert=True)
{
  mE = XX-vMu
  iT = dim(mE)[2]; ip = dim(mE)[1]
  if (!is.null(bcd)){mE = sqrt(bcd)%*%mE}
  mS = mE%*%t(mE)/iT
  rcov = rho*mT + (1-rho)*scfac*mS
  if(invert){
    if(ip>iT){
      nu = (1-rho)*scfac; mU = mE/sqrt(iT) 
      #if(target<2){
      inv_rcov = InvSMW(rho=rho, mT=mT, nu=nu, mU=mU)
      #} else {
      #  mB = rho*mT;
      #  inv_rcov=InvSMWg(nu=nu, mU=mU, mB=mB)
      #}
    }else{
      inv_rcov = chol2inv(chol(rcov))
    }
    return(list(rho=rho, mT=mT, cov=mS, rcov=rcov, inv_rcov=inv_rcov))
  }
  else return(list(rho=rho, mT=mT, cov=mS, rcov=rcov))
  
}

InvSMW <- function(rho, mT, nu, mU)
  #compute inverse of covariance matrix using Sherman-Morrison-Woodbury identity when dimension is larger than sample size 
  # input:
  #   rho: the regularization parameter
  #   mT: the target matrix
  #	  nu: the scaling factor multiplied with (1-rho)
  #   mU: the scaled data
  # output:
  #	  the inverse of the covariance matrix
{
  ip = dim(mT)[1]
  pp = dim(mU)
  vD = sqrt(diag(mT))
  imD = diag(vD^(-1))
  R = imD%*%mT%*%imD
  constcor = R[2,1]
  I = diag(1,ip)
  J = matrix(1,ip,ip)
  imR = 1/(1-constcor)*( I - constcor/(1+(ip-1)*constcor)*J )
  imB = (rho)^(-1)*imD%*%imR%*%imD
  Temp <- base::chol2inv( base::chol(diag(pp[2]) + nu * (t(mU) %*% (imB%*%mU))) )
  return( imB - (imB%*%mU) %*% (nu * Temp) %*% (t(mU)%*%imB) )
}

cstep_mrcd <- function(mX, rho=NULL, mT=NULL, target=1, vMu=NULL, mIS=NULL, h=NULL, alpha=.75, index=NULL, maxit=50)
  # Apply generalized C-steps to obtain optimal subset
  # input:
  #   mX: the p by T matrix of the residuals or data, not necessarily demeaned
  #   rho: the regularization parameter
  #   mT: the target matrix
  #   target: structure of the robust positive definite target matrix (default=1)
  #           0: identity matrix  
  #           1: non-diagonal matrix with an equicorrelation structure (see (17) in paper)
  #   vMu: the initial mean (as vector)
  #   mIS: the p by p matrix of the initial inverted covariance
  #   h: the size of the subset OR alpha: the proportion of the contamination
#   maxit: the maximal number of iteration of the C-step algorithm
#   index: the initial subset H_0
# output (a list)
#   index: the optimal h-subset
#   numit: the number of iterations
#   mu: the vector with means
#   cov: the regularized covariance estimate
#   icov: the inverse of the regularized covariance matrix
#   rho: the regularization parameter
#   mT: the target matrix
#   dist: the Mahalanobis distances using the MRCD estimates
#   scfac: the scaling factor
{
  iT = dim(mX)[2]; ip = dim(mX)[1]
  if(is.null(h)) h = floor(alpha*iT)
  
  # scaling factor
  scfac = scfactor(alpha=alpha, ip=ip)
  
  # random choice
  if(is.null(index)) index = sample(1:iT,h) #if no index is given we sample one...
  XX = mX[,index] # p x h
  
  if(is.null(vMu)) vMu = rowMeans(XX)
  
  if(is.null(mIS)){
    ret = RCOV(XX=XX, vMu=vMu, rho=rho, mT=mT, scfac=scfac, target=target, invert=T)
    mIS = ret$inv_rcov
  }
  
  vdst = diag(t(mX-vMu) %*% (mIS %*% (mX-vMu))) #vdst = apply(mX-vMu,2,ftmp)
  index = sort(sort.int(vdst,index.return=T)$ix[1:h])
  
  iter = 1
  while(iter<maxit){
    XX = mX[,index]; vMu = rowMeans(XX)
    ret = RCOV(XX=XX, vMu=vMu, rho=rho, mT=mT, target=target, scfac=scfac, invert=T)
    mIS = ret$inv_rcov
    
    vdst = diag(t(mX-vMu) %*% (mIS %*% (mX-vMu))) 
    nndex = sort(sort.int(vdst,index.return=T)$ix[1:h])
    if(all(nndex==index)) break
    
    index = nndex
    iter = iter+1
  }
  
  return(list(index=index, numit=iter, mu=vMu, cov=ret$rcov,
              icov=ret$inv_rcov, rho=ret$rho, mT=ret$mT, dist=vdst, scfac=scfac))
}


###############################################################################################################
# Qn based correlation estimation

Qncor = function(x,y){
  out = (Qn(x+y)-Qn(x-y))/4
}

Qncormat = function(X){
  p = ncol(X)
  ret <- diag(p)
  for (i in 1:(p-1) ) {
    for (j in (i + 1):p) ret[i, j] <- ret[j, i] <- Qncor(X[,i],X[,j]) 
  }
  return(ret)
}


##################################################################################################
# r6pack

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
  sort.list(mahalanobisD(centeredx, FALSE, lambda))#[1:h]
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