####################################################################################################
##  name: MRCD
##  Minimum Regularized Covariance Determinant
####################################################################################################

require(robustbase)
require(tawny)
require(Matrix)
require(pcaPP) #Kendall's tau

####################################################################################################
##  Auxiliary functions
####################################################################################################

# source("helperfunctionsMRCD.R")
# source("r6pack.R")
####################################################################################################
##  Function mrcd()
####################################################################################################


#' mrcd ...
#' 
#' @param mX ...
#' @param target ...
#' @param h ...
#' @param alpha ...
#' @param rho ...
#' @param rescale ...
#' @param rescaleSVD ...
#' @param bc ...
#' @param maxcond ...
#' @param minscale ...
#' @param mindet ...
#' @param objective ...
#' @param maxit ...
#' @param initrho ...
#' @export
mrcd <- function(mX, target=0, h=NULL, alpha=.75, rho=NULL, rescale = TRUE, rescaleSVD=TRUE, bc=FALSE,
                 maxcond=1000, minscale=0.001, mindet=0, objective='geom', maxit=200, initrho=0.1)

  # compute the Minimum Regularized Covariance Determinant (MRCD) estimator
  # input:
  #   mX: the p by n or T matrix of the data (note that p is the dimension and n or T the size of the data)
  #   target: structure of the robust positive definite target matrix (default=1)
  #           0: target matrix is diagonal matrix with robustly estimated univariate scales on the diagonal  
  #           1: non-diagonal target matrix that incorporates an equicorrelation structure (see (17) in paper) 
  #           When rescale = TRUE, target matrix is the equicorrelation matrix. 
  #   h: the size of the subset (between ceiling(n/2) and n) OR alpha: the proportion of the contamination (between 0.5 and 1)
  #   rescale: if true, the p variables are first (robustly) standardized (default=TRUE)
  #   rescaleSVD: if true, use the singular value decomposition of the target matrix (default=TRUE)
  #   bc: if true, standardization is used to ensure proper correlation matrix (ones on diagonal). 
#   maxcond: maximum condition number allowed (see step 3.4 in algorithm 1)
#   minscale: minimum scale allowed
#   mindet: minimum determinant allowed for target matrix ()
#   objective: objective function to determine optimal subset, see (3) in paper (default='geom')
#           'det': typically one minimizes the determinant of the sample covariance based on the subset
#           'geom': p-th root of determinant or standardized generalized variance (for numerical reasons)
#   maxit: maximum number of generalized C-steps for each initial subset 

# output (a list):
#   mu: MRCD mean estimate, see (11) in paper
#   cov: MRCD scatter estimate, see (12) in paper
#   icov: MRCD precision matrix (inverse of covariance estimate), see (14) in paper
#   rho: estimated regularization parameter
#   index: optimal h-subset yielding lowest determinant
#   dist: Mahalanobis distances calculated using MRCD mean and precision matrix
#   mT: estimated target matrix

{
  # KB
  iT = dim(mX)[2]; ip = dim(mX)[1]
  
  #choose objective function to determine optimal subset
  if (objective=='det'){
    obj<-function(x){det(x)}
  }else if (objective=='geom'){
    obj<-function(x){det(x)^(1/ip)} #geometric mean of eigenvalues
  }
  
  ## 1. Standardize the p variables: compute standardized observations u_i, see (6) in paper, using median and Qn estimator
  if(rescale){
    vmx = apply(mX,1,"median")
    vsd = apply(mX,1,"Qn")
    # print(c("",mean(vsd<minscale),"percent of the sds"))
    vsd[vsd<minscale] <- minscale
    Dx = diag(vsd)
    mU = scale(t(mX),center=vmx,scale=vsd)
    mX = t(mU)
    mT=TargetCorr(mX,target=target,mindet=mindet);
  }else{
    mT=TargetCov(mX,target=target);
    mU = t(mX)
  }
  
  ## 2. Perform singular value decomposition of target matrix and compute observations w_i
  if(rescaleSVD){
    mTeigen=eigen(mT)
    mQ=mTeigen$vectors
    mL=diag(mTeigen$values)
    msqL=diag(sqrt(mTeigen$values))
    misqL=diag(sqrt(mTeigen$values)^(-1))
    mW=mU%*%mQ%*%misqL
    mX=t(mW)
    mT=diag(ip)
  }
  #vMu=rowMeans(mX)
  
  ## 3.1-3.2 Follow Hubert et al. (2012) to obtain 6 initial scatter matrices (if scatter matrix is not invertible, use its regularized version)
  ## 3.3 Determine subsets with lowest Mahalanobis distance  
  
  mind = r6pack(x=t(mX), h=h, full.h=T, scaled=F) #initial 6 h-subsets from DetMCD
 # temp = kbdetmcd(x=t(mX),h=h)
 # mind = kbdetmcd(x=t(mX),h=h)$mind
  #mind = kbdetmcd(x=t(mX),h=h)$best
  if(!is.null(h)){alpha = h/iT
  }else if(is.null(h)){h = floor(alpha*iT)}
  
  mind=mind[1:h,] 
  scfac=scfactor(alpha=alpha,ip=ip)
  
  ## 3.4 Determine smallest value of rho_i for each subset
  
  rho6pack=c()
  condnr=c()
  nsets <- ncol(mind) 
  if(is.null(rho)){
    for (k in 1:nsets){
      mXsubset = mX[ , mind[,k]]
      vMusubset=rowMeans(mXsubset)
      mE = mXsubset-vMusubset
      mS = mE%*%t(mE)/h    
      
      if( all(mT == diag(ip))){
        veigen = eigen( scfac*mS )$values
        e1 = min(veigen)
        ep = max(veigen)
        fncond = function(rho){
          condnr =  (rho+(1-rho)*ep)/(rho+(1-rho)*e1 )
          return( condnr-maxcond  )
        }        
      }else{
        fncond = function(rho){
          rcov = rho*mT + (1-rho)*scfac*mS
          temp = eigen(rcov)$values
          condnr = max(temp)/min(temp)
          return( condnr-maxcond  )
        }      
      }
      out = try(uniroot(f=fncond,lower=0.00001, upper=0.99), silent=TRUE)
      if(class(out)!="try-error"){
        rho6pack[k]=out$root
      }else{
        grid=c(0.000001,seq(0.001,0.99,by=0.001),0.999999)    
        if( all(mT == diag(ip))){
          objgrid = abs(fncond(grid))
          irho = min( grid[objgrid==min(objgrid)])          
        }else{
          objgrid = abs(apply(as.matrix(grid),1,"fncond"))
          irho = min( grid[objgrid==min(objgrid)])           
        }
        rho6pack[k]=irho
      }
      
    }
    
    ## 3.5 Set rho as max of the rho_i's obtained for each subset in previous step
    cutoffrho= max( c(0.1,median(rho6pack)))
    rho=max(rho6pack[rho6pack<=cutoffrho])
    Vselection=seq(1, nsets )
    Vselection[rho6pack>cutoffrho]=NA
    if (sum(!is.na(Vselection))==0){stop("None of the initial subsets is well-conditioned")}
    initV=min(Vselection, na.rm=TRUE)
    setsV=Vselection[!is.na(Vselection)]
    setsV=setsV[-1]
  }
  
  ## 3.6 For each of the six initial subsets, repeat the generalized C-steps (from Theorem 1) until convergence
  ## 3.7 Choose final subset that has lowest determinant among the ones obtained from the six initial subsets
  ret = cstep_mrcd(mX=mX, rho=rho, mT=mT, h=h, alpha=alpha, index=mind[,initV], maxit=maxit) 
  objret = obj(ret$cov)
  hindex = ret$index
  best6pack = initV
  for(k in setsV){ 
    tmp = cstep_mrcd(mX=mX, rho=rho, mT=mT, target=target, h=h, alpha=alpha, index=mind[,k], maxit=maxit) 
    objtmp = obj(tmp$cov)
    if(objtmp<=objret){ret=tmp; objret=objtmp; hindex=tmp$index; best6pack = k} 
  }
  c_alpha = ret$scfac #scaling factor  
  XX = mX[,hindex]
  mE = XX-ret$mu
  weightedScov =  mE%*%t(mE)/h
  
  ## 4 Standardization to ensure proper correlation matrix, see (13) in paper
  if(bc){
    vd = diag(mT)/diag(weightedScov)
    D = diag(vd)
    weightedScov = sqrt(D)%*%weightedScov%*%sqrt(D)
    scfac = 1
    
  }else{
    D  =  (c_alpha)*diag(1,ip) 
    scfac = 1
  }
  
  # MRCD estimates of the standardized data W (inner part of (12) in paper)
  MRCDmu = rowMeans(mX[,hindex])
  MRCDcov = rho*mT+(1-rho)*weightedScov
  
  # Computing inverse of scaled covariance matrix, using SMW identity if data is fat (inner part of (14) and (15) in paper). 
  if ( (ip>iT) & (target<=1) ){
    # !!!! formula InvSMW  used is when T is equicorrelation
    
    # KB correct
    #mE=mX-vMu;
    if (bc){mE=sqrt(D)%*%mE}
    nu=(1-rho)*scfac
    #mU=mE/sqrt(iT)
    mU=mE/sqrt(h)
    iMRCDcov=InvSMW(rho=rho,mT=mT,nu=nu,mU=mU)
  }else {
    iMRCDcov=chol2inv(chol(MRCDcov))
  }
  
  # /Backtransforming the rescaling steps that we applied in the beginning (outer part of (11) and (12) in paper)
  if(rescaleSVD){
    MRCDcov = mQ%*%msqL%*%MRCDcov%*%msqL%*%t(mQ)  
    iMRCDcov = mQ%*%misqL%*%iMRCDcov%*%misqL%*%t(mQ)
    mT = mQ%*%msqL%*%mT%*%msqL%*%t(mQ)
  }
  
  if(rescale){
    mX = t(t(mX)%*%Dx)+vmx
    MRCDmu = Dx%*%MRCDmu+vmx
    MRCDcov = Dx%*%MRCDcov%*%Dx
    mT = Dx%*%mT%*%Dx
    iDx = diag(1/diag(Dx))
    iMRCDcov = iDx%*%iMRCDcov%*%iDx
  }
  
  
  ## Compute the Mahalanobis distances based on MRCD estimates
  dist = sqrt(mahalanobis(t(mX), center=MRCDmu, cov=iMRCDcov, inverted=TRUE))
  

  return(list(mu=MRCDmu, 
              cov=MRCDcov,
              icov=iMRCDcov,
              rho=rho,
              index=hindex, 
              dist=dist,
              mT=mT,
              mind = mind
  )
  )
  
  
}
