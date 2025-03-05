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

### source("helperfunctionsMRCD.R")

####################################################################################################
##  Function mrcd()
####################################################################################################

mrcd.fixedset <- function(mX, target=0, h=NULL, alpha=.75, rho=NULL, rescale = TRUE, rescaleSVD=TRUE, bc=FALSE,
                 maxcond=1000, minscale=0.001, mindet=0, objective='geom', maxit=200, seed=NULL,initrho=0.1,
                 fixedset)
  # KB changed defaults: target=0 and bc=FALSE
  # Corrected formula inverse
  
  # compute the MRCD covariance (and much more)
  # input:
  #   mX, the p by n matrix of the data (pxn)!!!
  #   target: structure of the target matrix (0,'QnDiag', 1, 'equi')
  #   h: the size of the subset OR alpha: the proportion of the contamination (between 0.5 and 1)
  #   rescale: decomposition of covariance matrix (default=TRUE)
  #   rescaleSVD: singular value decomposition of target (default=TRUE)
  #   bc: standardization to ensure proper correlation matrix (ones on diagonal)
  #   maxcond: maximum condition number allowed
  #   minscale: minimum scale allowed
  #   mindet: minimum determinant allowed for correlation matrix
  #   objective: objective function (determinant or pth root of determinant)
  #   seed: the seed for the random number generator
  #   maxit: maximum number of C-steps

  # output:
  #   a list
{
  
  iT = dim(mX)[2]; ip = dim(mX)[1]
  
  #choose objective function (normally det, but due to small numbers pth root of det might be more interesting)
  if (objective=='det'){
    obj<-function(x){det(x)}
  }else if (objective=='geom'){
    obj<-function(x){det(x)^(1/ip)} #geometric mean of eigenvalues
  }else if (objective=='ldet'){
    obj<-function(x){log(det(x))}
  }
  
  if(!is.null(seed)) set.seed(seed)
  
  ## 1. Decomposition of covariance matrix; compute standardized observations u_i (11) using median and Qn 
  if(rescale){
    vmx = apply(mX,1,"median")
    vsd = apply(mX,1,"Qn")
    print(c("correcting",mean(vsd<minscale),"percent of the sds"))
    vsd[vsd<minscale] <- minscale
    Dx = diag(vsd)
    mU = scale(t(mX),center=vmx,scale=vsd)
    mX = t(mU)
    mT=TargetCorr(mX,target=target,mindet=mindet); #tv: corr en mindet als parameters meegegeven
  }else{
    mT=TargetCov(mX,target=target);
    mU = t(mX)# kb: in rescaleSVD mU is required
  }
  
  ##mX2=mX (zie voor distances later #Tim)
  
  ## 2. Rewrite regularized covariance matrix using singular value decomposition of target; compute observations w_i
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
  
  vMu=rowMeans(mX)
   
  ogk<-covOGK(t(mX), sigmamu=s_Qn)
  
 
  # kb
  
  mind=matrix(fixedset,ncol=1)  #Compute distances based on h_0=n/2 and then select h observations!
  h = length(fixedset)
  mindorig=mind
  
  scfac=scfactor(alpha=alpha,ip=ip)
  
  ## 3.2 Determine smallest value of rho_i for each subset
  #we bepalen rho dan op subset met slechts 0.5 punten.
  k=1
  rho6pack=c()
  condnr=c()
  if(is.null(rho)){
 
      #kb
      mXsubset = mX[ , mind]
      vMusubset=rowMeans(mXsubset)
      mE = mXsubset-vMusubset
      mS = mE%*%t(mE)/h    

      if( all(mT == diag(ip))){
        veigen = eigen( scfac*mS )$values
        e1 = min(veigen)
        ep = max(veigen)
        fncond = function(rho){
          # argument can be a vector
          #rveigen = rho + (1-rho)*veigen
          condnr =  (rho+(1-rho)*ep)/(rho+(1-rho)*e1 )
          return( condnr-maxcond  )
        }        
      }else{
        fncond = function(rho){
          rcov = rho*mT + (1-rho)*scfac*mS
          #condnr = kappa(rcov)
          temp = eigen(rcov)$values
          condnr = max(temp)/min(temp)
          return( condnr-maxcond  )
        }      
      }
      out = try(uniroot(f=fncond,lower=0.00001, upper=0.99))
      if(class(out)!="try-error"){
        rho6pack[k]=out$root
      }else{
        grid=c(0.000001,seq(0.001,0.99,by=0.001),0.999999)    
        if( all(mT == diag(ip))){
          objgrid = abs(fncond(grid))
          irho = min( grid[objgrid==min(objgrid)])          
        }else{
          # less relevant
          objgrid = abs(apply(as.matrix(grid),1,"fncond"))
          irho = min( grid[objgrid==min(objgrid)])           
        }
        rho6pack[k]=irho
 

    }
  
 
    rho=max(rho6pack)
  }
  
  #mind = r6pack_full(x=t(mX), hsets=mind, rho=rho, mT=mT, index=mind[,1], target=target) #based on h_0=n/2 sets, order all observations! 
  
  #mind=mind[1:h,]  #Compute distances based on h_0=n/2 and then select h observations!
  
  
  ## 3.4 Apply C-steps from MRCD paper
  ret = cstep_mrcd(mX=mX, rho=rho, mT=mT, h=h, alpha=alpha, index=mind , maxit=maxit) #maximal iterations for each possible starting subset.
  objret = obj(ret$cov)
  print("which subset",1)
  print(c("obj detMRCD", objret))
  print(c("number of iterations - csteps",ret$numit))
  hindex = ret$index
  best6pack = 1
  print(c("obj detMRCD trials", objret))
 
  print(c("MRCD optimal objective function", objret ))
  print(c("best subset", best6pack))
  c_alpha = ret$scfac #scaling factor  
  # kb: eenvoudiger opgeschreven, notatie volgend van rcov
  XX = mX[,hindex] # p x h
  mE = XX-ret$mu
  weightedScov =  mE%*%t(mE)/h
  #weights = rep(0,iT)
  #weights[hindex] = 1
  #weightedScov =(mX-ret$mu)%*%diag(weights)%*%t(mX-ret$mu)/(sum(weights)) 
  #print(c("obj Scov", obj(weightedScov) ))
  
  
  # Re-estimate target matrix on optimal h-subset
  # if(optarget){
  #   if(rescale){
  #     mT=TargetCorr(mX[,hindex],target=target,mindet=mindet)
  #   }else {
  #     mT=TargetCov(mX[,hindex],target=target);
  #   }
  # }
  
  ## 4. Backtransforming for computing final MRCD mean and covariance matrix
  
  ## 4.1  Standardization to ensure proper correlation matrix (18)
  if(bc){
    vd = diag(mT)/diag(weightedScov)
    D = diag(vd)
    weightedScov = sqrt(D)%*%weightedScov%*%sqrt(D) #S^star
    #print(c("obj weigthedScov", obj(weightedScov) ))# should be product det(D)*det(weightedScov)
    scfac = 1
    
  }else{
    D  =  (c_alpha)*diag(1,ip) 
    scfac = 1
  }
  
  MRCDmu = rowMeans(mX[,hindex])
  MRCDcov = rho*mT+(1-rho)*weightedScov
  
  #plot(sqrt(mahalanobis(t(mX), center=MRCDmu,cov=MRCDcov))) #OK
  
  #condnr berekenen
  temp=eigen(MRCDcov)$values
  condnr = max(temp)/min(temp)
  
  if (ip>iT){
    # TBC
    #mE=mX-vMu;
    if (bc){mE=sqrt(D)%*%mE}
    nu=(1-rho)*scfac
    mU=mE/sqrt(h)
    #mU=mE/sqrt(iT)
    if (target<2){
      iMRCDcov=InvSMW(rho=rho,mT=mT,nu=nu,mU=mU)
    } else{
      mB=rho*mT
      iMRCDcov=InvSMWg(nu=nu,mU=mU,mB=mB)
    }
  }else {
    iMRCDcov=chol2inv(chol(MRCDcov))
  }
  
  #plot(sqrt(mahalanobis(t(mX), center=MRCDmu,cov=iMRCDcov,inverted=TRUE))) # OK
  
  ###info voor ons
  MRCDcovscaled=MRCDcov
  mTscaled=mT
  weightedScovscaled=weightedScov
  iMRCDcovscaled=iMRCDcov
  ###
  
  
  if(rescaleSVD){
    MRCDcov = mQ%*%msqL%*%MRCDcov%*%msqL%*%t(mQ)  
    iMRCDcov = mQ%*%misqL%*%iMRCDcov%*%misqL%*%t(mQ)
    mT = mQ%*%msqL%*%mT%*%msqL%*%t(mQ)
  }
  
  #plot(sqrt(mahalanobis(t(mX2),center=MRCDmu,cov=MRCDcov)))
  #plot(sqrt(mahalanobis(t(mX2),center=MRCDmu,cov=iMRCDcov,intverted=TRUE)))
  
    
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
  
  
  ## Apply reweighting
  
  return(list(mu=MRCDmu, 
              cov=MRCDcov,
              icov=iMRCDcov,
              rho=rho,
              index=hindex, #index of observations in final subset
              dist=dist,
              weights=weights,
              weightedScov=weightedScov,
              mT=mT,
              vT=diag(mT),
              condnr=condnr,
              numit=ret$numit, #number of iterations
              best6pack=best6pack, #best initial subset (if all are equal, last one is chosen)
              c_alpha=c_alpha,
              MRCDcovscaled=MRCDcovscaled, ### info voor ons
              mTscaled=mTscaled,### info voor ons
              weightedScovscaled=weightedScovscaled,### info voor ons
              iMRCDcovscaled=iMRCDcovscaled,###info voor ons
              rho6pack=rho6pack,###info voor ons
              mind=mind,
              mindorig=mindorig))
  
  # Re-estimate target matrix on observations with non-zero weight
  
}