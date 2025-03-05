

optimalh <- function(mX, seth=seth, target=0, rho=NULL, rescale = TRUE, rescaleSVD=TRUE, bc=FALSE,
                 maxcond=1000, minscale=0.001, mindet=0, objective='geom', maxit=200, seed=NULL,initrho=0.1)
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
  n = iT;  
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
  

 # mindorig = r6pack(x=t(mX), h=h, full.h=T, scaled=F) #initial 6 h-subsets from DetMCD
  
  #mindorig = cbind(mind.1,mind.2,mind.3,mind.4,mind.5,mind.6)
  
  # plot(dist.ogk[mind])
  # length(intersect(mind,contindex))
  ## 3.1 Obtain 6 initial subsets (based on DetMCD)
  #mind = r6pack_new(x=t(mX), h=ceiling(iT/2), full.h=FALSE, scaled=TRUE) #initial 6 h-subsets from DetMCD 

  
  i=1;objh=c();bestinitial=c();rrr={};ccc={};
  for (h in seth){
      
      if(!is.null(h)){alpha = h/iT
      }else if(is.null(h)){h = floor(alpha*iT)}
      
     # mind=mindorig
      
      mind = r6pack(x=t(mX), h=h, full.h=T, scaled=F) #initial 6 h-subsets from DetMCD
      
      mind=mind[1:h,]  #Compute distances based on h_0=n/2 and then select h observations!
      
      scfac=scfactor(alpha=alpha,ip=ip)
      
      ## 3.2 Determine smallest value of rho_i for each subset
      #we bepalen rho dan op subset met slechts 0.5 punten.
      
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
      
      #mind = r6pack_full(x=t(mX), hsets=mind, rho=rho, mT=mT, index=mind[,1], target=target) #based on h_0=n/2 sets, order all observations! 
      
      #mind=mind[1:h,]  #Compute distances based on h_0=n/2 and then select h observations!
      
      
      ## 3.4 Apply C-steps from MRCD paper
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
      
      
      # we don't rescale: therefore its on standardized data
      
      
      objh[i]=objret
      rrr[[i]]=ret
      ccc[[i]]=ret$cov
      bestinitial[i]=best6pack
      i=i+1
  }
  
  
  return(list(objh=objh,
              rrr=rrr,
              ccc=ccc,
              seth=seth,
              bestinitial=bestinitial)) #best initial subset (if all are equal, last one is chosen)
  
}