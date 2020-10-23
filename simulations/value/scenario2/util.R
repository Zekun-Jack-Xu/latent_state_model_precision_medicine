
require(Matrix)

cleanup <- function(oldlist, nullid){
  newlist <- vector(mode="list",length=length(oldlist)-length(nullid))
  ii <- 1
  for(jj in 1:length(oldlist)){
    if(!jj%in%nullid){
      newlist[[ii]] <- oldlist[[jj]]
      ii <- ii + 1
    }
  }
  return(newlist)
}
###################################################################
#create formatted data  
format5 <- function(trtlist, foolist,  ulist, timelist,xlist){
  
  nsubj <- length(trtlist) 
  #remove the last treatment assignment since
  #the corresponding utility is not observed
  aug_trtlist <- lapply(1:nsubj, function(k)
                 trtlist[[k]][-length(trtlist[[k]])])
  trtvec <- Reduce("c",aug_trtlist)
  
  timedifflist <- lapply(1:nsubj,function(k)diff(timelist[[k]][-length(timelist[[k]])]))
  timediffvec <- Reduce("c",timedifflist)
  
  terminallist <- lapply(1:nsubj,function(k){
                      l <- length(aug_trtlist[[k]])
                      c(rep(0,l-1),1)   })
  terminalvec <- Reduce("c", terminallist)
  
  aug_probt <- lapply(1:nsubj, function(k)
                foolist[[k]][-c(1,nrow(foolist[[k]])),,drop=FALSE])
  probt <- Reduce("rbind", aug_probt)
  aug_probtp1 <- lapply(1:nsubj, function(k)
               foolist[[k]][-c(1,2),,drop=FALSE])
  probtp1 <- Reduce('rbind',aug_probtp1)
  
  aug_ulist <- lapply(1:nsubj,function(k)
              ulist[[k]][-1])
  uvec <- Reduce("c",aug_ulist)
  
  nrep <- sapply(1:nsubj, function(k) length(aug_ulist[[k]]))
  subjvec <- rep(1:nsubj, nrep)
  
  aug_xlist <- lapply(1:nsubj, function(k)
    rbind(xlist[[k]][-nrow(xlist[[k]]),,drop=FALSE])
    )
  aug_x <- Reduce("rbind",aug_xlist)
  
  aug_xlisttp1 <- lapply(1:nsubj, function(k)
    rbind(xlist[[k]][-1,,drop=FALSE])
  )
  aug_xtp1 <- Reduce("rbind",aug_xlisttp1)
  
  xnames <- paste(c(rep("x",ncol(aug_x)),
                    rep("xnext",ncol(aug_x))),
                        1:ncol(aug_x),sep="")
  
 
    augdat <- data.frame(cbind(subjvec,trtvec,probt,probtp1,
                               uvec, terminalvec, timediffvec,
                               aug_x,aug_xtp1))
    
        colnames(augdat) <- c("subjid","trt","p1t","p2t","p3t","p4t","p5t",
                              "p1tp1","p2tp1","p3tp1","p4tp1","p5tp1",
                              "u", "final",
                              "timediff",xnames)

  return(augdat)
}

predict_trt <- function(xvec, pp){
  logits <- c(pp%*%xvec,0)
  exps <- exp(logits)
  return(exps/sum(exps))
  #logp <- logits - max(logits)
  #return(exp(logp))
}





##############################################
###################################################
#valueparm | policyparm: argmin estimating equations
Lambda_func_tot <- function(valueparm, policyparm, nsubj, data, 
                            valuename1, valuename2,  policyname,
                            termprob,ppsmat,ktrt, tuning=0.001){
  
  
  #g/pi * (U + r*Phi(h^{k+1})'theta - Phi(h^k)'theta ) = 0
  p2 <- length(policyname)+1
  policyparmmat <- matrix(policyparm,nrow=ktrt-1,ncol=p2,byrow=TRUE)
  pivec <- sapply(1:nrow(data), function(k) ppsmat[k,data$trt[k]])
  gvec <- sapply(1:nrow(data), function(k){
    xvec <- c(1,data.matrix(data[k,policyname]))
    thistrt <- data$trt[k]
    predict_trt(xvec,policyparmmat)[thistrt]
  } )  
  
  y <- data$u * gvec/pivec
  x1 <- cbind(1,data[,valuename1])
  x2 <- termprob * cbind(1, data[,valuename2])
  x <- data.matrix((x1 - x2)) * gvec/pivec
  
  penalty <- tuning * Matrix::Diagonal(length(valuename1)+1) 
  result <- solve(t(x)%*%x + penalty, t(x)%*%y)
  return(result)
  
}

##############################
#policyparm | valueparm: argmax empirical value
V_func_tot <- function(policyparm, valueparm,
                       nsubj, data, 
                       valuename1, valuename2,  policyname,
                       termprob,ppsmat,ktrt,
                       tausq=0.0625, 
                       kappa=c(0,0.25,0.5,0.75,1),
                       type=1, tuning=0.001){
  
  
  utility <- data$u
  trt <- data$trt
  #pps <- data$pps
  valuemat1 <- data.matrix(data)[,valuename1,drop=FALSE]
  valuemat2 <- data.matrix(data)[,valuename2,drop=FALSE]
  policymat <- data.matrix(data)[,policyname,drop=FALSE]
  p2 <- length(policyname)+1
  policyparmmat <- matrix(policyparm,nrow=ktrt-1,ncol=p2,byrow=TRUE)
  
  val <- vEEtot(nrow(data), valueparm, valuemat1, valuemat2,
                policyparmmat, policymat,utility, trt,
                ppsmat, termprob, tausq,kappa,type,tuning)
  
  return(-val)
}

#########################################
####################################
fit.tot <- function(valueparm, policyparm,
                    nsubj, augdata, valuename1, valuename2, policyname,
                    termprob, ppsmat,ktrt, tausq=0.0625, 
                    kappa=c(0,0.25,0.5,0.75,1),
                    type=1,  tuning1=1, tuning2=1,
                    maxit=30,eps=1e-3,
                    print=TRUE){
  
  #start with value parm all zeros
  
  temppol <- optim(policyparm, V_func_tot,
                   valueparm=valueparm,
                   nsubj=nsubj,data=augdata, 
                   valuename1=valuename1, valuename2=valuename2,
                   policyname=policyname,
                   termprob=termprob,ppsmat=ppsmat,
                   ktrt=ktrt,tausq=tausq,
                   kappa=kappa, type=type, tuning=tuning2,
                   method="BFGS")
  
  policyparm <- temppol$par
  old_par_pol <- temppol$par
  old_nllk_pol <- temppol$val
  
  for(k in 1:maxit){
    if(print==TRUE) cat("iteration:", k,"\n")
    
    tempv <- Lambda_func_tot(valueparm, policyparm, nsubj, augdata, 
                             valuename1, valuename2,  policyname,
                             termprob, ppsmat, ktrt,tuning=tuning1) 
    
    valueparm <- as.vector(tempv)

    temppol <- optim(policyparm, V_func_tot,
                     valueparm=valueparm,
                     nsubj=nsubj,data=augdata, 
                     valuename1=valuename1, valuename2=valuename2,
                     policyname=policyname,
                     termprob=termprob,ppsmat=ppsmat,
                     ktrt=ktrt,tausq=tausq,
                     kappa=kappa, type=type, tuning=tuning2, 
                     method="BFGS")
    
    policyparm <- temppol$par
    new_par_pol <- temppol$par
    new_nllk_pol <- temppol$value
    
    tol1 <- sum(abs(new_par_pol-old_par_pol))/sum(abs(old_par_pol))
    tol2 <- abs(new_nllk_pol - old_nllk_pol)/abs(old_nllk_pol)
    tol3 <- abs(new_nllk_pol - old_nllk_pol)
    if( tol1 < eps & (tol2 < eps | tol3 < eps)){
      break
    }else{
      old_nllk_pol <- new_nllk_pol
      old_par_pol <- new_par_pol
    }
    if(print==TRUE) {
      cat("reldiff_pol:", tol1,"; -value",old_nllk_pol,"\n")
      print(valueparm)
      print(policyparm)
    }
  }
  
  
  return(list(valueparm=valueparm,policyparm=policyparm))
}

##############################################
###################################################
#valueparm | policyparm: argmin estimating equations
Lambda_func_avg <- function(augvalueparm, policyparm, nsubj, data,
valuename1, valuename2,  policyname,
ppsmat, ktrt, tausq=0.0625,
kappa=c(0,0.25,0.5,0.75,1),
type=1, tuning=0.001){
    
    
    #separate datamat into valuemat and policymat
    #separate parms into valueparm and policyparm
    utility <- data$u
    trt <- data$trt
    #pps <- data$pps
    subjden <- augdata$subjden
    
    valueparm <- augvalueparm[-length(augvalueparm)]
    vavg <- augvalueparm[length(augvalueparm)]
    
    valuemat1 <- data.matrix(data)[,valuename1,drop=FALSE]
    valuemat2 <- data.matrix(data)[,valuename2,drop=FALSE]
    policymat <- data.matrix(data)[,policyname,drop=FALSE]
    p2 <- length(policyname)+1
    policyparmmat <- matrix(policyparm,nrow=ktrt-1,ncol=p2,byrow=TRUE)
    
    val <- lambdaEEavg(nsubj, valueparm, valuemat1, valuemat2,
    policyparmmat, policymat,utility, trt,
    ppsmat, subjden,tausq,kappa,type,tuning, vavg)
    
    return(val)
    
}

##############################
#policyparm | valueparm: argmax empirical value
V_func_avg <- function(policyparm, augvalueparm,
nsubj, data,
valuename1, valuename2,  policyname,
ppsmat,ktrt,tausq=0.0625,
kappa=c(0,0.25,0.5,0.75,1),
type=1, tuning=0.001){
    
    m <- m
    utility <- data$u
    trt <- data$trt
    #pps <- data$pps
    subjden <- augdata$subjden
    
    timediff <- data$timediff
    thisstate <- data$thisstate
    valuemat1 <- data.matrix(data)[,valuename1,drop=FALSE]
    valuemat2 <- data.matrix(data)[,valuename2,drop=FALSE]
    policymat <- data.matrix(data)[,policyname,drop=FALSE]
    valueparm <- augvalueparm[-length(augvalueparm)]
    
    p2 <- length(policyname)+1
    policyparmmat <- matrix(policyparm,nrow=ktrt-1,ncol=p2,byrow=TRUE)
    
    val <- vEEavg(nsubj, valueparm, valuemat1, valuemat2,
    policyparmmat, policymat,utility, trt,
    ppsmat, subjden, tausq,kappa,type,tuning)
    
    return(-val)
}

#########################################
####################################
fit.avg <- function(augvalueparm, policyparm,
nsubj, augdata, valuename1, valuename2, policyname,
ppsmat, ktrt, tausq=0.0625,
kappa=c(0,0.25,0.5,0.75,1),
type=1,  tuning1=1, tuning2=1,
maxit=30,eps=1e-3,
print=TRUE,...){
    
    #start with value parm all zeros
    
    temppol <- optim(policyparm, V_func_avg,
    augvalueparm=augvalueparm,
    nsubj=nsubj,data=augdata,
    valuename1=valuename1, valuename2=valuename2,
    policyname=policyname, ppsmat=ppsmat,
    ktrt=ktrt,tausq=tausq,
    kappa=kappa, type=type, tuning=tuning2,
    method="BFGS",...)
    
    policyparm <- temppol$par
    old_par_pol <- temppol$par
    old_nllk_pol <- temppol$val
    
    for(k in 1:maxit){
        if(print==TRUE) cat("iteration:", k,"\n")
        
        tempv <- optim(augvalueparm, Lambda_func_avg,
        policyparm=policyparm,
        nsubj=nsubj,ppsmat=ppsmat,data=augdata,
        valuename1=valuename1, valuename2=valuename2,
        policyname=policyname, ktrt=ktrt,tausq=tausq,
        kappa=kappa, type=type, tuning=tuning1,
        method="BFGS")
        
        augvalueparm <- tempv$par
        
        temppol <- optim(policyparm, V_func_avg,
        augvalueparm=augvalueparm,
        nsubj=nsubj,data=augdata,
        valuename1=valuename1, valuename2=valuename2,
        policyname=policyname, ppsmat=ppsmat,
        ktrt=ktrt,tausq=tausq,
        kappa=kappa, type=type, tuning=tuning2,
        method="BFGS",...)
        
        policyparm <- temppol$par
        new_par_pol <- temppol$par
        new_nllk_pol <- temppol$value
        
        tol1 <- sum(abs(new_par_pol-old_par_pol))/sum(abs(old_par_pol))
        tol2 <- abs(new_nllk_pol - old_nllk_pol)/abs(old_nllk_pol)
        tol3 <- abs(new_nllk_pol - old_nllk_pol)
        if( tol1 < eps & (tol2 < eps | tol3 < eps)){
            break
        }else{
            old_nllk_pol <- new_nllk_pol
            old_par_pol <- new_par_pol
        }
        if(print==TRUE) {
            cat("reldiff_pol:", tol1,"; -value",-augvalueparm[length(augvalueparm)],"\n")
            print(augvalueparm)
            print(policyparm)
        }
    }
    
    
    return(list(augvalueparm=augvalueparm,policyparm=policyparm))
}
