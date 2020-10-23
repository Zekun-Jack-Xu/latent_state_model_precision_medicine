

rm(list=ls())
library(rarhsmm)
library(nnet)

library(rpart)
require(Matrix)
source("batchsim.R")
source("util.R")

nboot <- 500
bmean <- matrix(NA,nrow=nboot,ncol=10)
colnames(bmean) <- c("obs","stopomlin","stopomquad","stomdplin","stomdpquad",
                     "detpomlin","detpomquad","detmdplin","detmdpquad", "opt")

nsubj <- 200
m <- 5
maxtime <- 364
maxstage <- 50

##############################3
b <- 1
while(b<=nboot){
    print(b)
    
    ############################################
    trtfun0 <- function(obj){
        vec <- c(1, obj$x[1],obj$x[2],obj$x[3])
        coef1v3 <- c(-0.2,0.1,-0.1,0.1)
        coef2v3 <- c(-0.2,-0.1,0.1,-0.1)
        lp1v3 <- coef1v3 %*% vec
        lp2v3 <- coef2v3 %*% vec
        p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
        
        trt <- gen_multinomial(1,3,c(p1,p2,p3),c(1,2,3))
        #trt3 as the reference
        if(trt==1) trtcontrast <- c(1,0)
        else if(trt==2) trtcontrast <- c(0,1)
        else trtcontrast <- c(0,0)
        
        return(list(trt=trt,trtcontrast=trtcontrast))
    }
    trtfun <- trtfun0
    
    #simulate one step back utility
    utilfun0 <- function(obj){
        vec <- c(1, obj$lastobj$trt,
        obj$x[1],obj$x[2],obj$x[3],
        obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],obj$delta[5],#10
        obj$lastobj$delta[5])
    
        x1 <- vec[3]
        x2 <- vec[4]
        x3 <- vec[5]
        p5 <- vec[10]
        
        trt <- vec[2] 
        
        u <- 2 - abs(x1)  - abs(x3) 
        #state <- which.max(vec[6:10])
        #u <- ifelse(state==5,1,-1)
        return(u)
    }
    
    
    utilfun <- utilfun0
    
    sim1 <- batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
    
    trtfun, utilfun, only_u=FALSE)
    
    xlist <- sim1$xlist
    #trtlist <- sim1$trtlist
    trtcontrastlist <- lapply(1:nsubj,function(k) cbind(1,sim1$trtcontrastlist[[k]]))
    
    timelist <- sim1$timelist
    statelist <- sim1$statelist
    problist <- sim1$problist
    initiallist <- sim1$initiallist
    statelist <- sim1$statelist
    ulist <- sim1$ulist
    mod <- sim1$mod
    
    
    
    ##########################################
    
    ulist1 <- lapply(1:length(ulist),function(k) ulist[[k]][-1])
 
    usumdis <- try(sapply(1:nsubj,function(k){
        rr <- length(ulist1[[k]])
        coef <- 0.9^seq(0,rr-1,by=1)
        t(coef)%*%ulist1[[k]] }),silent=TRUE)
    if(class(usumdis)=="try-error") next
    obsmean <- mean(usumdis,na.rm=TRUE)
    

    #remove those with only 0 or 1 decision stage
    nullid <- NULL
    for(j in 1:nsubj){
        if(is.null(nrow(xlist[[j]]))) nullid <- c(nullid,j)
        else if(nrow(xlist[[j]])<=3) nullid <- c(nullid,j)
    }
    
    if(length(nullid)==0) nullid <- -1
    
 
    hmmtimelist <- lapply(1:nsubj, function(k){
        if(length(timelist[[k]])==1) return(0)
        else return(timelist[[k]][-1])
    } )
    hmminitiallist <- lapply(1:nsubj,function(k) {
        if(is.null(nrow(problist[[k]]))) return(problist[[k]])
        if(nrow(problist[[k]])<2) return(problist[[k]])
        else return(problist[[k]][2,,drop=FALSE])})
    hmmxlist <- lapply(1:nsubj, function(k) {
        if(is.null(nrow(xlist[[k]]))) return(0)
        else return(xlist[[k]][-1,,drop=FALSE])})
    hmmtrtcontrastlist <- trtcontrastlist
    
    fitted <- try(ar1.smp.fit2(mod,hmmtimelist, hmminitiallist,
    hmmxlist, hmmtrtcontrastlist,
    nullid, mfactor=0, stepsize=1e-4,maxit=4,eps=1e-4,
    ar_diag=TRUE,cov_diag=TRUE,print=FALSE),
    silent=TRUE)
    if(class(fitted)=="try-error") next
    #proc.time() - time
    
    natparm <- fitted$natparm
    
    
    #################################################
    foolist <- vector(mode="list",length=nsubj)
    
    for(k in 1:nsubj){
        if(nrow(trtcontrastlist[[k]])==1) foolist[[k]] <-
        matrix(initiallist[[k]],nrow=1)
        else foolist[[k]] <- ar1.smp.filter1(timelist[[k]], natparm,
        initiallist[[k]],trtcontrastlist[[k]],
        xlist[[k]], posterior=FALSE)
    }
    
    
    ####################################
    trtlist <- sim1$trtlist
 
    nullid2 <- NULL
    for(j in 1:nsubj){
        if(is.null(nrow(xlist[[j]]))) nullid2 <- c(nullid2,j)
        else if(nrow(xlist[[j]])<=1) nullid2 <- c(nullid2,j)
    }
    
    
    trimtrtlist <- cleanup(trtlist,nullid2)
    trimfoolist <- cleanup(foolist,nullid2)
    trimulist <- cleanup(ulist,nullid2)
    trimtimelist <- cleanup(timelist,nullid2)
    trimxlist <- cleanup(xlist,nullid2)
    
    
    augdata <- format5(trimtrtlist, trimfoolist, trimulist, trimtimelist,
    trimxlist)
    
    augdata <- transform(augdata, p1tq=p1t^2,p2tq=p2t^2,p3tq=p3t^2,
    p4tq=p4t^2,p1tp1q=p1tp1^2,p2tp1q=p2tp1^2,
    p3tp1q=p3tp1^2,p4tp1q=p4tp1^2,x1q=x1^2,
    x2q=x2^2,x3q=x3^2,xnext1q=xnext1^2,
    xnext2q=xnext2^2,xnext3q=xnext3^2)
    
    #################################################
    ###################################################################
    #propensity score
    
    ppsmod <- multinom(trt ~ x1+x2+x3, data=augdata, trace=FALSE)
    ppsmat <- predict(ppsmod, augdata[,c("x1","x2","x3")],type="probs")
    
    
    disfactor <- rep(0.9, nrow(augdata))
    
    #######################################
    #gather data
    augdata1 <- transform(augdata, x1=x1/10,x2=x2/10,x3=x3/10,
    xnext1=xnext1/10,xnext2=xnext2/10,xnext3=xnext3/10,
    x1q=x1q/10,x2q=x2q/10,x3q=x3q/10,
    xnext1q=xnext1q/10,xnext2q=xnext2q/10,xnext3q=xnext3q/10)
    nsubjrow <- length(trimulist)
    
    
    valuename1 <- c("p1t", "p2t","p3t","p4t","x1","x2","x3")
    valuename2 <- c("p1tp1", "p2tp1","p3tp1","p4tp1","xnext1","xnext2","xnext3")
    policyname <- valuename1
    
    tausq=0.09
    kappa=c(0,0.25,0.5,0.75,1)
    
    tuning1=1e-4
    tuning2=1e-2
    type <- 1
    
    ktrt=3
    valueparm <- rep(0, length(valuename1)+1)
    policyparm <- rep(0, (length(valuename1)+1)*(ktrt-1))
    
    
    print("fit pomdp lin")
    plin <- try(fit.sto.dis(valueparm, policyparm,
    nsubjrow, augdata1, valuename1, valuename2, policyname,
    disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
    tausq=tausq, kappa=kappa, type=type,
    tuning1=tuning1, tuning2=tuning2, print=FALSE),
    silent=TRUE)
    #print(plin)
    if(class(plin)=="try-error") next
    
    #############################################################3
    #linear stochastic
    trtfun2 <- function(obj){
        vec <- c(1,
        obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],
        obj$x[1]/10,obj$x[2]/10,obj$x[3]/10)
        coef1v3 <- plin$policyparm[1:8]
        coef2v3 <- plin$policyparm[9:16]
        lp1v3 <- coef1v3 %*% vec
        lp2v3 <- coef2v3 %*% vec
        p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
        
        trt <- gen_multinomial(1,3,c(p1,p2,p3),c(1,2,3))
        #trt3 as the reference
        if(trt==1) trtcontrast <- c(1,0)
        else if(trt==2) trtcontrast <- c(0,1)
        else trtcontrast <- c(0,0)
        
        return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun2
  
    nsims <- 1
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this2 <- numeric(nsims)
    
    print("sim pomdp linear")
    while(j<=nsims){
        temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
        trtfun, utilfun,only_u=TRUE),
        silent=TRUE)
        if(class(temp)=="try-error") next
        else {
            temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
            this2list <- try(sapply(1:nsubj,function(k){
                rr <- length(temp[[k]])
                coef <- 0.9^seq(0,rr-1,by=1)
                t(coef)%*%temp[[k]]}),silent=TRUE)
            if(class(this2list)=="try-error") next
            this2[j] <- mean(this2list,na.rm=TRUE)
            j <- j+1
        }
    }
    
    
    #quadratic
    
    valuename1 <- c("p1t", "p2t","p3t","p4t","x1","x2","x3",
    "p1tq","p2tq","p3tq","p4tq","x1q","x2q","x3q")
    valuename2 <- c("p1tp1", "p2tp1","p3tp1","p4tp1","xnext1","xnext2","xnext3",
    "p1tp1q", "p2tp1q","p3tp1q","p4tp1q","xnext1q","xnext2q","xnext3q")
    policyname <- c("p1t", "p2t","p3t","p4t","x1","x2","x3")
    
    valueparm <- rep(0, length(valuename1)+1)
    
    print("fit pomdp quad")
    pquad <- try(fit.sto.dis(valueparm, policyparm,
    nsubjrow, augdata1, valuename1, valuename2, policyname,
    disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
    tausq=tausq, kappa=kappa, type=type,
    tuning1=tuning1, tuning2=tuning2, print=FALSE),
    silent=TRUE)
    if(class(pquad)=="try-error") next
    #####################################
    ####################
    #quadratic
    
    trtfun3 <- function(obj){
        vec <- c(1,
        obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],
        obj$x[1]/10,obj$x[2]/10,obj$x[3]/10)
        coef1v3 <- pquad$policyparm[1:8]
        coef2v3 <- pquad$policyparm[9:16]
        lp1v3 <- coef1v3 %*% vec
        lp2v3 <- coef2v3 %*% vec
        p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
        
        trt <- gen_multinomial(1,3,c(p1,p2,p3),c(1,2,3))
        #trt3 as the reference
        if(trt==1) trtcontrast <- c(1,0)
        else if(trt==2) trtcontrast <- c(0,1)
        else trtcontrast <- c(0,0)
        
        return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun3
    
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this3 <- numeric(nsims)
    
    print("sim pomdp quad")
    while(j<=nsims){

        temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
        trtfun, utilfun, only_u=TRUE),
        silent=TRUE)
        if(class(temp)=="try-error") next
        else {
            temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
            this3list <- try(sapply(1:nsubj,function(k){
                rr <- length(temp[[k]])
                coef <- 0.9^seq(0,rr-1,by=1)
                t(coef)%*%temp[[k]]}),silent=TRUE)
            if(class(this3list)=="try-error") next
            this3[j] <- mean(this3list,na.rm=TRUE)
            j <- j+1
        }
    }
     
    ######################################
    # mdp
    valuename1 <- c("x1","x2","x3")
    valuename2 <- c("xnext1","xnext2","xnext3")
    policyname <- valuename1
    
    
    tuning1=1e-6
    tuning2=1e-6
    type <- 1
    
    valueparm <- rep(0, length(valuename1)+1)
    policyparm <- rep(0, (length(valuename1)+1)*(ktrt-1))
    ktrt=3
    
    print("fit mdp lin")
    plin_mdp <- try(fit.sto.dis(valueparm, policyparm,
    nsubjrow, augdata1, valuename1, valuename2, policyname,
    disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
    tausq=tausq, kappa=kappa, type=type,
    tuning1=tuning1, tuning2=tuning2, print=FALSE),
    silent=TRUE)
    if(class(plin_mdp)=="try-error") next
    
    #mdp quadratic
    valuename1 <- c("x1","x2","x3", "x1q","x2q","x3q")
    valuename2 <- c("xnext1","xnext2","xnext3","xnext1q","xnext2q","xnext3q")
    policyname <- c("x1","x2","x3")
    
    valueparm <- rep(0, length(valuename1)+1)
    print("fit mdp quad")
    pquad_mdp <- try(fit.sto.dis(valueparm, policyparm,
    nsubjrow, augdata1, valuename1, valuename2, policyname,
    disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
    tausq=tausq, kappa=kappa, type=type,
    tuning1=tuning1, tuning2=tuning2, print=FALSE),
    silent=TRUE)
    if(class(pquad_mdp)=="try-error") next
    #####################################################

    ###################################################

    trtfun4 <- function(obj){
        vec <- c(1, obj$x[1],obj$x[2],obj$x[3])
        coef1v3 <- plin_mdp$policyparm[1:4]
        coef2v3 <- plin_mdp$policyparm[5:8]
        lp1v3 <- coef1v3 %*% vec
        lp2v3 <- coef2v3 %*% vec
        p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
        
        trt <- gen_multinomial(1,3,c(p1,p2,p3),c(1,2,3))
        #trt3 as the reference
        if(trt==1) trtcontrast <- c(1,0)
        else if(trt==2) trtcontrast <- c(0,1)
        else trtcontrast <- c(0,0)
        
        return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    
    trtfun <- trtfun4
    
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this4 <- numeric(nsims)

    while(j<=nsims){
        temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
        trtfun, utilfun, only_u=TRUE),
        silent=TRUE)
        if(class(temp)=="try-error") next
        else {
            temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
            this4list <- try(sapply(1:nsubj,function(k){
                rr <- length(temp[[k]])
                coef <- 0.9^seq(0,rr-1,by=1)
                t(coef)%*%temp[[k]]}),silent=TRUE)
            if(class(this4list)=="try-error") next
            this4[j] <- mean(this4list,na.rm=TRUE)
            j <- j+1
        }
    }
 
    trtfun5 <- function(obj){
        vec <- c(1, obj$x[1],obj$x[2],obj$x[3])
        coef1v3 <- pquad_mdp$policyparm[1:4]
        coef2v3 <- pquad_mdp$policyparm[5:8]
        lp1v3 <- coef1v3 %*% vec
        lp2v3 <- coef2v3 %*% vec
        p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
        p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
        
        trt <- gen_multinomial(1,3,c(p1,p2,p3),c(1,2,3))
        #trt3 as the reference
        if(trt==1) trtcontrast <- c(1,0)
        else if(trt==2) trtcontrast <- c(0,1)
        else trtcontrast <- c(0,0)
        return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun5
    
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this5 <- numeric(nsims)
    
    while(j<=nsims){

        temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
        trtfun, utilfun, only_u=TRUE),
        silent=TRUE)
        if(class(temp)=="try-error") next
        else {
            temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
            this5list <- try(sapply(1:nsubj,function(k){
                rr <- length(temp[[k]])
                coef <- 0.9^seq(0,rr-1,by=1)
                t(coef)%*%temp[[k]]}),silent=TRUE)
            if(class(this5list)=="try-error") next
            this5[j] <- mean(this5list,na.rm=TRUE)
            j <- j+1
        }
    }
 
    
    #################################################3
    #deterministic
    trtfun6 <- function(obj){
      vec <- c(1,
               obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],
               obj$x[1]/10,obj$x[2]/10,obj$x[3]/10)
      coef1v3 <- plin$policyparm[1:8]
      coef2v3 <- plin$policyparm[9:16]
      lp1v3 <- coef1v3 %*% vec
      lp2v3 <- coef2v3 %*% vec
      p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
      
      trt <- which.max(c(p1,p2,p3))
      #trt3 as the reference
      if(trt==1) trtcontrast <- c(1,0)
      else if(trt==2) trtcontrast <- c(0,1)
      else trtcontrast <- c(0,0)
      
      return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun6
    
    
    nsims <- 1
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this6 <- numeric(nsims)

    while(j<=nsims){
      temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
                                       trtfun, utilfun,only_u=TRUE),
                  silent=TRUE)
      if(class(temp)=="try-error") next
      else {
        temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
        this6list <- try(sapply(1:nsubj,function(k){
          rr <- length(temp[[k]])
          coef <- 0.9^seq(0,rr-1,by=1)
          t(coef)%*%temp[[k]]}),silent=TRUE)
        if(class(this2list)=="try-error") next
        this6[j] <- mean(this6list,na.rm=TRUE)
        j <- j+1
      }
    }
    
    #
    trtfun7 <- function(obj){
      vec <- c(1,
               obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],
               obj$x[1]/10,obj$x[2]/10,obj$x[3]/10)
      coef1v3 <- pquad$policyparm[1:8]
      coef2v3 <- pquad$policyparm[9:16]
      lp1v3 <- coef1v3 %*% vec
      lp2v3 <- coef2v3 %*% vec
      p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
      
      trt <- which.max(c(p1,p2,p3))
      #trt3 as the reference
      if(trt==1) trtcontrast <- c(1,0)
      else if(trt==2) trtcontrast <- c(0,1)
      else trtcontrast <- c(0,0)
      
      return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun7

    j <- 1
    templist <- vector(mode="list",length=nsims)
    this7 <- numeric(nsims)

    while(j<=nsims){
      
      temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
                                       trtfun, utilfun, only_u=TRUE),
                  silent=TRUE)
      if(class(temp)=="try-error") next
      else {
        temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
        this7list <- try(sapply(1:nsubj,function(k){
          rr <- length(temp[[k]])
          coef <- 0.9^seq(0,rr-1,by=1)
          t(coef)%*%temp[[k]]}),silent=TRUE)
        if(class(this7list)=="try-error") next
        this7[j] <- mean(this7list,na.rm=TRUE)
        j <- j+1
      }
    }
    
    #######3
    trtfun8 <- function(obj){
      vec <- c(1, obj$x[1],obj$x[2],obj$x[3])
      coef1v3 <- plin_mdp$policyparm[1:4]
      coef2v3 <- plin_mdp$policyparm[5:8]
      lp1v3 <- coef1v3 %*% vec
      lp2v3 <- coef2v3 %*% vec
      p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
      
      trt <- which.max(c(p1,p2,p3))
      #trt3 as the reference
      if(trt==1) trtcontrast <- c(1,0)
      else if(trt==2) trtcontrast <- c(0,1)
      else trtcontrast <- c(0,0)
      
      return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    
    trtfun <- trtfun8
    
    #linear v-learning
    
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this8 <- numeric(nsims)
    
    while(j<=nsims){
      temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
                                       trtfun, utilfun, only_u=TRUE),
                  silent=TRUE)
      if(class(temp)=="try-error") next
      else {
        temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
        this8list <- try(sapply(1:nsubj,function(k){
          rr <- length(temp[[k]])
          coef <- 0.9^seq(0,rr-1,by=1)
          t(coef)%*%temp[[k]]}),silent=TRUE)
        if(class(this8list)=="try-error") next
        this8[j] <- mean(this8list,na.rm=TRUE)
        j <- j+1
      }
    }
    
    trtfun9 <- function(obj){
      vec <- c(1, obj$x[1],obj$x[2],obj$x[3])
      coef1v3 <- pquad_mdp$policyparm[1:4]
      coef2v3 <- pquad_mdp$policyparm[5:8]
      lp1v3 <- coef1v3 %*% vec
      lp2v3 <- coef2v3 %*% vec
      p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
      
      trt <- which.max(c(p1,p2,p3))
      #trt3 as the reference
      if(trt==1) trtcontrast <- c(1,0)
      else if(trt==2) trtcontrast <- c(0,1)
      else trtcontrast <- c(0,0)
      return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun9
    
    #linear v-learning
    
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this9 <- numeric(nsims)
    
    while(j<=nsims){
      
      temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
                                       trtfun, utilfun, only_u=TRUE),
                  silent=TRUE)
      if(class(temp)=="try-error") next
      else {
        temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
        this9list <- try(sapply(1:nsubj,function(k){
          rr <- length(temp[[k]])
          coef <- 0.9^seq(0,rr-1,by=1)
          t(coef)%*%temp[[k]]}),silent=TRUE)
        if(class(this9list)=="try-error") next
        this9[j] <- mean(this9list,na.rm=TRUE)
        j <- j+1
      }
    }
    
   
    
    
    ########################################
    #optimal deterministic by plugging in the true state probabilities
    trimtrtlist <- cleanup(trtlist,nullid2)
    trimfoolist <- cleanup(problist,nullid2)
    trimulist <- cleanup(ulist,nullid2)
    trimtimelist <- cleanup(timelist,nullid2)
    trimxlist <- cleanup(xlist,nullid2)
    
    
    augdata <- format5(trimtrtlist, trimfoolist, trimulist, trimtimelist,
                       trimxlist)
    
    augdata <- transform(augdata, p1tq=p1t^2,p2tq=p2t^2,p3tq=p3t^2,
                         p4tq=p4t^2,p1tp1q=p1tp1^2,p2tp1q=p2tp1^2,
                         p3tp1q=p3tp1^2,p4tp1q=p4tp1^2,x1q=x1^2,
                         x2q=x2^2,x3q=x3^2,xnext1q=xnext1^2,
                         xnext2q=xnext2^2,xnext3q=xnext3^2)
  
    valuename1 <- c("p1t", "p2t","p3t","p4t","x1","x2","x3")
    valuename2 <- c("p1tp1", "p2tp1","p3tp1","p4tp1","xnext1","xnext2","xnext3")
    policyname <- valuename1
    
    tausq=0.09
    kappa=c(0,0.25,0.5,0.75,1)
    
    tuning1=1e-4
    tuning2=1e-2
    type <- 1
    
    ktrt=3
    valueparm <- rep(0, length(valuename1)+1)
    policyparm <- rep(0, (length(valuename1)+1)*(ktrt-1))
    
    
    print("fit pomdp lin")
    plin <- try(fit.sto.dis(valueparm, policyparm,
                            nsubjrow, augdata1, valuename1, valuename2, policyname,
                            disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
                            tausq=tausq, kappa=kappa, type=type,
                            tuning1=tuning1, tuning2=tuning2, print=FALSE),
                silent=TRUE)
    #print(plin)
    if(class(plin)=="try-error") next
    
    #############################################################3
    
    trtfun10 <- function(obj){
      vec <- c(1,
               obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],
               obj$x[1]/10,obj$x[2]/10,obj$x[3]/10)
      coef1v3 <- plin$policyparm[1:8]
      coef2v3 <- plin$policyparm[9:16]
      lp1v3 <- coef1v3 %*% vec
      lp2v3 <- coef2v3 %*% vec
      p1 <- exp(lp1v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p2 <- exp(lp2v3) / (1 + exp(lp1v3) + exp(lp2v3))
      p3 <- 1 / (1 + exp(lp1v3) + exp(lp2v3))
      
      #trt <- gen_multinomial(1,3,c(p1,p2,p3),c(1,2,3))
      trt <- which.max(c(p1, p2, p3))
      #trt3 as the reference
      if(trt==1) trtcontrast <- c(1,0)
      else if(trt==2) trtcontrast <- c(0,1)
      else trtcontrast <- c(0,0)
      
      return(list(trt=trt,trtcontrast=trtcontrast))
    }
    
    trtfun <- trtfun10
    #linear v-learning
    nsims <- 1
    j <- 1
    templist <- vector(mode="list",length=nsims)
    this10 <- numeric(nsims)
    
    while(j<=nsims){
      temp <- try(batchsim5_fuzzy_inf( m,nsubj,maxtime,maxstage,
                                       trtfun, utilfun,only_u=TRUE),
                  silent=TRUE)
      if(class(temp)=="try-error") next
      else {
        temp <- lapply(1:length(temp),function(k) temp[[k]][-1])
        this10list <- try(sapply(1:nsubj,function(k){
          rr <- length(temp[[k]])
          coef <- 0.9^seq(0,rr-1,by=1)
          t(coef)%*%temp[[k]]}),silent=TRUE)
        if(class(this10list)=="try-error") next
        this10[j] <- mean(this10list,na.rm=TRUE)
        j <- j+1
      }
    }
    
 
    bmean[b,] <- c(obsmean, mean(this2),mean(this3),
    mean(this4),mean(this5),mean(this6),mean(this7),
    mean(this8),mean(this9),mean(this10))
    if(b<30) print(bmean[b,])
    b <- b+1
    
}

save(list=c("bmean"), file="dis09n200_1.RData")








