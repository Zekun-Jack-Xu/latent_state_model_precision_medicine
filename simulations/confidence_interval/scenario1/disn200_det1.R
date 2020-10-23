library(rarhsmm)
library(nnet)

library(rpart)
require(Matrix)

library(numDeriv)

source("batchsim.R")
source("util.R")

nboot <- 500
nsubj <- 200

m <- 5
maxtime <- 364
maxstage <- 50

widthvec <- NULL
covervec <- NULL
biasvec <- NULL

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
  
  
  nullid <- NULL
  for(j in 1:nsubj){
    if(is.null(nrow(xlist[[j]]))) nullid <- c(nullid,j)
    else if(nrow(xlist[[j]])<=3) nullid <- c(nullid,j)
  }
  
  if(length(nullid)==0) nullid <- -1
  
  #
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
  
  #
  nullid2 <- NULL
  for(j in 1:nsubj){
    if(is.null(nrow(xlist[[j]]))) nullid2 <- c(nullid2,j)
    else if(nrow(xlist[[j]])<=1) nullid2 <- c(nullid2,j)
  }
  
  #
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
  #estimating equations
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
  tuning2=1e-1
  type <- 1
  
  ktrt=3
  #valueparm <- rep(0, length(valuename1)+1)
  #policyparm <- rep(0, (length(valuename1)+1)*(ktrt-1))
  valueparm <- c(1.5, -1.5, -1.5, -1, -1.3, -0.5, -0.4, -0.4)
  policyparm <- c(16, 1, -1, -1, 1.5, 0, -0.2, 0.1,
  16, -1, 1, 1, -1.5, 0, 0.2, 0.1)
  print("fit pomdp lin")
  plin <- try(fit.det.dis(valueparm, policyparm,
  nsubj, augdata1, valuename1, valuename2, policyname,
  disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
  tausq=tausq, kappa=kappa, type=type,
  tuning1=tuning1, tuning2=tuning2,print=FALSE),
  silent=TRUE)
  #print(plin)
  if(class(plin)=="try-error") next
  
  
  
  ########################################
  #optimal deterministic
   
  augdata <- format5(trimtrtlist, problist, trimulist, trimtimelist,
                     trimxlist)
  
  augdata <- transform(augdata, p1tq=p1t^2,p2tq=p2t^2,p3tq=p3t^2,
                       p4tq=p4t^2,p1tp1q=p1tp1^2,p2tp1q=p2tp1^2,
                       p3tp1q=p3tp1^2,p4tp1q=p4tp1^2,x1q=x1^2,
                       x2q=x2^2,x3q=x3^2,xnext1q=xnext1^2,
                       xnext2q=xnext2^2,xnext3q=xnext3^2) 
  
  tuning1=1e-8
  tuning2=4e-2
  valueparm <- rep(0, length(valuename1)+1)
  policyparm <- rep(0, (length(valuename1)+1)*(ktrt-1))
  
  print("fit pomdp lin")
  popt <- try(fit.det.dis(valueparm, policyparm,
                          nsubj, augdata1, valuename1, valuename2, policyname,
                          disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
                          tausq=tausq, kappa=kappa, type=type,
                          tuning1=tuning1, tuning2=tuning2, print=FALSE),
              silent=TRUE)
  #print(plin)
  if(class(popt)=="try-error") next
  
  #############################################################3
  #optimal policy: linear
  trtfun10 <- function(obj){
    vec <- c(1,
             obj$delta[1],obj$delta[2],obj$delta[3],obj$delta[4],
             obj$x[1]/10,obj$x[2]/10,obj$x[3]/10)
    coef1v3 <- popt$policyparm[1:8]
    coef2v3 <- popt$policyparm[9:16]
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
  nsims <- 10
  j <- 1
  templist <- vector(mode="list",length=nsims)
  this10 <- numeric(nsims)
  
  print("sim pomdp linear")
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
  
  value_opt <- mean(this10)
  ############################
  # asym. for policy parameters
  
  time <- proc.time()
  gd <- numDeriv::grad(value_dis, plin$policyparm, valueparm=plin$valueparm,
                       nsubj=nsubj, data=augdata1,
                       valuename1=valuename1, valuename2=valuename2,
                       policyname=policyname, disfactor=disfactor,
                       ppsmat=ppsmat, ktrt=3, type="deterministic")
  proc.time() - time
  time <- proc.time()
  hs <- numDeriv::hessian(value_dis, plin$policyparm, valueparm=plin$valueparm,
                          nsubj=nsubj, data=augdata1,
                          valuename1=valuename1, valuename2=valuename2,
                          policyname=policyname, disfactor=disfactor,
                          ppsmat=ppsmat, ktrt=3, type="deterministic")
  proc.time() - time
   
  inv2 <- try(solve(hs),silent=TRUE)
  if(class(inv2)=="try-error") next
  mid2 <- gd%*%t(gd)
  cov_policy <- inv2%*%mid2%*%inv2/nsubj
  se_policy <- sqrt(diag(cov_policy))
  if(any(is.nan(se_policy))) next
  #sampling neighborhood pointwise
  
  ub <- rep(NA, 50)
  lb <- rep(NA, 50)
  
  
  jj <- 1
  while(jj <= 50){ 
    newparm <- runif(length(plin$policyparm),
                     plin$policyparm - 4,
                     plin$policyparm + 4)
    stat <- nsubj * t(newparm - plin$policyparm)%*%cov_policy%*%(newparm - plin$policyparm)
    #print(stat)
    if(stat > qchisq(0.975,df=16)) next
    
    #
    meanselist <- try(get_meanse_dis(newparm, augdata1, ppsmat, nsubj,
                                     valuename1, valuename2,
                                     policyname, ktrt=3, discount=0.9, type="deterministic"), silent=TRUE)
  
    if(class(meanselist)=="try-error") next
    if(is.na(meanselist$value_se) | meanselist$value_se>50) next
    value_est <- try(value_dis(plin$policyparm, plin$valueparm, nsubj, augdata1,
    valuename1, valuename2,  policyname,
    disfactor,ppsmat,ktrt, type="deterministic"), silent=TRUE)
    ub[jj] <- value_est + meanselist$value_se * qnorm(1-0.05/4)/ sqrt(nsubj)
    lb[jj] <- value_est + meanselist$value_se * qnorm(0.05/4) / sqrt(nsubj)
    jj <- jj + 1
  }
  
  ub <- max(ub, na.rm=TRUE) 
  lb <- min(lb, na.rm=TRUE)
  width <- (ub - lb)/2
  
  cover <- ifelse(ub>value_opt & lb<value_opt, 1, 0)
  if(cover==0){
      if(value_opt < 1 | value_est < 1) next # not converging
  }
  
  covervec <- c(covervec, cover)
  widthvec <- c(widthvec, width)
  cat(value_opt, "," ,width, ",",value_est,",",sum(covervec),"\n")
  
  b <- b+1
  
}


#save coverage, width for value_est
save(list=c("covervec","widthvec"), file="disn200_det1.RData")








