 
########################################################################
batchsim5_fuzzy_inf <- function(m,nsubj,maxtime,maxstage,
#mumat0,covcube0,wvec,
trtfun, utilfun,only_u=FALSE){
    
    trtlist <- vector(mode="list",length=nsubj)
    timelist <- vector(mode="list",length=nsubj)
    xlist <- vector(mode="list",length=nsubj)
    statelist <- vector(mode="list",length=nsubj)
    problist <- vector(mode="list",length=nsubj)
    baselist <- vector(mode="list",length=nsubj)
    initiallist <- vector(mode="list",length=nsubj)
    ulist <- vector(mode="list",length=nsubj)
    trtcontrastlist <- vector(mode="list",length=nsubj)
    
    error <- 0
    i <- 1
    
    while(i <= nsubj){
        if(error > 200){ return("error!")}
        
        #gaussian <- gaussianmix(mumat0,covcube0,wvec)
        #base4y <- gaussian$obs
        #delta <- gaussian$probs
        delta <- rep(1/5,5)
        #uu <- runif(1)
        #delta <- c(uu, 1-uu)
        
        
        scalemat <- matrix(c( runif(1,-7,-6),0,0,  runif(1,-7,-6),0,0,
        runif(1,-7,-6),0,5,  runif(1,-7,-6),5,0,
        
        runif(1,-7,-6),0,0,  runif(1,-7,-6),5,0,
        runif(1,-7,-6),0,0,  runif(1,-7,-6),0,5,
        
        runif(1,-7,-6),0,0,  runif(1,-7,-6),5,0,
        runif(1,-7,-6),0,0,  runif(1,-7,-6),0,5,
        
        runif(1,-7,-6),0,5,  runif(1,-7,-6),0,0,
        runif(1,-7,-6),0,0,  runif(1,-7,-6),5,0,
        
        runif(1,-7,-6),2,2,  runif(1,-7,-6),2,2,
        runif(1,-7,-6),2,2,  runif(1,-7,-6),2,2 ),
        byrow=TRUE, nrow=20)
        
        #mumat <- matrix(c(-4,-3,-2,-1,-2,-3,1,2,3,4,3,2,0,0,0),
        #byrow=TRUE,nrow=5,ncol=3)
        mumat <- matrix(c(2,2,2, 2,1,-2, -2,1,2, -2,-2,-2,0,0,0),
        byrow=TRUE,nrow=5,ncol=3)
        
        covcube <- array(0, dim=c(3,3,5))
        arcube <- array(0, dim=c(3,3,5))
        arcube[,,1] <- diag(0.1,3)
        arcube[,,2] <- diag(0.1,3)
        arcube[,,3] <- diag(-0.1,3)
        arcube[,,4] <- diag(-0.1,3)
        covcube[,,1] <- diag(0.2-0.1,3) + matrix(0.1,3,3)
        covcube[,,2] <- diag(0.2-0.1,3) + matrix(0.1,3,3)
        covcube[,,3] <- diag(0.2+0.1,3) + matrix(-0.1,3,3)
        covcube[,,4] <- diag(0.2+0.1,3) + matrix(-0.1,3,3)
        covcube[,,5] <- diag(0.2,3)
        
        
        mod <- list(m=m,scalemat=scalemat,
        mumat=mumat, covcube=covcube, arcube=arcube)
        u5 <- runif(1, -5, -4)
        u6 <- runif(1, -3, -2)
        zeroparm <- c(u5, 0.1,0.1,-0.1,-0.1)#int, p[1:m-1]   -0.1,0.1
        rateparm <- c(u6, 0.1,0.1,-0.1,-0.1)#             0.1,-0.1
        
        
        result <- ar1.smp.sim3(maxtime,maxstage,mod,delta,zeroparm,rateparm,
        trtfun, utilfun)
        
        
        #gather the results
        timelist[[i]] <- result$times
        trtlist[[i]] <- result$trts
        baselist[[i]] <- NULL
        problist[[i]] <- result$stateprobs
        statelist[[i]] <- result$states
        initiallist[[i]] <- delta
        ulist[[i]] <- result$utility
        if(length(ulist[[i]])<2) next
        trtcontrastlist[[i]] <- result$trtcontrasts
        
        if(i==nsubj & is.null(result$series)){
            xlist[[i+1]] <- 1
            xlist[[i+1]] <- NULL
        }else{
            xlist[[i]] <- result$series
        }
        
        i <- i+1
    }
    
    if(only_u) return(ulist)
    else return(list(xlist=xlist,baselist=baselist,timelist=timelist,
    trtlist=trtlist,statelist=statelist,
    initiallist=initiallist,trtcontrastlist=trtcontrastlist,
    problist=problist,ulist=ulist,mod=mod))
}
