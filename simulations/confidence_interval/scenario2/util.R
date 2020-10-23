
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
get_meanse_dis <- function(xi, data, ppsmat, nsubj,
valuename1, valuename2,
policyname, ktrt=3, discount=0.9, type="stochastic"){
    
    p2 <- length(policyname)+1
    policyparmmat <- matrix(xi, nrow=ktrt-1, ncol=p2,byrow=TRUE)
    pivec <- sapply(1:nrow(data), function(k) ppsmat[k, data$trt[k]])
    
    if(type=="stochastic"){
        gvec <- sapply(1:nrow(data), function(k){
            xvec <- c(1,data.matrix(data[k,policyname]))
            thistrt <- data$trt[k]
            predict_trt(xvec, policyparmmat)[thistrt]
        } )
    }else{
        gvec <- sapply(1:nrow(data), function(k){
            xvec <- c(1,data.matrix(data[k,policyname]))
            thistrt <- data$trt[k]
            probs <- predict_trt(xvec, policyparmmat)
            ifelse(which.max(probs)==thistrt, 1, 0)
        } )
    }
    
    y <- data$u
    x1 <- data.matrix(cbind(1,data[,valuename1]))
    x2 <- data.matrix(discount * cbind(1, data[,valuename2]))
    
    dx <- (x1 - x2) * gvec/pivec
    dy <- y * gvec / pivec
    alpha <- solve(t(x1)%*%dx, t(x1)%*%dy)
    
    
    c1list <- lapply(1:nrow(data), function(k){
        gvec[k] / pivec[k] * x1[k,] %*% t(x1[k,] - discount*x2[k,])
    })
    c1 <- Reduce("+", c1list) / nsubj
    
    c2list <- lapply(1:nrow(data), function(k){
        mult <- gvec[k] / pivec[k] * (y[k] + x2[k,]%*%alpha - x1[k,]%*%alpha)
        as.numeric(mult^2) * x1[k,]%*%t(x1[k,])
    })
    c2 <- Reduce("+", c2list) / nsubj
    
    c3list <- lapply(1:nrow(data), function(k){
        result <- matrix(0, nrow(c2), ncol(c2))
        xvec <- c(1,data.matrix(data[k,policyname]))
        #thistrt <- data$trt[k]
        #probs <- predict_trt(xvec, policyparmmat)
        #mult <- ifelse(which.max(probs)==thistrt, 1, 0)
        bvec <- x1[k, 2:5]
        varb <- diag(bvec) - outer(bvec, bvec)
        mult <- gvec[k] / pivec[k] * y[k]
        grad1 <- as.numeric(mult) * bvec %*% t(bvec)
        grad2 <- as.numeric(mult) * bvec %*% t(xi[2:5])
        gbmat <- (grad1 + grad2)%*%varb%*%t(grad1 + grad2)
        result[2:5, 2:5] <- gbmat
        return(result)
    })
    c3 <- Reduce("+", c3list) / nsubj
    inv <- solve(c1)
    c1c2=inv%*%(c2+c3)%*%t(inv)
    pnphis=data.matrix(apply(x1,2,mean))
    
    value_mean <- t(pnphis) %*% alpha
    value_se <- sqrt(t(pnphis) %*% c1c2 %*% pnphis)
    return(list(value_mean=value_mean, value_se=value_se))
}


##############################################
###################################################
get_meanse_avg <- function(xi, data, ppsmat, nsubj,
valuename1, valuename2,
policyname, ktrt=3, type="stochastic"){
    n <- nrow(data)
    p2 <- length(policyname)+1
    policyparmmat <- matrix(xi, nrow=ktrt-1, ncol=p2,byrow=TRUE)
    pivec <- sapply(1:nrow(data), function(k) ppsmat[k, data$trt[k]])
    
    if(type=="stochastic"){
        gvec <- sapply(1:nrow(data), function(k){
            xvec <- c(1,data.matrix(data[k,policyname]))
            thistrt <- data$trt[k]
            probs <- predict_trt(xvec, policyparmmat)
            ifelse(which.max(probs)==thistrt, 1, 0)
        } )
    }else{
        gvec <- sapply(1:nrow(data), function(k){
            xvec <- c(1,data.matrix(data[k,policyname]))
            thistrt <- data$trt[k]
            probs <- predict_trt(xvec, policyparmmat)
            ifelse(which.max(probs)==thistrt, 1, 0)
        } )
    }
    
    y <- data$u
    x1 <- data.matrix(cbind(1, data[,valuename1]))
    x2 <- data.matrix(cbind(0, data[,valuename2]))
    
    dx <- (x1 - x2) * gvec/pivec
    dy <- y * gvec / pivec
    alpha <- solve(t(x1)%*%dx, t(x1)%*%dy)
    
    
    c1list <- lapply(1:nrow(data), function(k){
        gvec[k] / pivec[k] * x1[k,] %*% t(x1[k,] - x2[k,])
    })
    c1 <- Reduce("+", c1list) / nsubj
    
    c2list <- lapply(1:nrow(data), function(k){
        mult <- gvec[k] / pivec[k] * (y[k] + x2[k,]%*%alpha - x1[k,]%*%alpha)
        as.numeric(mult^2) * x1[k,]%*%t(x1[k,])
    })
    c2 <- Reduce("+", c2list) / nsubj
    
    c3list <- lapply(1:nrow(data), function(k){
        result <- matrix(0, nrow(c2), ncol(c2))
        xvec <- c(1,data.matrix(data[k,policyname]))
        #thistrt <- data$trt[k]
        #probs <- predict_trt(xvec, policyparmmat)
        #mult <- ifelse(which.max(probs)==thistrt, 1, 0)
        bvec <- x1[k, 2:5]
        varb <- diag(bvec) - outer(bvec, bvec)
        mult <- gvec[k] / pivec[k] * y[k]
        grad1 <- as.numeric(mult) * bvec %*% t(bvec)
        grad2 <- as.numeric(mult) * bvec %*% t(xi[2:5])
        gbmat <- (grad1 + grad2)%*%varb%*%t(grad1 + grad2)
        result[2:5, 2:5] <- gbmat
        return(result)
    })
    c3 <- Reduce("+", c3list) / nsubj
    
    inv <- solve(c1)
    c1c2=inv%*%c2%*%t(inv)
    c1c3=inv%*%c3%*%t(inv)
    value_mean <- alpha[1]
    value_se <- sqrt(abs(c1c2[1,1])+abs(c1c3[1,1]))
    return(list(value_mean=value_mean, value_se=value_se))
}

 
