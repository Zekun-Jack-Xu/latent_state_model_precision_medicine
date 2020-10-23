library(rarhsmm)

load("stepbd_mi5_SAD.RData")
source("util.R")

###################################################
#DEALING WITH THE MMM imputed DATA SET
MMM <- 2
 
rawcmf3$CGI <- ifelse(rawcmf3$CGI<0, 1, rawcmf3$CGI)
rawcmf3$DEPRESD <- rawcmf3$DEPRESD / 100
rawcmf3$LESSINT <- rawcmf3$LESSINT / 100
rawcmf3$MOODELEV <- rawcmf3$MOODELEV / 100

rawcmf1 <- subset(rawcmf3,.imp==MMM)
table(rawcmf1$trt)
 

####################################
cmflist <- split(rawcmf1,rawcmf1$STEPID)
for(i in length(cmflist):1){
  if(nrow(cmflist[[i]])==1) cmflist[[i]] <- NULL
}
nsubj <- length(cmflist)  #1073 

#################################################
initiallist <- lapply(1:nsubj, function(k) rep(0.2,5))


xlist <- lapply(1:nsubj,function(k) data.matrix(cmflist[[k]][,c("tempsumd","tempsumm",
                                                                "DEPRESD","LESSINT",
                                                                "MOODELEV"),
                                                             drop=FALSE]) )

timelist <- lapply(1:nsubj,function(k) c(0,cmflist[[k]][,c("DAYSCONS")])  )


#15 possible trt: 
#1-3: anti
#4-6: mood
#7-9: antilow + mood
#10-12: antimedian + mood
#13-15: antihigh + mood
#designmatrix: anti, antilow, antimedian, mood, moodlow, moodmedian
trtcontrastlist <- lapply(1:nsubj,function(k){
  trtnum <- cmflist[[k]]$trt
  result <- rep(NULL, 6)
  for(i in 1:length(trtnum)){
    if(trtnum[i]==1){
      temprow <- c(1,1,0,0,0,0)
    }else if(trtnum[i]==2){
      temprow <- c(1,0,1,0,0,0)
    }else if(trtnum[i]==3){
      temprow <- c(1,0,0,0,0,0)
    }else if(trtnum[i]==4){
      temprow <- c(0,0,0,1,1,0)
    }else if(trtnum[i]==5){
      temprow <- c(0,0,0,1,0,1)
    }else if(trtnum[i]==6){
      temprow <- c(0,0,0,1,0,0)
    }else if(trtnum[i]==7){
      temprow <- c(1,1,0,1,1,0)
    }else if(trtnum[i]==8){
      temprow <- c(1,1,0,1,0,1)
    }else if(trtnum[i]==9){
      temprow <- c(1,1,0,1,0,0)
    }else if(trtnum[i]==10){
      temprow <- c(1,0,1,1,1,0)
    }else if(trtnum[i]==11){
      temprow <- c(1,0,1,1,0,1)
    }else if(trtnum[i]==12){
      temprow <- c(1,0,1,1,0,0)
    }else if(trtnum[i]==13){
      temprow <- c(1,0,0,1,1,0)
    }else if(trtnum[i]==14){
      temprow <- c(1,0,0,1,0,1)
    }else if(trtnum[i]==15){
      temprow <- c(1,0,0,1,0,0)
    }
    result <- rbind(result,temprow)
  }
  return(result)
} )



##########3
m <- 5
scalemat <- matrix(rep(c(-3,0,0,-3,0,0),20), 
                   byrow=TRUE,nrow=20,ncol=6)

mumat <- matrix(c(7,1,0.8,0.8,0.1,
                  3,6,0.1,0.1,0.7,
                  8,5,0.7,0.7,0.5,
                  3,5,0.1,0.1,0.7,
                  2,1,0.1,0.1,0.1),
                byrow=TRUE,nrow=5,ncol=5)

covcube <- array(0, dim=c(5,5,5))
arcube <- array(0, dim=c(5,5,5))
covcube[,,1] <- diag(c(6,1,0.1,0.1,0.1))
covcube[,,2] <- diag(c(7,3,0.1,0.1,0.1))
covcube[,,3] <- diag(c(8,2,0.1,0.1,0.1))
covcube[,,4] <- diag(c(6,2,0.1,0.1,0.1))
covcube[,,5] <- diag(c(2,1,0.1,0.1,0.1))


mod <- list(m=m,scalemat=scalemat,
            mumat=mumat, covcube=covcube, arcube=arcube)

#########################################################
#start fitting

nullid <- -1

time <- proc.time()
fitted <- ar1.smp.fit2(mod,timelist, initiallist, xlist, trtcontrastlist, 
                       nullid, mfactor=0, stepsize=0.01,maxit=30,eps=1e-4,
                       ar_diag=TRUE,cov_diag=TRUE, print=FALSE)
print(proc.time() - time)


print(fitted$nllk)
natparm1 <- fitted$natparm
save(natparm1,initiallist, file="hmmparm2.RData")


time <- proc.time()
fitted <- ar1.smp.fit2(mod,timelist, initiallist, xlist, trtcontrastlist,
nullid, mfactor=0, stepsize=0.01,maxit=30,eps=1e-4,
ar_diag=TRUE,cov_diag=TRUE, print=FALSE)
print(proc.time() - time)


print(fitted$nllk)
natparm2 <- fitted$natparm
save(natparm1,natparm2, initiallist, file="hmmparm2.RData")


time <- proc.time()
fitted <- ar1.smp.fit2(mod,timelist, initiallist, xlist, trtcontrastlist,
nullid, mfactor=0, stepsize=0.01,maxit=30,eps=1e-4,
ar_diag=TRUE,cov_diag=TRUE, print=FALSE)
print(proc.time() - time)


print(fitted$nllk)
natparm3 <- fitted$natparm
save(natparm1,natparm2,natparm3, initiallist, file="hmmparm2.RData")
