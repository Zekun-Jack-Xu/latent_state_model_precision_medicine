
rm(list=ls())
library(rarhsmm)

load("stepbd_mi5_SAD.RData")
source("util.R")
##########################################

load("hmmparm4.RData")
natparm <- natparm1
MMM <- 4

#################################################
m=5
rawcmf3$CGI <- ifelse(rawcmf3$CGI<0, 1, rawcmf3$CGI)
rawcmf3$DEPRESD <- rawcmf3$DEPRESD / 100
rawcmf3$LESSINT <- rawcmf3$LESSINT / 100
rawcmf3$MOODELEV <- rawcmf3$MOODELEV / 100

###################################################
#DEALING WITH THE MMM imputed DATA SET
rawcmf1 <- subset(rawcmf3,.imp==MMM )

cmflist <- split(rawcmf1,rawcmf1$STEPID)
for(i in length(cmflist):1){
  if(nrow(cmflist[[i]])==1) cmflist[[i]] <- NULL
}
nsubj <- length(cmflist)   #3009

#####################################


xlist <- lapply(1:nsubj,function(k) data.matrix(cmflist[[k]][,c("tempsumd","tempsumm",
                                                                "DEPRESD","LESSINT",
                                                                "MOODELEV"),
                                                             drop=FALSE]) )

timelist <- lapply(1:nsubj,function(k) c(0,cmflist[[k]][,c("DAYSCONS")])  )

ulist <- lapply(1:nsubj,function(k)
  2-cmflist[[k]][,"tempsumd"]/22-cmflist[[k]][,"tempsumm"]/16 )


#1. antilow,2.antimedian, 3.antihigh,
#4. moodlow, 5. moodmedian, 6. moodhigh,
#7-9: antilow + mood
#10-12: antimedian + mood
#13-15: antihigh + mood
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



foolist <- vector(mode="list",length=nsubj)

for(k in 1:nsubj){
  if(nrow(trtcontrastlist[[k]])==1) foolist[[k]] <- 
      matrix(initiallist[[k]],nrow=1)
  else foolist[[k]] <- ar1.smp.filter1(timelist[[k]], natparm, 
                                       initiallist[[k]],trtcontrastlist[[k]],xlist[[k]],
                                       posterior=FALSE)
}


p1 <- unlist(sapply(1:nsubj,function(k) foolist[[k]][-1,1]))
p2 <- unlist(sapply(1:nsubj,function(k) foolist[[k]][-1,2]))
p3 <- unlist(sapply(1:nsubj,function(k) foolist[[k]][-1,3]))
p4 <- unlist(sapply(1:nsubj,function(k) foolist[[k]][-1,4]))
p5 <- unlist(sapply(1:nsubj,function(k) foolist[[k]][-1,5]))
clinstat <- unlist(sapply(1:nsubj,function(k)cmflist[[k]]$CLINSTAT))
aggregate(p1~clinstat,FUN=mean)
aggregate(p2~clinstat,FUN=mean)
aggregate(p3~clinstat,FUN=mean)
aggregate(p4~clinstat,FUN=mean)
aggregate(p5~clinstat,FUN=mean)



trtlist <- lapply(1:nsubj,function(k) cmflist[[k]]$trt)

#unction to get the data together
augdata <- format5(trtlist, foolist, ulist,timelist,xlist)

augdata$x1 <- augdata$x1 / 22
augdata$x2 <- augdata$x2 / 16
augdata$xnext1 <- augdata$xnext1 / 22
augdata$xnext2 <- augdata$xnext2 / 16


#################################
# integrate age, bipolar type

demographic <- lapply(1:nsubj, function(k) cmflist[[k]][-1,c("AGE","LIFEDX1")])
demographic <- Reduce("rbind", demographic)
augdata <- data.frame(cbind(augdata, demographic))

augdata$AGE <- (augdata$AGE - min(augdata$AGE)) / (max(augdata$AGE)-min(augdata$AGE))
augdata$LIFEDX1 <- ifelse(augdata$LIFEDX1==1, 1, 0)


######################################
#propensity score  
library(nnet)

ppsmod <- multinom(trt~x1+x2+x3+x4+x5+AGE+LIFEDX1,data=augdata)
ppsmat <-  predict(ppsmod, data.frame(
  x1=augdata$x1,x2=augdata$x2,x3=augdata$x3,
  x4=augdata$x4, x5=augdata$x5, AGE=augdata$AGE, 
  LIFEDX1=augdata$LIFEDX1),type="probs")

for(r in 1:nrow(ppsmat)){
  for(c in 1:ncol(ppsmat)){
    if(ppsmat[r,c] < 0.05) {
      ppsmat[r,c] <- 0.05
    }
  }
}


##########3
#interactions
augdata$p1tdx <- augdata$p1t * augdata$LIFEDX1
augdata$p2tdx <- augdata$p2t * augdata$LIFEDX1
augdata$p3tdx <- augdata$p3t * augdata$LIFEDX1
augdata$p4tdx <- augdata$p4t * augdata$LIFEDX1
augdata$p1tp1dx <- augdata$p1tp1 * augdata$LIFEDX1
augdata$p2tp1dx <- augdata$p2tp1 * augdata$LIFEDX1
augdata$p3tp1dx <- augdata$p3tp1 * augdata$LIFEDX1
augdata$p4tp1dx <- augdata$p4tp1 * augdata$LIFEDX1

augdata$x1dx <- augdata$x1 * augdata$LIFEDX1
augdata$x2dx <- augdata$x2 * augdata$LIFEDX1
augdata$x3dx <- augdata$x3 * augdata$LIFEDX1
augdata$x4dx <- augdata$x4 * augdata$LIFEDX1
augdata$x5dx <- augdata$x5 * augdata$LIFEDX1
augdata$xnext1dx <- augdata$xnext1 * augdata$LIFEDX1
augdata$xnext2dx <- augdata$xnext2 * augdata$LIFEDX1
augdata$xnext3dx <- augdata$xnext3 * augdata$LIFEDX1
augdata$xnext4dx <- augdata$xnext4 * augdata$LIFEDX1
augdata$xnext5dx <- augdata$xnext5 * augdata$LIFEDX1


############################################################
#1. expected total discounted(0.95) reward
disfactor <- rep(0.95, nrow(augdata))

valuename1 <- c("p1t", "p2t","p3t","p4t","x1","x2","x3","x4","x5",
                "p1tdx","p2tdx","p3tdx","p4tdx",
                "x1dx","x2dx","x3dx","x4dx","x5dx",
                "AGE","LIFEDX1")
valuename2 <- c("p1tp1", "p2tp1","p3tp1","p4tp1",
                "xnext1","xnext2","xnext3","xnext4","xnext5",
                "p1tp1dx", "p2tp1dx","p3tp1dx","p4tp1dx",
                "xnext1dx","xnext2dx","xnext3dx","xnext4dx","xnext5dx",
                "AGE","LIFEDX1")
policyname <-  c("p1t", "p2t","p3t","p4t","x1","x2","x3","x4","x5",
                 "p1tdx","p2tdx","p3tdx","p4tdx",
                 "x1dx","x2dx","x3dx","x4dx","x5dx",
                 "AGE","LIFEDX1")

#linear 
tuning1=1e-8
tuning2=1e-4

maxit=100
eps=1e-3

valueparm <- rep(0, 21)
policyparm <- rep(0, 21*14)
ktrt=15

ss <- proc.time()
plin_dis095 <- fit.sto.dis(valueparm, policyparm,
                           nsubj, augdata, valuename1, valuename2, policyname,
                           disfactor=disfactor,ppsmat=ppsmat,ktrt=ktrt,
                           tuning1=tuning1, tuning2=tuning2,
                           maxit=maxit,eps=eps, print=TRUE)
proc.time() - ss


save(plin_dis095,file="cvlearnparm4int.RData")

####################################################################
#2. expected average reward
runlength <- rle(augdata$subjid)
subjden <- rep(runlength$lengths, runlength$lengths)
augdata$subjden <- subjden


tuning1=1e-8
tuning2=1e-4

maxit=100
eps=1e-3

augvalueparm <- rep(0, 22)
policyparm <- rep(0, 21*14)
ktrt=15

#policygradient
plin_avg <- fit.sto.avg(augvalueparm, policyparm,
                        nsubj, augdata, valuename1, valuename2, policyname,
                        ppsmat, ktrt,
                        tuning1=tuning1, tuning2=tuning2,
                        maxit=maxit,eps=eps,print=TRUE)

save(plin_dis095,plin_avg,file="cvlearnparm4int.RData")
