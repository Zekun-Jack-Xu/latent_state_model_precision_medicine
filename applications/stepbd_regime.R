
rm(list=ls())
library(rarhsmm)
source("util.R")

readdata <- function(cvname, hmmname){
  load(cvname)
  plin_dis09_pol <- plin_dis095$policyparm
  plin_avg_pol <- plin_avg$policyparm
  
  load("stepbd_mi5_SAD.RData")
  load(hmmname)
  natparm <- natparm1
  
  rawcmf3$CGI <- ifelse(rawcmf3$CGI<0, 1, rawcmf3$CGI)
  rawcmf3$DEPRESD <- rawcmf3$DEPRESD / 100
  rawcmf3$LESSINT <- rawcmf3$LESSINT / 100
  rawcmf3$MOODELEV <- rawcmf3$MOODELEV / 100
  
  rawcmf1 <- subset(rawcmf3, .imp==1)
  
  cmflist <- split(rawcmf1,rawcmf1$STEPID)
  for(i in length(cmflist):1){
    if(nrow(cmflist[[i]])==1) cmflist[[i]] <- NULL
  }
  nsubj <- length(cmflist)  
  
  xlist <- lapply(1:nsubj,function(k) data.matrix(cmflist[[k]][,c("tempsumd","tempsumm",
                                                                  "DEPRESD","LESSINT",
                                                                  "MOODELEV"),
                                                               drop=FALSE]) )
  
  timelist <- lapply(1:nsubj,function(k) c(0,cmflist[[k]][,c("DAYSCONS")])  )
  
  ulist <- lapply(1:nsubj,function(k) 
    2- cmflist[[k]][,"tempsumd"]/22- cmflist[[k]][,"tempsumm"]/16)  
  
  
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
  
  trtlist <- lapply(1:nsubj,function(k) cmflist[[k]]$trt)
  augdata <- format5(trtlist, foolist, ulist,timelist,xlist)
  augdata$x1scale <- augdata$x1 / 22
  augdata$x2scale <- augdata$x2 / 16
  augdata$xnext1scale <- augdata$xnext1 / 22
  augdata$xnext2scale <- augdata$xnext2 / 16
  
  ################################################
  # integrate demographics
  demographic <- lapply(1:nsubj, function(k) cmflist[[k]][-1,c("AGE","LIFEDX1")])
  demographic <- Reduce("rbind", demographic)
  augdata <- data.frame(cbind(augdata, demographic))
  
  head(augdata,10)
  
  #normalize the age, dichotomize gender and type
  augdata$AGE <- (augdata$AGE - min(augdata$AGE)) / (max(augdata$AGE)-min(augdata$AGE)) 
  augdata$LIFEDX1 <- ifelse(augdata$LIFEDX1==1, 1, 0)
  
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
  
  return(list(plin_dis09_pol=plin_dis09_pol,plin_avg_pol=plin_avg_pol,
              augdata=augdata))
}


list1 <- readdata("cvlearnparm1int.RData","hmmparm1.RData")  
list2 <- readdata("cvlearnparm2int.RData","hmmparm2.RData")  
list3 <- readdata("cvlearnparm3int.RData","hmmparm3.RData")  
list4 <- readdata("cvlearnparm4int.RData","hmmparm4.RData")  
list5 <- readdata("cvlearnparm5int.RData","hmmparm5.RData")  

temp1 <- rbind(list1$plin_dis09_pol,list2$plin_dis09_pol,list3$plin_dis09_pol,
               list4$plin_dis09_pol,list5$plin_dis09_pol)
plin_dis09_pol <- apply(temp1, 2, mean)

temp2 <- rbind(list1$plin_avg_pol,list2$plin_avg_pol,list3$plin_avg_pol,
               list4$plin_avg_pol,list5$plin_avg_pol)
plin_avg_pol <- apply(temp2, 2, mean)

augdata <- list1$augdata
augdata$p1t <- apply(rbind(list1$augdata$p1t,list2$augdata$p1t,list3$augdata$p1t,
                           list4$augdata$p1t,list5$augdata$p1t),2,mean)
augdata$p2t <- apply(rbind(list1$augdata$p2t,list2$augdata$p2t,list3$augdata$p2t,
                           list4$augdata$p2t,list5$augdata$p2t),2,mean)
augdata$p3t <- apply(rbind(list1$augdata$p3t,list2$augdata$p3t,list3$augdata$p3t,
                           list4$augdata$p3t,list5$augdata$p3t),2,mean)
augdata$p4t <- apply(rbind(list1$augdata$p4t,list2$augdata$p4t,list3$augdata$p4t,
                           list4$augdata$p4t,list5$augdata$p4t),2,mean)


augdata$p1tdx <- apply(rbind(list1$augdata$p1tdx,list2$augdata$p1tdx,list3$augdata$p1tdx,
                           list4$augdata$p1tdx,list5$augdata$p1tdx),2,mean)
augdata$p2tdx <- apply(rbind(list1$augdata$p2tdx,list2$augdata$p2tdx,list3$augdata$p2tdx,
                           list4$augdata$p2tdx,list5$augdata$p2tdx),2,mean)
augdata$p3tdx <- apply(rbind(list1$augdata$p3tdx,list2$augdata$p3tdx,list3$augdata$p3tdx,
                           list4$augdata$p3tdx,list5$augdata$p3tdx),2,mean)
augdata$p4tdx <- apply(rbind(list1$augdata$p4tdx,list2$augdata$p4tdx,list3$augdata$p4tdx,
                           list4$augdata$p4tdx,list5$augdata$p4tdx),2,mean)



labels <- sapply(1:nrow(augdata), function(k) {
  ex <- t(cbind(1, 
                data.matrix(augdata[k,c("p1t", "p2t","p3t","p4t","x1","x2","x3","x4","x5",
                                        "p1tdx","p2tdx","p3tdx","p4tdx",
                                        "x1dx","x2dx","x3dx","x4dx","x5dx",
                                        "AGE","LIFEDX1")]) ) )
  pp <- matrix(plin_dis09_pol,14,21,byrow=TRUE) 
  return(which.max(predict_trt(ex, pp)))
})
augdata$labels <- as.character(labels)


savedata <- augdata[,c("subjid","trt","p1t","p2t","p3t","p4t",
                       "x1","x2","x3","x4","x5","x1scale","x2scale",
                       "AGE","LIFEDX1","labels")]
save(savedata, file="stepbd_regimes_cart_dis_int.RData")


#########################################


################################################
labels <- sapply(1:nrow(augdata), function(k) {
  ex <- t(cbind(1, 
                data.matrix(augdata[k,c("p1t", "p2t","p3t","p4t","x1","x2","x3","x4","x5",
                                        "p1tdx","p2tdx","p3tdx","p4tdx",
                                        "x1dx","x2dx","x3dx","x4dx","x5dx",
                                        "AGE","LIFEDX1")]) ) )
  pp <- matrix(plin_avg_pol,14,21,byrow=TRUE) 
  return(which.max(predict_trt(ex, pp)))
})
augdata$labels <- as.character(labels)


savedata <- augdata[,c("subjid","trt","p1t","p2t","p3t","p4t",
                       "x1","x2","x3","x4","x5","x1scale","x2scale",
                       "AGE","LIFEDX1","labels")]
save(savedata, file="stepbd_regimes_cart_avg_int.RData")

