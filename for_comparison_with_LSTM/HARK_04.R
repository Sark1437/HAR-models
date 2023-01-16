# Description: HARK model
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 04 Trimming dataset, QLIKE correction.

#Importing data
library(HARModel); library(xts); data(SP500RM);
setwd("/Users/sarkaajsingh/Desktop/")
snp <- read.csv("snp_daily_returns.csv",sep=" "); 
snp$returns <- (log(snp$GSPC.Close)-log(snp$GSPC.Open))*100
snp$Index <- strptime(as.character(snp$Index), "%Y-%m-%d") #Parsing the column as a date
snp <- xts(snp$returns,order.by = snp[,1]) #Converting df into ts
n <- intersect(as.character(index(SP500RM$RV)),as.character(index(snp)))
RV <- as.vector(SP500RM$RV[n]); RQ <- as.vector(SP500RM$RQ[n]);snp <- snp[n]
RV <- RV[1:3874]; RQ <- RQ[1:3874]; lRV <- log(RV)

#MLE for HARK model
Like <- function(v){
  b0 <- v[1]; b1 <- v[2]; b2 <- v[3]; b3 <- v[4]; q <- v[5]; H <- RQ*(2/78)/(RV^2)
  n <- 22; S <- 0
  P <- matrix(0,nrow=22,ncol=22); P <- P+diag(22)*var(lRV)
  A <- c(rep(mean(lRV),22));
  C <- c(b0,rep(0,21)); Z <- t(c(1,rep(0,21)));
  vec <- c(b1,1/4*b2,1/4*b2,1/4*b2,1/4*b2,
           1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,
           1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3)
  T <- diag(21); T <- cbind(T,rep(0,21)); T <- rbind(vec,T);
  Q <- matrix(0,nrow=22,ncol=22); Q[1,1] <- q
  
  for (i in 23:length(lRV)){
    Ap <- C + T%*%A
    Pp <- T%*%P%*%t(T)+Q
    F <- (Z%*%Pp%*%t(Z)+H[i]);
    K <- Pp%*%t(Z)%*%(1/F)
    V <- lRV[i]-Z%*%Ap
    A <- Ap+K%*%V
    P <- Pp-K%*%Z%*%Pp
    S <- n/2*log(2*pi)+1/2*(log(F)+V*(1/F)*V)+S
  }
  return(c(S))
}

Like(c(0.2,0.2,0.2,0.2,0.2))

start_time <- Sys.time()
result <- optim(c(2,2,2,2,2),Like, method="BFGS", hessian = TRUE)
results.HARK <- result$par[1:4]
Sys.time()-start_time

# library(optextras)
# optsp$deps <- 1e-2
# grad_forward <- function(v){grfwd(v,Like,env = optsp)}
# result <- suppressWarnings(optim(c(0,0,0,0,0),gr = function(...) grad_forward(...),Like, method = "BFGS"))

#Prediction w/ MGF method
b0 <- result$par[1]; b1 <- result$par[2]; b2 <- result$par[3]; b3 <- result$par[4]; q <- result$par[5]; h <- RQ*(2/78)/(RV^2)
lRV.check <- 0; lRV.hat.HARK <- 0; RV.hat.HARK <- 0
P <- matrix(0,nrow=22,ncol=22); P <- P+diag(22)*var(lRV)
A <- c(rep(mean(lRV),22));
C <- c(b0,rep(0,21)); H <- h; Z <- t(c(1,rep(0,21))) 
vec <- c(b1,1/4*b2,1/4*b2,1/4*b2,1/4*b2,
         1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,
         1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3)
T <- diag(21); T <- cbind(T,rep(0,21)); T <- rbind(vec,T)
Q <- matrix(0,nrow=22,ncol=22); Q[1,1] <- q; At <- list(); Pt <- list()

for (i in 23:length(RV)){
  Ap <- C + T%*%A
  Pp <- T%*%P%*%t(T)+Q
  F <- Z%*%Pp%*%t(Z)+H[i]
  K <- Pp%*%t(Z)%*%solve(F)
  V <- lRV[i]-Z%*%Ap
  A <- Ap+K%*%V; At[[i]] <- A
  P <- Pp-K%*%Z%*%Pp; Pt[[i]] <- P
  RV.hat.HARK[i] <- exp(Z%*%Ap+1/2*Z%*%Pp%*%t(Z))
  lRV.hat.HARK[i] <- Z%*%Ap
}

#Residuals of lRV
res.HARK <- lRV[1:length(RV)]-lRV.hat.HARK[1:length(RV)]
(MSE.HARK <- mean(res.HARK^2))

par(mfrow=c(2,1),mai=c(0.15,0.4,0.1,0.1),mgp=c(0.1,0.5,0.1))
plot(res.HARK,type="l",xlab="Observations",ylab="Residuals",main="HARK")
acf(res.HARK)

#Residuals
res.HARK <- RV[1:length(RV)]-RV.hat.HARK[1:length(RV)]
(MSE.HARK <- mean((res.HARK^2)[24:length(res.HARK)])) #MSE and QLIKE are computed here from for RV[24:length(RV)] for consistency in model comparison
(QLIKE.HARK <- mean(RV[24:length(RV)]/RV.hat.HARK[24:length(RV)]-log(RV[24:length(RV)]/RV.hat.HARK[24:length(RV)])-1))

#Estimation error
variance <- diag(solve(result$hessian))
(SE.HARK <- sqrt(variance[1:4]))
results.HARK
2*pnorm(abs(results.HARK/SE.HARK[1:4]), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
results.HARK-1.96*SE.HARK[1:4]; results.HARK+1.96*SE.HARK[1:4] # 95% confidence interval 

#Plots
plot(RV[1:length(RV)],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.HARK[1:length(RV)],col="red")

################################################################################

#Out-of-sample
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV); RQ <- as.vector(SP500RM$RQ[n])

#Prediction w/ MGF method
b0 <- result$par[1]; b1 <- result$par[2]; b2 <- result$par[3]; b3 <- result$par[4]; q <- result$par[5]; h <- RQ*(2/78)/(RV^2)
lRV.check <- 0; lRV.hat.HARK <- 0; RV.hat.HARK <- 0
P <- matrix(0,nrow=22,ncol=22); P <- P+diag(22)*var(lRV)
A <- c(rep(mean(lRV),22));
C <- c(b0,rep(0,21)); H <- h; Z <- t(c(1,rep(0,21))) 
vec <- c(b1,1/4*b2,1/4*b2,1/4*b2,1/4*b2,
         1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,
         1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3)
T <- diag(21); T <- cbind(T,rep(0,21)); T <- rbind(vec,T)
Q <- matrix(0,nrow=22,ncol=22); Q[1,1] <- q; At <- list(); Pt <- list()

for (i in 23:length(RV)){
  Ap <- C + T%*%A
  Pp <- T%*%P%*%t(T)+Q
  F <- Z%*%Pp%*%t(Z)+H[i]
  K <- Pp%*%t(Z)%*%solve(F)
  V <- lRV[i]-Z%*%Ap
  A <- Ap+K%*%V; At[[i]] <- A
  P <- Pp-K%*%Z%*%Pp; Pt[[i]] <- P
  RV.hat.HARK[i] <- exp(Z%*%Ap+1/2*Z%*%Pp%*%t(Z))
  lRV.hat.HARK[i] <- Z%*%Ap
}

#Residuals
res.SHARK <- RV[3897:4096]-RV.hat.HARK[3897:4096]
(MSE.SHARK <- mean(res.SHARK^2)) 
(QLIKE.SHARK <- mean(RV[3897:4096]/RV.hat.HARK[3897:4096]-log(RV[3897:4096]/RV.hat.HARK[3897:4096])-1))

#Plots
plot(RV[3897:4096],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.HARK[3897:4096],col="red")

################################################################################

#VaR backtesting
library("rugarch")

var95 <- -1.52*sqrt(RV.hat.HARK[3897:4096]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[3897:4096])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[3897:4096],type="l",ylim=c(-5,5))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$cc.LRp
