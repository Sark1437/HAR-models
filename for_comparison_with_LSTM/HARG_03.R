# Description: logHAR-GARCH(1,1) model using conditional MLE
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 03 Trimming dataset, QLIKE correction.

#Importing data
library(HARModel); library(xts); data(SP500RM);
setwd("/Users/sarkaajsingh/Desktop/")
snp <- read.csv("snp_daily_returns.csv",sep=" "); 
snp$returns <- (log(snp$GSPC.Close)-log(snp$GSPC.Open))*100
snp$Index <- strptime(as.character(snp$Index), "%Y-%m-%d") #Parsing the column as a date
snp <- xts(snp$returns,order.by = snp[,1]) #Converting df into ts
n <- intersect(as.character(index(SP500RM$RV)),as.character(index(snp)))
RV <- as.vector(SP500RM$RV[n]); snp <- snp[n]
RV <- RV[1:3874]; lRV <- log(RV)
#Daily lRV averaged from past 22 days
lRV.m <- c(); for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}
#Daily lRV averaged from past 5 days
lRV.w <- c(); for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}
#Conditional Likelihood function for HAR-GARCH(1,1)
Like <- function(v){
  b0 <- v[1]; b1 <- v[2]; b2 <- v[3]; b3 <- v[4]
  omega <- exp(v[5]); alpha <- exp(v[6]); beta <- exp(v[7])
  q <- var(lRV); S <- 0
  
  for (i in 24:length(lRV)){
    q <- omega+alpha*(lRV[i-1]-(b0+b1*lRV[i-2]+b2*lRV.w[i-2]+b3*lRV.m[i-2]))^2+beta*q
    V <- lRV[i]-(b0+b1*lRV[i-1]+b2*lRV.w[i-1]+b3*lRV.m[i-1])
    S <- 1/2*(log(2*pi)+log(q)+V^2/q)+S
  }
  return(S)
}

Like(c(0.2,0.2,0.2,0.2,log(0.2),log(0.2),log(0.2)))
result <- optim(c(0.1,0.1,0.1,0.1,log(0.1),log(0.1),log(0.1)),Like, method="BFGS",hessian = TRUE)
results <- c(result$par[1:4],exp(result$par[5:7]))
results.HARG <- results

#Computing time-varying variance, qt
b0 <- results[1]; b1 <- results[2]; b2 <- results[3]; b3 <- results[4]
omega <- results[5]; alpha <- (results[6]); beta <-results[7]
q <- var(lRV); S <- 0; qt <- c(0)

for (i in 24:length(lRV)){
  q <- omega+alpha*(lRV[i-1]-(b0+b1*lRV[i-2]+b2*lRV.w[i-2]+b3*lRV.m[i-2]))^2+beta*q
  V <- lRV[i]-(b0+b1*lRV[i-1]+b2*lRV.w[i-1]+b3*lRV.m[i-1])
  S <- 1/2*(log(2*pi)+log(q)+V^2/q)+S
  qt[i] <- q
}

#Prediction of lRV
lRV.hat.HARG <-c(0)
for (i in 2:length(23:length(lRV))){
  lRV.hat.HARG[i+22] <- results[1]+results[2]*lRV[i+21]+results[3]*mean(lRV[(17+i):(21+i)])
                           +results[4]*mean(lRV[(i):(21+i)])
}

#Residuals of lRV
res.HARG <- lRV[24:length(lRV)]-lRV.hat.HARG[24:length(lRV)]
(MSE.HARG <- mean(res.HARG^2))

par(mfrow=c(2,1),mai=c(0.15,0.4,0.1,0.1),mgp=c(0.1,0.5,0.1))
plot(res.HARG,type="l",xlab="Observations",ylab="Residuals",main="HARG")
acf(res.HARG)

#Prediction w/ MGF method
RV.hat.HARG <-c(0)
for (i in 2:length(23:length(lRV))){
  RV.hat.HARG[i+22] <- exp(results[1]+results[2]*lRV[i+21]+results[3]*mean(lRV[(17+i):(21+i)])
                         +results[4]*mean(lRV[(i):(21+i)])+qt[i+22]/2)
}

#Residuals of RV
res.HARG <- RV[24:length(lRV)]-RV.hat.HARG[24:length(lRV)]
(MSE.HARG <- mean(res.HARG^2))
(QLIKE.HARG <- mean(RV[24:length(lRV)]/RV.hat.HARG[24:length(lRV)]-log(RV[24:length(lRV)]/RV.hat.HARG[24:length(lRV)])-1))

#Estimation error (correct way to account for exp in the likelihood)
Like <- function(v){
  b0 <- v[1]; b1 <- v[2]; b2 <- v[3]; b3 <- v[4]
  omega <- (v[5]); alpha <- (v[6]); beta <- (v[7])
  q <- var(lRV); S <- 0
  
  for (i in 24:length(lRV)){
    q <- omega+alpha*(lRV[i-1]-(b0+b1*lRV[i-2]+b2*lRV.w[i-2]+b3*lRV.m[i-2]))^2+beta*q
    V <- lRV[i]-(b0+b1*lRV[i-1]+b2*lRV.w[i-1]+b3*lRV.m[i-1])
    S <- 1/2*(log(2*pi)+log(q)+V^2/q)+S
  }
  return(S)
}
library(numDeriv)
hess <- hessian(Like,results) #Run Like without exp first before this step
variance <- diag(solve(hess))
(SE.HARG <- sqrt(variance))
results.HARG
2*pnorm(abs(results.HARG/SE.HARG), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
results.HARG-1.96*SE.HARG; results.HARG+1.96*SE.HARG # 95% confidence interval 

#Plots
plot(RV[24:length(lRV)],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.HARG[24:length(lRV)],col="red")

################################################################################
#Out-of-sample
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV)

#Daily lRV averaged from past 22 days
lRV.m <- c(); for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}

#Daily lRV averaged from past 5 days
lRV.w <- c(); for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}

#Computing time-varying variance, qt
b0 <- results[1]; b1 <- results[2]; b2 <- results[3]; b3 <- results[4]
omega <- results[5]; alpha <- (results[6]); beta <-results[7]
q <- var(lRV); S <- 0; qt <- c(0)

for (i in 24:length(lRV)){
  q <- omega+alpha*(lRV[i-1]-(b0+b1*lRV[i-2]+b2*lRV.w[i-2]+b3*lRV.m[i-2]))^2+beta*q
  V <- lRV[i]-(b0+b1*lRV[i-1]+b2*lRV.w[i-1]+b3*lRV.m[i-1])
  S <- 1/2*(log(2*pi)+log(q)+V^2/q)+S
  qt[i] <- q
}

#Prediction
RV.hat.HARG <-c(0)
for (i in 3875:4074){
  RV.hat.HARG[i+22] <- exp(results[1]+results[2]*lRV[i+21]+results[3]*mean(lRV[(17+i):(21+i)])
                           +results[4]*mean(lRV[(i):(21+i)])+qt[i+22]/2)
}

#Residuals of RV
res.HARG <- RV[3897:4096]-RV.hat.HARG[3897:4096]
(MSE.HARG <- mean(res.HARG^2))
(QLIKE.HARG <- mean(RV[3897:4096]/RV.hat.HARG[3897:4096]-log(RV[3897:4096]/RV.hat.HARG[3897:4096])-1))

#Plots
plot(RV[3897:4096],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.HARG[3897:4096],col="red")

################################################################################
#VaR backtesting
library("rugarch")

var95 <- -1.52*sqrt(RV.hat.HARG[3897:4096]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[3897:4096])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[3897:4096],type="l",ylim=c(-5,5))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$cc.LRp
