# Description: HAR model - log(RV)
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 01 Trimming dataset, QLIKE correction.

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
#Daily RV averaged from past 22 days
lRV.m <- c(); for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}
#Daily RV averaged from past 5 days
lRV.w <- c(); for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}
#OLS, estimation over entire sample size
y <- c(lRV[23:length(RV)])
x1 <- c(lRV[22:(length(RV)-1)])
x2 <- c(lRV.w[22:(length(RV)-1)])
x3 <- c(lRV.m[22:(length(RV)-1)])
X <- cbind(1,x1,x2,x3)
Beta.HAR <- solve((t(X)%*%X))%*%t(X)%*%y; results.HAR <- as.vector(Beta.HAR)
#Prediction
lRV.hat.HAR <- X%*%Beta.HAR
res.lHAR <- lRV[23:length(lRV)]-lRV.hat.HAR
q <- 1/(length(res.lHAR)-4)*(t(res.lHAR)%*%res.lHAR);
#ALT Prediction
RV.hat.HAR <-c(0)
for (i in 1:length(22:(length(RV)-1))){
  RV.hat.HAR[i+22] <- exp(results.HAR[1]+results.HAR[2]*lRV[i+21]+results.HAR[3]*mean(lRV[(17+i):(21+i)])
                       +results.HAR[4]*mean(lRV[(i):(21+i)])+q/2)
}
RV.hat.HAR <- RV.hat.HAR[23:length(RV)]
#Residuals
res.HAR <- RV[23:length(RV)]-RV.hat.HAR
(MSE.HAR <- mean(res.HAR^2)) #MSE and QLIKE are computed here from for RV[24:length(RV)] for consistency in model comparison
(QLIKE.HAR <- mean(RV[23:length(RV)]/RV.hat.HAR[1:length(RV.hat.HAR)]-log(RV[23:length(RV)]/RV.hat.HAR[1:length(RV.hat.HAR)])-1))
#Estimation error standard method
sample.variance.of.error.term <- 1/(length(res.lHAR)-4)*(t(res.lHAR)%*%res.lHAR)
covariance.matrix.of.Beta <- sample.variance.of.error.term[1]*solve(t(X)%*%X)
(SE.HAR <- sqrt(diag(covariance.matrix.of.Beta)))
Beta.HAR

2*pnorm(abs(Beta.HAR/SE.HAR), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
Beta.HAR-1.96*SE.HAR; Beta.HAR+1.96*SE.HAR # 95% confidence interval 

#Plots
plot(RV[23:length(RV)],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.HAR[1:length(RV)],col="red")

################################################################################
#Out-of-sample
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV)

#Prediction
RV.hat.HAR <-c(0)
for (i in 3875:4074){
  RV.hat.HAR[i+22] <- exp(results.HAR[1]+results.HAR[2]*lRV[i+21]+results.HAR[3]*mean(lRV[(17+i):(21+i)])
                          +results.HAR[4]*mean(lRV[(i):(21+i)])+q/2)
}

#Residuals
res.HAR <- RV[3897:4096]-RV.hat.HAR[3897:4096]
(MSE.HAR <- mean(res.HAR^2))
(QLIKE.HAR <- mean(RV[3897:4096]/RV.hat.HAR[3897:4096]-log(RV[3897:4096]/RV.hat.HAR[3897:4096])-1))

#Plots
plot(RV[3897:4096],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.HAR[3897:4096],col="red")

################################################################################

#VaR backtesting
library("rugarch")

var95 <- -1.52*sqrt(RV.hat.HAR[3897:4096]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[3897:4096])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[3897:4096],type="l",ylim=c(-5,5))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$cc.LRp

