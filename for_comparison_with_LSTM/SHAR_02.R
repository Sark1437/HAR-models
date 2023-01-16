# Description: SHAR model using conditional MLE (uses log(RV))
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 02 Trimming dataset, QLIKE correction.

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
lRV.m <- c()
for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}

#Daily lRV averaged from past 5 days
lRV.w <- c()
for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}

#Pre-sample estimation for initial b0,b1,b2,b3 values
#OLS, estimation over entire sample size
y <- c(lRV[23:522])
x1 <- c(lRV[22:521])
x2 <- c(lRV.w[22:521])
x3 <- c(lRV.m[22:521])
X <- cbind(1,x1,x2,x3)

Beta <- solve((t(X)%*%X))%*%t(X)%*%y

#load "pracma for calculation of pseudoinverse of Fisher information matrix
library("pracma")

#Conditional Likelihood function for SHAR
Like <- function(v){
  b0 <- Beta[1]; b1 <- Beta[2]; b2 <- Beta[3]; b3 <- Beta[4] 
  omega<-(0);a1<-(v[1]);a2<-(v[2]);a3<-(v[3]);a4<-(v[4]);a5<-(v[5])
  q <- var(lRV); S <- 0;
  A <- c(a1,a2,a3,a4,a5); B <- 1; f <- c(b0,b1,b2,b3,log(q))
  cat("Parameters:", v)
  for (i in 23:length(lRV)){
    if (exp(f[5])==Inf | is.nan(exp(f[5])) | exp(f[5])==0 |  1/exp(f[5])==Inf){S=Inf; break}
    V <- lRV[i]-(f[1]+f[2]*lRV[i-1]+f[3]*lRV.w[i-1]+f[4]*lRV.m[i-1])
    S <- 1/2*(log(2*pi)+log(exp(f[5]))+V^2/exp(f[5]))+S
    score <- c(V/exp(f[5]),V*lRV[i-1]/exp(f[5]),V*lRV.w[i-1]/exp(f[5]),V*lRV.m[i-1]/exp(f[5]),(-1/2+V^2/(2*exp(f[5]))))
    I <- 1/exp(f[5])*(cbind(c(1,lRV[i-1],lRV.w[i-1],lRV.m[i-1],0),
                            c(lRV[i-1],lRV[i-1]^2,lRV[i-1]*lRV.w[i-1],lRV[i-1]*lRV.m[i-1],0),
                            c(lRV.w[i-1],lRV[i-1]*lRV.w[i-1],lRV.w[i-1]^2,lRV.w[i-1]*lRV.m[i-1],0),
                            c(lRV.m[i-1],lRV[i-1]*lRV.m[i-1],lRV.w[i-1]*lRV.m[i-1],lRV.m[i-1]^2,0),
                            c(0,0,0,0,exp(f[5])/2)))
    s <- pinv(I)%*%score
    f <- omega+A*s+B*f
  }
  cat(", Objective:", S, "\n"); return(S)
}
Like(c(0.02651402, 0.01075011, -0.01490841, 0.005618987, 0.03354273))

#Optimization
start <- Sys.time()
result_SHAR_BFGS <- optim(c(0.008011510,  0.008780363, -0.022348912,  0.008218722,  0.030742818),Like,method = "BFGS",hessian=TRUE) 
(result_SHAR_BFGS$time <- Sys.time()-start)

# library(optextras)
# optsp$deps <- 1e-2
# grad_forward <- function(v){grfwd(v,Like,env = optsp)}
# start <- Sys.time()
# result_SHAR_BFGS_for <- suppressWarnings(optim(c(0,0,0,0,0),gr = function(...) grad_forward(...),Like, method = "BFGS"))
# (result_SHAR_BFGS_for$time <- Sys.time()-start)
# 
# start <- Sys.time()
# result_SHAR_NM <- optim(c(0,0,0,0,0),Like,method = "Nelder-Mead")
# (result_SHAR_NM$time <- Sys.time()-start)

results <- result_SHAR_BFGS$par
###############################################################################
#Computing TVP using GAS function
omega <- 0
A <- c(results[1:5])
B <- 1
F.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
s.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
score.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
F.SHAR[,22] <- c(Beta[1],Beta[2],Beta[3],Beta[4],log(var(lRV)))
f <-  c(Beta[1],Beta[2],Beta[3],Beta[4],log(var(lRV)))

for (i in 23:length(lRV)){
  V <- lRV[i]-(f[1]+f[2]*lRV[i-1]+f[3]*lRV.w[i-1]+f[4]*lRV.m[i-1])
  score <- c(V/exp(f[5]),V*lRV[i-1]/exp(f[5]),V*lRV.w[i-1]/exp(f[5]),V*lRV.m[i-1]/exp(f[5]),(-1/2+V^2/(2*exp(f[5]))))
  I <- 1/exp(f[5])*(cbind(c(1,lRV[i-1],lRV.w[i-1],lRV.m[i-1],0),
                          c(lRV[i-1],lRV[i-1]^2,lRV[i-1]*lRV.w[i-1],lRV[i-1]*lRV.m[i-1],0),
                          c(lRV.w[i-1],lRV[i-1]*lRV.w[i-1],lRV.w[i-1]^2,lRV.w[i-1]*lRV.m[i-1],0),
                          c(lRV.m[i-1],lRV[i-1]*lRV.m[i-1],lRV.w[i-1]*lRV.m[i-1],lRV.m[i-1]^2,0),
                          c(0,0,0,0,exp(f[5])/2)))
  s <- pinv(I)%*%score
  f <- omega+A*s+B*f
  F.SHAR[,i] <- f
  s.SHAR[,i] <- s
  score.SHAR[,i] <- score
}

#Prediction of lRV
lRV.hat.SHAR <-c(0)
for (i in 1:length(23:(length(RV)-1))){
  lRV.hat.SHAR[i+23] <- F.SHAR[1,(i+22)]+F.SHAR[2,(i+22)]*lRV[i+22]+F.SHAR[3,(i+22)]*mean(lRV[(18+i):(22+i)])
  +F.SHAR[4,(i+22)]*mean(lRV[(01+i):(22+i)])
}

#Residuals of lRV
res.SHAR <- lRV[24:length(RV)]-lRV.hat.SHAR[24:length(RV)]
(MSE.SHAR <- mean(res.SHAR^2))

par(mfrow=c(2,1),mai=c(0.15,0.4,0.1,0.1),mgp=c(0.1,0.5,0.1))
plot(res.SHAR,type="l",xlab="Observations",ylab="Residuals",main="SHAR")
acf(res.SHAR)

#Prediction w/ MGF method
RV.hat.SHAR <-c(0)
for (i in 1:length(23:(length(RV)-1))){
  RV.hat.SHAR[i+23] <- exp(F.SHAR[1,(i+22)]+F.SHAR[2,(i+22)]*lRV[i+22]+F.SHAR[3,(i+22)]*mean(lRV[(18+i):(22+i)])
                           +F.SHAR[4,(i+22)]*mean(lRV[(01+i):(22+i)])+exp(F.SHAR[5,(i+22)])/2)
}

#Residuals
res.SHAR <- RV[24:length(RV)]-RV.hat.SHAR[24:length(RV)]
(MSE.SHAR <- mean(res.SHAR^2))
(QLIKE.SHAR <- mean(RV[24:length(RV)]/RV.hat.SHAR[24:length(RV)]-log(RV[24:length(RV)]/RV.hat.SHAR[24:length(RV)])-1))

#Estimation error
variance <- diag(solve(result_SHAR_BFGS$hessian))
(SE.SHAR <- sqrt(variance))
results
2*pnorm(abs(results/SE.SHAR), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
results-1.96*SE.SHAR; results+1.96*SE.SHAR # 95% confidence interval 

#Plots
plot(RV[24:length(RV)],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.SHAR[24:length(RV)],col="red")

################################################################################

#Out-of-sample
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV)

#Daily lRV averaged from past 22 days
lRV.m <- c()
for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}

#Daily lRV averaged from past 5 days
lRV.w <- c()
for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}

#Computing TVP using GAS function
omega <- 0
A <- c(results[1:5])
B <- 1
F.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
s.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
score.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
F.SHAR[,22] <- c(Beta[1],Beta[2],Beta[3],Beta[4],log(var(lRV)))
f <-  c(Beta[1],Beta[2],Beta[3],Beta[4],log(var(lRV)))

for (i in 23:length(lRV)){
  V <- lRV[i]-(f[1]+f[2]*lRV[i-1]+f[3]*lRV.w[i-1]+f[4]*lRV.m[i-1])
  score <- c(V/exp(f[5]),V*lRV[i-1]/exp(f[5]),V*lRV.w[i-1]/exp(f[5]),V*lRV.m[i-1]/exp(f[5]),(-1/2+V^2/(2*exp(f[5]))))
  I <- 1/exp(f[5])*(cbind(c(1,lRV[i-1],lRV.w[i-1],lRV.m[i-1],0),
                          c(lRV[i-1],lRV[i-1]^2,lRV[i-1]*lRV.w[i-1],lRV[i-1]*lRV.m[i-1],0),
                          c(lRV.w[i-1],lRV[i-1]*lRV.w[i-1],lRV.w[i-1]^2,lRV.w[i-1]*lRV.m[i-1],0),
                          c(lRV.m[i-1],lRV[i-1]*lRV.m[i-1],lRV.w[i-1]*lRV.m[i-1],lRV.m[i-1]^2,0),
                          c(0,0,0,0,exp(f[5])/2)))
  s <- pinv(I)%*%score
  f <- omega+A*s+B*f
  F.SHAR[,i] <- f
  s.SHAR[,i] <- s
  score.SHAR[,i] <- score
}

#Prediction w/ MGF method
RV.hat.SHAR <-c(0)
for (i in 1:length(23:(length(RV)-1))){
  RV.hat.SHAR[i+23] <- exp(F.SHAR[1,(i+22)]+F.SHAR[2,(i+22)]*lRV[i+22]+F.SHAR[3,(i+22)]*mean(lRV[(18+i):(22+i)])
                           +F.SHAR[4,(i+22)]*mean(lRV[(01+i):(22+i)])+exp(F.SHAR[5,(i+22)])/2)
}

#Residuals
res.SHARK <- RV[3897:4096]-RV.hat.SHAR[3897:4096]
(MSE.SHARK <- mean(res.SHARK^2)) 
(QLIKE.SHARK <- mean(RV[3897:4096]/RV.hat.SHAR[3897:4096]-log(RV[3897:4096]/RV.hat.SHAR[3897:4096])-1))

#Plots
plot(RV[3897:4096],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(RV.hat.SHAR[3897:4096],col="red")

################################################################################

#VaR backtesting
library("rugarch")

var95 <- -1.52*sqrt(RV.hat.SHAR[3897:4096]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[3897:4096])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[3897:4096],type="l",ylim=c(-5,5))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[3897:4096]), VaR=var95, conf.level = 0.95)$cc.LRp


