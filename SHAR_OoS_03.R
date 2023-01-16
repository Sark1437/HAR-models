# Description: SHAR model using conditional MLE (uses log(RV)) Out-of-sample forecast
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 01 Trimming dataset, QLIKE correction.
#        02 VaR addition

#Importing data
library(HARModel); library(xts); data(SP500RM); 
setwd("C:/Users/sarka/Desktop")
snp <- read.csv("snp_daily_returns.csv",sep=""); 
snp$returns <- (log(snp$GSPC.Close)-log(snp$GSPC.Open))*100
snp$Index <- strptime(as.character(snp$Index), "%Y-%m-%d") #Parsing the column as a date
snp <- xts(snp$returns,order.by = snp[,1]) #Converting df into ts
n <- intersect(as.character(index(SP500RM$RV)),as.character(index(snp)))
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV); snp <- snp[n]

#Daily lRV averaged from past 22 days
lRV.m <- c()
for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}

#Daily lRV averaged from past 5 days
lRV.w <- c()
for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}

#load "pracma for calculation of pseudoinverse of Fisher information matrix
library("pracma")

rw.results.SHAR <- list(); rw.Beta.SHAR <- list()
rw.start <- 1243; rw.end <- 1442; window <- 500 #forecasting 12/04/2011 to 26/01/2012 (500days estimation) rolling from 21/04/2009 to 03/02/2010
points <- rw.start:rw.end #!!!Run this for continuous Oos
points <- 24:(length(RV)-window) #!!!Run this for entire sample
points <- 24:(length(RV)-window); points <- points[(length(points)-199):length(points)] #!!!Run this for last 200 samples
for (rw in points){ 
  cat("rw:", which(points==rw), "\n")
  #Pre-sample estimation for initial b0,b1,b2,b3 values
  #OLS, estimation over entire sample size
  y <- c(lRV[(rw-1):(rw+499-1)])
  x1 <- c(lRV[(rw-2):(rw+499-2)])
  x2 <- c(lRV.w[(rw-2):(rw+499-2)])
  x3 <- c(lRV.m[(rw-2):(rw+499-2)])
  X <- cbind(1,x1,x2,x3)
  
  Beta <- solve((t(X)%*%X))%*%t(X)%*%y
  
  #Conditional Likelihood function for SHAR
  Like <- function(v){
    b0 <- Beta[1]; b1 <- Beta[2]; b2 <- Beta[3]; b3 <- Beta[4] 
    omega<-(0);a1<-(v[1]);a2<-(v[2]);a3<-(v[3]);a4<-(v[4]);a5<-(v[5])
    q <- var(lRV); S <- 0;
    A <- c(a1,a2,a3,a4,a5); B <- 1; f <- c(b0,b1,b2,b3,log(q))
    #cat("Parameters:", v)
    for (i in rw:(rw+window-1)){
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
    #cat(", Objective:", S, "\n"); 
    return(S)
  }
  
  #Optimization
  results <- optim(c(0.0,0.0,0.0,0.0,0.0),Like,method = "BFGS") 
  rw.results.SHAR[[rw]] <- results$par
  rw.Beta.SHAR[[rw]] <- Beta
}
################################################################################
#Computing TVP using GAS function
rw.RV.hat.SHAR <-c(0)
for (rw in points){
  omega <- 0
  A <- c(rw.results.SHAR[[rw]][1:5])
  B <- 1
  F.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
  s.SHAR <- matrix(0,nrow=5,ncol=length(lRV))
  f <-  c(rw.Beta.SHAR[[rw]][1],rw.Beta.SHAR[[rw]][2],rw.Beta.SHAR[[rw]][3],rw.Beta.SHAR[[rw]][4],log(var(lRV)))
  F.SHAR[,(rw-1)] <- f
  
  for (i in rw:(rw+window-1)){
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
  }
  #Prediction w/ MGF method
  rw.RV.hat.SHAR[rw+window] <- exp(F.SHAR[1,(rw+window-1)]+F.SHAR[2,(rw+window-1)]*lRV[(rw+window-1)]+F.SHAR[3,(rw+window-1)]*mean(lRV[(rw+window-1-4):(rw+window-1)])
                                   +F.SHAR[4,(rw+window-1)]*mean(lRV[(rw+window-1-21):(rw+window-1)])+exp(F.SHAR[5,(rw+window-1)])/2)
}

#Residuals
res.SHAR <- RV[points+window]-rw.RV.hat.SHAR[points+window]
(MSE.SHAR <- mean(res.SHAR^2))
(QLIKE.SHAR <- mean(RV[points+window]/rw.RV.hat.SHAR[points+window]-log(RV[points+window]/rw.RV.hat.SHAR[points+window])-1))

#Plots
plot(RV[points+window],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(rw.RV.hat.SHAR[points+window],col="red")

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[points+window],RV[points+window]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[points+window],rw.RV.hat.SHAR[points+window]); colnames(y2) <- c("date","RV")

ggplot()+
  geom_line(data = y1, aes(date, RV), color = "black")+
  geom_line(data = y2, aes(date, RV), color = "red")+
  labs(y="RV", x="Year")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(size = 15))

################################################################################
#library(quantmod) #install.packages("quantmod") again to update if could not establish session error.
#snp <- getSymbols("^GSPC", auto.assign=FALSE, from = "1997-04-08", to = "2013-08-31", src="yahoo")
#write.zoo(snp,file="snp_daily_returns")
library("rugarch")

var95 <- -1.52*sqrt(rw.RV.hat.SHAR[points+window]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[points+window])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[points+window],type="l",ylim=c(-15,10))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$cc.LRp
