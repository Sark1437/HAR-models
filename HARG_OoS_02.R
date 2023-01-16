# Description: logHAR-GARCH(1,1) model using conditional MLE- Out-of-sample forecast
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 01 Trimming dataset, QLIKE correction.
#        02 VaR addition

#Importing data
library(HARModel); library(xts); data(SP500RM); 
setwd("/Users/sarkaajsingh/Desktop/")
snp <- read.csv("snp_daily_returns.csv",sep=""); 
snp$returns <- (log(snp$GSPC.Close)-log(snp$GSPC.Open))*100
snp$Index <- strptime(as.character(snp$Index), "%Y-%m-%d") #Parsing the column as a date
snp <- xts(snp$returns,order.by = snp[,1]) #Converting df into ts
n <- intersect(as.character(index(SP500RM$RV)),as.character(index(snp)))
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV); snp <- snp[n]

#Daily lRV averaged from past 22 days
lRV.m <- c(); for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}
#Daily lRV averaged from past 5 days
lRV.w <- c(); for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}
#Conditional Likelihood function for HAR-GARCH(1,1) 
rw.results.HARG <- list()
rw.start <- 721; rw.end <- 920; window <- 500 #forecasting 12/04/2011 to 26/01/2012 (500days estimation) rolling from 21/04/2009 to 03/02/2010
points <- rw.start:rw.end #!!!Run this for continuous Oos
points <- 24:(length(RV)-window) #!!!Run this for entire sample
points <- 24:(length(RV)-window); points <- points[(length(points)-199):length(points)] #!!!Run this for last 200 samples
for (rw in points){
  cat("rw:", which(points==rw), "\n")
  Like <- function(v){
    b0 <- v[1]; b1 <- v[2]; b2 <- v[3]; b3 <- v[4]
    omega <- exp(v[5]); alpha <- exp(v[6]); beta <-exp(v[7])
    q <- var(lRV); S <- 0
    
    for (i in rw:(rw+window-1)){
      q <- omega+alpha*(lRV[i-1]-(b0+b1*lRV[i-2]+b2*lRV.w[i-2]+b3*lRV.m[i-2]))^2+beta*q
      V <- lRV[i]-(b0+b1*lRV[i-1]+b2*lRV.w[i-1]+b3*lRV.m[i-1])
      S <- 1/2*(log(2*pi)+log(q)+V^2/q)+S
    }
    return(S)
  }
  
  result <- optim(c(0.1,0.1,0.1,0.1,log(0.1),log(0.1),log(0.1)),Like, method="BFGS")
  rw.results.HARG[[rw]] <- c(result$par[1:4],exp(result$par[5:7]))
}
#Computing time-varying variance, qt
rw.qt <- c(0)
for (rw in points){
  b0 <- rw.results.HARG[[rw]][1]; b1 <- rw.results.HARG[[rw]][2]; b2 <- rw.results.HARG[[rw]][3]; b3 <- rw.results.HARG[[rw]][4]
  omega <- rw.results.HARG[[rw]][5]; alpha <- rw.results.HARG[[rw]][6]; beta <-rw.results.HARG[[rw]][7]
  q <- var(lRV); S <- 0; qt <- c(0)
  
  for (i in rw:(rw+window)){
    q <- omega+alpha*(lRV[i-1]-(b0+b1*lRV[i-2]+b2*lRV.w[i-2]+b3*lRV.m[i-2]))^2+beta*q
    V <- lRV[i]-(b0+b1*lRV[i-1]+b2*lRV.w[i-1]+b3*lRV.m[i-1])
    S <- 1/2*(log(2*pi)+log(q)+V^2/q)+S
  }
 rw.qt[rw+window] <- q
} 

#Prediction w/ MGF method
rw.RV.hat.HARG <-c(0)
for (rw in points){
  rw.RV.hat.HARG[rw+window] <- exp(rw.results.HARG[[rw]][1]+rw.results.HARG[[rw]][2]*lRV[rw+window-1]+rw.results.HARG[[rw]][3]*mean(lRV[(window-1-4+rw):(window-1+rw)])
                              +rw.results.HARG[[rw]][4]*mean(lRV[(rw+window-1-21):(window-1+rw)])+rw.qt[rw+window]/2)
}

#Residuals
res.HARG <- RV[points+window]-rw.RV.hat.HARG[points+window]
(MSE.HARG <- mean(res.HARG^2))
(QLIKE.HARG <- mean(RV[points+window]/rw.RV.hat.HARG[points+window]-log(RV[points+window]/rw.RV.hat.HARG[points+window])-1))


#Plots
plot(RV[points+window],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(rw.RV.hat.HARG[points+window],col="red")

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[points+window],RV[points+window]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[points+window],rw.RV.hat.HARG[points+window]); colnames(y2) <- c("date","RV")

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

var95 <- -1.52*sqrt(rw.RV.hat.HARG[points+window]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[points+window])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[points+window],type="l",ylim=c(-15,10))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$cc.LRp
