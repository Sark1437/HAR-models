# Description: HAR model - Out-of-sample forecast - log(RV)
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
#Daily RV averaged from past 22 days
lRV.m <- c(); for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}
#Daily RV averaged from past 5 days
lRV.w <- c(); for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}
#OLS, estimation
rw.results <- list(); rw.X <- list(); rw.q <- list()
rw.start <- 721; rw.end <- 920; window <- 500 #forecasting 12/04/2011 to 26/01/2012 (500days estimation) rolling from 21/04/2009 to 03/02/2010
points <- rw.start:rw.end #!!!Run this for continuous Oos
points <- 24:(length(RV)-window) #!!!Run this for entire sample
points <- 24:(length(RV)-window); points <- points[(length(points)-199):length(points)] #!!!Run this for last 200 samples
for (rw in points){
  y <- c(lRV[(rw+1):(rw+window-1)])
  x1 <- c(lRV[(rw):(rw+window-2)])
  x2 <- c(lRV.w[(rw):(rw+window-2)])
  x3 <- c(lRV.m[(rw):(rw+window-2)])
  X <- cbind(1,x1,x2,x3)
  Beta.HAR <- solve((t(X)%*%X))%*%t(X)%*%y; rw.results[[rw]] <- as.vector(Beta.HAR); rw.X[[rw]] <- X
  lRV.hat.HAR <- X%*%Beta.HAR
  res.HAR <- lRV[(rw+1):(rw+window-1)]-lRV.hat.HAR
  q <- 1/(length(res.HAR)-4)*(t(res.HAR)%*%res.HAR); rw.q[[rw]] <- q
}
#Prediction
rw.RV.hat.HAR <-c(0)
for (rw in points){
  rw.RV.hat.HAR[rw+window] <- exp(rw.results[[rw]][1]+rw.results[[rw]][2]*lRV[rw+window-1]+rw.results[[rw]][3]*mean(lRV[(window-1-4+rw):(window-1+rw)])
                               +rw.results[[rw]][4]*mean(lRV[(rw+window-1-21):(window-1+rw)])+rw.q[[rw]]/2)
}
#Residuals
res.HAR <- RV[points+window]-rw.RV.hat.HAR[points+window]
(MSE.HAR <- mean(res.HAR^2))
(QLIKE.HAR <- mean(RV[points+window]/rw.RV.hat.HAR[points+window]-log(RV[points+window]/rw.RV.hat.HAR[points+window])-1))

#Plots
plot(RV[points+window],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(rw.RV.hat.HAR[points+window],col="red")

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[points+window],RV[points+window]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[points+window],rw.RV.hat.HAR[points+window]); colnames(y2) <- c("date","RV")

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

var95 <- -1.52*sqrt(rw.RV.hat.HAR[points+window]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[points+window])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[points+window],type="l",ylim=c(-15,10))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$cc.LRp
