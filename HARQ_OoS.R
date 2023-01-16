# Description: HARQ model - Out-of-sample forecast
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 11/10/2022
# Edits: 

#Importing data
library(HARModel); library(xts); data(SP500RM);
setwd("/Users/sarkaajsingh/Desktop/")
snp <- read.csv("snp_daily_returns.csv",sep=""); 
snp$returns <- (log(snp$GSPC.Close)-log(snp$GSPC.Open))*100
snp$Index <- strptime(as.character(snp$Index), "%Y-%m-%d") #Parsing the column as a date
snp <- xts(snp$returns,order.by = snp[,1]) #Converting df into ts
n <- intersect(as.character(index(SP500RM$RV)),as.character(index(snp)))
RV <- as.vector(SP500RM$RV[n]); RQ <- sqrt(as.vector(SP500RM$RQ[n])); snp <- snp[n]
#Daily RV averaged from past 22 days
RV.m <- c(); for (i in 1:(length(RV)-21)){RV.m[i+21] <- mean(RV[i:(i+21)])}
#Daily RV averaged from past 5 days
RV.w <- c(); for (i in 1:(length(RV)-4)){RV.w[i+4] <- mean(RV[i:(i+4)])}
#OLS, estimation
rw.results <- list(); rw.X <- list()
rw.start <- 721; rw.end <- 920; window <- 500 #forecasting 12/04/2011 to 26/01/2012 (500days estimation) rolling from 21/04/2009 to 03/02/2010
points <- rw.start:rw.end #!!!Run this for continuous Oos
points <- 24:(length(RV)-window) #!!!Run this for entire sample
points <- 24:(length(RV)-window); points <- points[(length(points)-199):length(points)] #!!!Run this for last 200 samples
for (rw in points){
  y <- c(RV[(rw+1):(rw+window-1)])
  x1 <- c(RV[(rw):(rw+window-2)]*RQ[(rw):(rw+window-2)])
  x2 <- c(RV[(rw):(rw+window-2)])
  x3 <- c(RV.w[(rw):(rw+window-2)])
  x4 <- c(RV.m[(rw):(rw+window-2)])
  X <- cbind(1,x1,x2,x3,x4)
  Beta.HARQ <- solve((t(X)%*%X))%*%t(X)%*%y; rw.results[[rw]] <- as.vector(Beta.HARQ); rw.X[[rw]] <- X
}
# #Prediction
# rw.RV.hat.HARQ <- c(0)
# for (rw in rw.start:rw.end){
#   rw.RV.hat.HARQ[rw+window] <- rw.X[[rw]]%*%rw.results[[rw]]
# }
#ALT Prediction
rw.RV.hat.HARQ <-c(0)
for (rw in points){
  rw.RV.hat.HARQ[rw+window] <- (rw.results[[rw]][1]+rw.results[[rw]][3]*RV[rw+window-1]+rw.results[[rw]][2]*RV[rw+window-1]*RQ[rw+window-1]+rw.results[[rw]][4]*mean(RV[(window-1-4+rw):(window-1+rw)])+rw.results[[rw]][5]*mean(RV[(rw+window-1-21):(window-1+rw)]))
}
rw.RV.hat.HARQ[(rw.RV.hat.HARQ<=0)] <- mean(rw.RV.hat.HARQ,na.rm=TRUE) #Insanity filter
#Residuals
res.HARQ <- RV[points+window]-rw.RV.hat.HARQ[points+window]
(MSE.HARQ <- mean(res.HARQ^2))
(QLIKE.HARQ <- mean(RV[points+window]/rw.RV.hat.HARQ[points+window]-log(RV[points+window]/rw.RV.hat.HARQ[points+window])-1))


#Plots
plot(RV[points+window],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(rw.RV.hat.HARQ[points+window],col="red")

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[points+window],RV[points+window]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[points+window],rw.RV.hat.HARQ[points+window]); colnames(y2) <- c("date","RV")

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

var95 <- -1.52*sqrt(rw.RV.hat.HARQ[points+window]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[points+window])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[points+window],type="l",ylim=c(-15,10))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$cc.LRp
