# Description: HARK model Out-of-sample forecast
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
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV); snp <- snp[n]; RQ <- as.vector(SP500RM$RQ[n])

#MLE for HARK model
rw.results.HARK <- list()
rw.start <- 721; rw.end <- 920; window <- 500 #forecasting 12/04/2011 to 26/01/2012 (500days estimation) rolling from 21/04/2009 to 03/02/2010
points <- rw.start:rw.end #!!!Run this for continuous Oos
points <- 24:(length(RV)-window) #!!!Run this for entire sample
points <- 24:(length(RV)-window); points <- points[(length(points)-199):length(points)] #!!!Run this for last 200 samples
for (rw in points){
  cat("rw:", which(points==rw), "\n")
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
    
    #cat("Parameters:", v)
    for (i in rw:(rw+window-1)){
      Ap <- C + T%*%A
      Pp <- T%*%P%*%t(T)+Q
      F <- (Z%*%Pp%*%t(Z)+H[i]);
      K <- Pp%*%t(Z)%*%(1/F)
      #if (abs(K[1])>1){S=Inf; break}
      V <- lRV[i]-Z%*%Ap
      A <- Ap+K%*%V
      P <- Pp-K%*%Z%*%Pp
      S <- n/2*log(2*pi)+1/2*(log(F)+V*(1/F)*V)+S
    }
    #cat(", Objective:", S, "\n"); 
    return(c(S))
  }
  # library(optextras)
  # optsp$deps <- 1e-2
  # grad_forward <- function(v){grfwd(v,Like,env = optsp)}
  # result <- suppressWarnings(optim(c(0,0,0,0,0),gr = function(...) grad_forward(...),Like, method = "BFGS"))
  # rw.results.HARK[[rw]] <- result$par[1:5]
  
  result <- optim(c(2,2,2,2,2),Like, method = "BFGS")
  rw.results.HARK[[rw]] <- result$par[1:5]
}
#Prediction w/ MGF method
rw.RV.hat.HARK <-c(0)
for (rw in points){
  b0 <- rw.results.HARK[[rw]][1]; b1 <- rw.results.HARK[[rw]][2]; b2 <- rw.results.HARK[[rw]][3]; b3 <- rw.results.HARK[[rw]][4]; q <- rw.results.HARK[[rw]][5]; h <- RQ*(2/78)/(RV^2)
  P <- matrix(0,nrow=22,ncol=22); P <- P+diag(22)*var(lRV)
  A <- c(rep(mean(lRV),22));
  C <- c(b0,rep(0,21)); H <- h; Z <- t(c(1,rep(0,21))) 
  vec <- c(b1,1/4*b2,1/4*b2,1/4*b2,1/4*b2,
           1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,
           1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3,1/17*b3)
  T <- diag(21); T <- cbind(T,rep(0,21)); T <- rbind(vec,T)
  Q <- matrix(0,nrow=22,ncol=22); Q[1,1] <- q;
  
  for (i in rw:(rw+window)){
    Ap <- C + T%*%A
    Pp <- T%*%P%*%t(T)+Q
    F <- Z%*%Pp%*%t(Z)+H[i]
    K <- Pp%*%t(Z)%*%solve(F)
    V <- lRV[i]-Z%*%Ap
    A <- Ap+K%*%V
    P <- Pp-K%*%Z%*%Pp
  }
  rw.RV.hat.HARK[rw+window] <- exp(Z%*%Ap+1/2*Z%*%Pp%*%t(Z))
}  
#Residuals
res.HARK <- RV[points+window]-rw.RV.hat.HARK[points+window]
(MSE.HARK <- mean(res.HARK^2))
(QLIKE.HARK <- mean(RV[points+window]/rw.RV.hat.HARK[points+window]-log(RV[points+window]/rw.RV.hat.HARK[points+window])-1))

#Plots
plot(RV[points+window],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(rw.RV.hat.HARK[points+window],col="red")

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[points+window],RV[points+window]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[points+window],rw.RV.hat.HARK[points+window]); colnames(y2) <- c("date","RV")

ggplot()+
  geom_line(data = y1, aes(date, RV), color = "black")+
  geom_line(data = y2, aes(date, RV), color = "red")+
  labs(y="RV", x="Year")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(size = 15))


###############################################################################

#library(quantmod) #install.packages("quantmod") again to update if could not establish session error.
#snp <- getSymbols("^GSPC", auto.assign=FALSE, from = "1997-04-08", to = "2013-08-31", src="yahoo")
#write.zoo(snp,file="snp_daily_returns")
library("rugarch")

var95 <- -1.52*sqrt(rw.RV.hat.HARK[points+window]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[points+window])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[points+window],type="l",ylim=c(-15,10))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$cc.LRp
