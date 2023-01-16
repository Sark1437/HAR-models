# Description: SHARK model using exact MLE (uses log(RV)) using KF from Delle Monache and Corsi 2021 Out-of-sample forecast
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 28/09/2022
# Edits: 03 Trimming dataset, QLIKE correction.
#        04 VaR addition

#Importing data
library(HARModel); library(xts); data(SP500RM); 
setwd("/Users/sarkaajsingh/Desktop")
snp <- read.csv("snp_daily_returns.csv",sep=""); 
snp$returns <- (log(snp$GSPC.Close)-log(snp$GSPC.Open))*100
snp$Index <- strptime(as.character(snp$Index), "%Y-%m-%d") #Parsing the column as a date
snp <- xts(snp$returns,order.by = snp[,1]) #Converting df into ts
n <- intersect(as.character(index(SP500RM$RV)),as.character(index(snp)))
RV <- as.vector(SP500RM$RV[n]); lRV <- log(RV); snp <- snp[n]; RQ <- as.vector(SP500RM$RQ[n])

#Daily lRV averaged from past 22 days
lRV.m <- c()
for (i in 1:(length(lRV)-21)){lRV.m[i+21] <- mean(lRV[i:(i+21)])}
#Daily lRV averaged from past 5 days
lRV.w <- c()
for (i in 1:(length(lRV)-4)){lRV.w[i+4] <- mean(lRV[i:(i+4)])}
#Load "pracma" for calculation of pseudoinverse of Fisher information matrix
library("pracma")
#Load "matrixcal" for commutation matrix
library("matrixcalc"); C_n_n <- commutation.matrix(22, 22); N <- diag(484)+C_n_n
#Load "fastmatrix" for fast kronecker product calculations
library("fastmatrix")

rw.results.SHARK <- list(); rw.Beta.SHARK <- list() #Skip this if rerunning code from saved data
rw.start <- 1243; rw.end <- 1442; window <- 500 #forecasting 12/04/2011 to 26/01/2012 (500days estimation) rolling from 21/04/2009 to 03/02/2010
points <- rw.start:rw.end #!!!Run this for continuous Oos
points <- 24:(length(RV)-window) #!!!Run this for entire sample
points <- 24:(length(RV)-window); points <- points[(length(points)-199):length(points)] #!!!Run this for last 200 samples
start <- Sys.time()
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
  
  ################################################################################
  Like <- function(v){
    #Initialization of parameters
    b0 <- Beta[1]; b1 <- Beta[2]; b2 <- Beta[3]; b3 <- Beta[4]; q <- var(lRV)
    omega<-0;a1<-(v[1]);a2<-(v[2]);a3<-(v[3]);a4<-(v[4]);a5<-(v[5])
    alpha <- c(a1,a2,a3,a4,a5); B <- 1; f <- c(b0,b1,b2,b3,log(q))
    #Variables for KF
    P <- matrix(0,nrow=22,ncol=22); P <- P+diag(22)*var(lRV)
    A <- c(rep(mean(lRV),22)); Z <- t(c(1,rep(0,21))) ;C <- c(f[1],rep(0,21))
    Q <- matrix(0,nrow=22,ncol=22); Q[1,1] <- exp(f[5]); H <- RQ*(2/78)/(RV^2)
    vec <- c(f[2],1/4*f[3],1/4*f[3],1/4*f[3],1/4*f[3],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4])
    T <- diag(21); T <- cbind(T,rep(0,21)); T <- rbind(vec,T);
    #Variables for score function and Fisher information
    dT <- cbind(c(rep(0,484)),c(1,rep(0,483)),1/4*c(rep(0,22),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),rep(0,374)),1/17*c(rep(0,110),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21)),c(rep(0,484)))
    dQ <- cbind(c(rep(0,484)),c(rep(0,484)),c(rep(0,484)),c(rep(0,484)),c(1,rep(0,483)))
    dC <- cbind(c(1,rep(0,21)),c(rep(0,22)),c(rep(0,22)),c(rep(0,22)),c(rep(0,22)))
    dP <- matrix(10000,nrow = 484,ncol=5); dA <- matrix(0,nrow = 22,ncol = 5); Z_kron_Z <- kronecker.prod(Z,Z);
    #Loop initialization
    #cat("Parameters:", v); 
    S <- 0;
    for (i in rw:(rw+window-1)){
      C[1] <- f[1]
      T[1,1] <- f[2]; T[1,2:5] <- 1/4*f[3]; T[1,6:22] <- 1/17*f[4]
      if (exp(f[5])==Inf){S=Inf; break}
      Q[1,1] <- exp(f[5])
      dQ[1,5] <- exp(f[5]) #This is done due to chain rule/jacobian
      
      V <- lRV[i]-Z%*%A
      F <- Z%*%P%*%t(Z)+H[i]
      K <- T%*%P%*%t(Z)%*%(1/F)
      if (abs(K[1])>1){S=Inf; break}
      
      dV <- -Z%*%dA
      dF <- Z_kron_Z%*%dP
      dK <- kronecker.prod((1/F[1])*Z%*%P,diag(22))%*%dT+kronecker.prod((1/F[1])*Z,T)%*%dP-1/F[1]*K%*%dF
      
      score <- 1/2*((t(dF)*1/F[1]^2)%*%vec(V^2-F)-(2*t(dV)*1/F[1])%*%V)
      I <- 1/2*((t(dF)*1/F[1]^2)%*%dF+2*(t(dV)*1/F[1])%*%dV)
      s <- pinv(I)%*%score
      
      dA <- kronecker.prod(t(A),diag(22))%*%dT+dC+T%*%dA+V[1]*dK+K%*%dV #Check is dC + or *
      dP <-N%*%(kronecker.prod(T%*%P,diag(22))%*%dT-kronecker.prod(K%*%F,diag(22))%*%dK)+kronecker.prod(T,T)%*%dP-kronecker.prod(K,K)%*%dF+dQ
      A <- C+T%*%A+K%*%V
      P <- T%*%P%*%t(T-K%*%Z)+Q
      f <- omega+alpha*s+B*f
      S <- 22/2*log(2*pi)+1/2*(log(F)+V*(1/F)*V)+S
    }
    #cat(", Objective:", S, "\n"); 
    return(S)
  }
  
  #Optimization
  # library(optimParallel)
  # cl <- makeCluster(8); setDefaultCluster(cl = cl)
  # clusterExport(cl,c("Beta","lRV","seq","lRV.w","lRV.m","RQ","RV","N"))
  # start <- Sys.time() #Need to run "optimParallel_source first to use "BFGS" instead of "L-BFGS-B"
  # result_SHARK_parallel <- optimParallel(c(0,0,0,0,0),Like, method = "BFGS", parallel = c(loginfo = T,forward = T, clusterEvalQ(cl,c(library("pracma"),library("matrixcalc"),library("fastmatrix")))))
  # (result_SHARK_parallel$time <- Sys.time()-start)
  
  # library(optextras)
  # optsp$deps <- 1e-2
  # grad_forward <- function(v){grfwd(v,Like,env = optsp)}
  # results <- suppressWarnings(optim(c(0,0,0,0,0),gr = function(...) grad_forward(...),Like, method = "BFGS"))
  # rw.results.SHARK[[rw]] <- results$par
  # rw.Beta.SHARK[[rw]] <- Beta
  # if (which(points==rw)%%10==0){save.image("C:/Users/sarka/Desktop/OneDrive/Temp/checkpoint.RData")} #Checkpoint creation
  # if (which(points==rw)%%10==0){save.image("/Users/sarkaajsingh/Library/CloudStorage/OneDrive-Personal/Temp/checkpoint.RData")} #Checkpoint creation
  # print(Sys.time()-start)
  
  results <- optim(c(0.0,0.0,0.0,0.0,0.0),Like, method = "BFGS")
  rw.results.SHARK[[rw]] <- results$par
  rw.Beta.SHARK[[rw]] <- Beta
  if (which(points==rw)%%10==0){save.image("/Users/sarkaajsingh/Desktopcheckpoint.RData")} #Checkpoint creation
  print(Sys.time()-start)
  
  # results <- optim(c(0,0,0,0,0),Like,method = "Nelder-Mead",control = c(reltol = 1e-6))
  # rw.results.SHARK[[rw]] <- results$par
  # rw.Beta.SHARK[[rw]] <- Beta
  # 
  # start <- Sys.time()
  # result_SHARK_NM_full <- optim(c(0,0,0,0,0),Like,method = "Nelder-Mead")
  # (result_SHARK_NM_full$time <- Sys.time()-start)
}

###############################################################################
#Computing TVP using GAS function
#Initialization of parameters
rw.RV.hat.SHARK <-c(0)
for (rw in points){
  omega <- 0; alpha <- c(rw.results.SHARK[[rw]][1:5]); B <- 1
  F.SHARK <- matrix(0,nrow=5,ncol=length(lRV))
  s.SHARK <- matrix(0,nrow=5,ncol=length(lRV))
  At <- matrix(0,nrow=22,ncol=length(lRV)); Pt <- list()
  f <-  c(rw.Beta.SHARK[[rw]][1],rw.Beta.SHARK[[rw]][2],rw.Beta.SHARK[[rw]][3],rw.Beta.SHARK[[rw]][4],log(var(lRV)))
  F.SHARK[,(rw-1)] <- f
  #Variables for KF
  P <- matrix(0,nrow=22,ncol=22); P <- P+diag(22)*var(lRV)
  A <- c(rep(mean(lRV),22)); Z <- t(c(1,rep(0,21))) ;C <- c(f[1],rep(0,21))
  Q <- matrix(0,nrow=22,ncol=22); Q[1,1] <- exp(f[5]); H <- RQ*(2/78)/(RV^2)
  vec <- c(f[2],1/4*f[3],1/4*f[3],1/4*f[3],1/4*f[3],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4],1/17*f[4])
  T <- diag(21); T <- cbind(T,rep(0,21)); T <- rbind(vec,T); Tt <- list(); Ct <- list(); Qt <- list(); dAt <- list(); dPt <- list()
  #Variables for score function and Fisher information
  dT <- cbind(c(rep(0,484)),c(1,rep(0,483)),1/4*c(rep(0,22),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),rep(0,374)),1/17*c(rep(0,110),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21),1,rep(0,21)),c(rep(0,484)))
  dQ <- cbind(c(rep(0,484)),c(rep(0,484)),c(rep(0,484)),c(rep(0,484)),c(1,rep(0,483)))
  dC <- cbind(c(1,rep(0,21)),c(rep(0,22)),c(rep(0,22)),c(rep(0,22)),c(rep(0,22)))
  dP <- matrix(10000,nrow = 484,ncol=5); dA <- matrix(0,nrow = 22,ncol = 5); Z_kron_Z <- kronecker(Z,Z);
  
  for (i in rw:(rw+window-1)){
    C[1] <- f[1]
    T[1,1] <- f[2]; T[1,2:5] <- 1/4*f[3]; T[1,6:22] <- 1/17*f[4]
    if (exp(f[5])==Inf){S=Inf; break}
    Q[1,1] <- exp(f[5])
    dQ[1,5] <- exp(f[5]) #This is done due to chain rule/jacobian
    
    V <- lRV[i]-Z%*%A
    F <- Z%*%P%*%t(Z)+H[i]
    K <- T%*%P%*%t(Z)%*%(1/F)
    if (abs(K[1])>1){S=Inf; break}
    
    dV <- -Z%*%dA
    dF <- Z_kron_Z%*%dP
    dK <- kronecker.prod((1/F[1])*Z%*%P,diag(22))%*%dT+kronecker.prod((1/F[1])*Z,T)%*%dP-1/F[1]*K%*%dF
    
    score <- 1/2*((t(dF)*1/F[1]^2)%*%vec(V^2-F)-(2*t(dV)*1/F[1])%*%V)
    I <- 1/2*((t(dF)*1/F[1]^2)%*%dF+2*(t(dV)*1/F[1])%*%dV)
    s <- pinv(I)%*%score
    
    dA <- kronecker.prod(t(A),diag(22))%*%dT+dC+T%*%dA+V[1]*dK+K%*%dV #Check is dC + or *
    dP <-N%*%(kronecker.prod(T%*%P,diag(22))%*%dT-kronecker.prod(K%*%F,diag(22))%*%dK)+kronecker.prod(T,T)%*%dP-kronecker.prod(K,K)%*%dF+dQ
    A <- C+T%*%A+K%*%V
    P <- T%*%P%*%t(T-K%*%Z)+Q
    f <- omega+alpha*s+B*f
    Tt[[i]] <- T; Ct[[i]] <- C; At[,i] <- A; Pt[[i]] <- P; dAt[[i]] <- dA; dPt[[i]] <- dP
    F.SHARK[,i] <- f
    s.SHARK[,i] <- s
  }
  
  #Prediction w/ MGF method
  rw.RV.hat.SHARK[rw+window] <- exp(Z%*%(At[,(rw+window-1)])+1/2*Z%*%Pt[[rw+window-1]]%*%t(Z))
}

#Residuals
res.SHARK <- RV[points+window]-rw.RV.hat.SHARK[points+window]
(MSE.SHARK <- mean(res.SHARK^2))
(QLIKE.SHARK <- mean(RV[points+window]/rw.RV.hat.SHARK[points+window]-log(RV[points+window]/rw.RV.hat.SHARK[points+window])-1))

#Plots
plot(RV[points+window],type="l",xlab="Observations",ylab="RV")
title("S&P500 (from 1997-04-08 to 2013-08-30)")
lines(rw.RV.hat.SHARK[points+window],col="red")

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[points+window],RV[points+window]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[points+window],rw.RV.hat.SHARK[points+window]); colnames(y2) <- c("date","RV")

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

var95 <- -1.52*sqrt(rw.RV.hat.SHARK[points+window]) #0.05 quantile of standardized t with dof=4.17
(violations <- sum(as.vector(snp[points+window])<var95))
(violations/length(var95)*100)

plot(as.vector(snp)[points+window],type="l",ylim=c(-15,10))
lines(var95,col="gold")

VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$uc.LRp
VaRTest(alpha=0.05, actual=as.vector(snp[points+window]), VaR=var95, conf.level = 0.95)$cc.LRp
