# Description: SHARK model using exact MLE (uses log(RV)) using KF from Delle Monache and Corsi 2021
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 28/09/2022
# Edits: 03 Trimming dataset, QLIKE correction.

#Importing data
library(HARModel); library(xts); data(SP500RM); RV <- as.vector(SP500RM$RV); lRV <- log(RV); RQ <- as.vector(SP500RM$RQ)
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
#Load "pracma" for calculation of pseudoinverse of Fisher information matrix
library("pracma")
#Load "matrixcal" for commutation matrix
library("matrixcalc"); C_n_n <- commutation.matrix(22, 22); N <- diag(484)+C_n_n
#Load "fastmatrix" for fast kronecker product calculations
library("fastmatrix")
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
  cat("Parameters:", v); S <- 0;
  for (i in 23:length(RV)){
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
    S <- 1/2*log(2*pi)+1/2*(log(F)+V*(1/F)*V)+S #should technically be 22/2*log(2*pi) but we use 1/2*log(2*pi) to ensure consistent LL values accross models for comparison
  }
  cat(", Objective:", S, "\n"); return(S)
}
start <- Sys.time()
Like(c(-0.020577134, -0.008657636,  0.015427116, -0.008596719, -0.026981821))
Sys.time()-start

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
# start <- Sys.time()
# result_SHARK_BFGS <- suppressWarnings(optim(c(0,0,0,0,0),gr = function(...) grad_forward(...),Like, method = "BFGS", hessian = TRUE))
# (result_SHARK_BFGS$time <- Sys.time()-start)

start <- Sys.time()
result_SHARK_BFGS_full <- optim(c(0,0,0,0,0),Like, method = "BFGS",hessian = TRUE)
(result_SHARK_BFGS_full$time <- Sys.time()-start)

# start <- Sys.time()
# result_SHARK_NM <- optim(c(0,0,0,0,0),Like,method = "Nelder-Mead",control = c(reltol = 1e-6))
# (result_SHARK_NM$time <- Sys.time()-start)

# start <- Sys.time()
# result_SHARK_NM_full <- optim(c(0,0,0,0,0),Like,method = "Nelder-Mead")
# (result_SHARK_NM_full$time <- Sys.time()-start)
###############################################################################
results <- result_SHARK_BFGS_full$par
###############################################################################
#Computing TVP using GAS function
#Initialization of parameters
omega <- 0; alpha <- c(results[1:5]); B <- 1
F.SHARK <- matrix(0,nrow=5,ncol=length(lRV))
s.SHARK <- matrix(0,nrow=5,ncol=length(lRV))
At <- matrix(0,nrow=22,ncol=length(lRV)); Pt <- list()
F.SHARK[,22] <- c(Beta[1],Beta[2],Beta[3],Beta[4],log(var(lRV)))
f <-  c(Beta[1],Beta[2],Beta[3],Beta[4],log(var(lRV)))
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

for (i in 23:length(RV)){
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

#Prediction of lRV
lRV.hat.SHARK <-c(0)
for (i in 1:length(23:(length(RV)-1))){
  lRV.hat.SHARK[i+23] <- (F.SHARK[1,(i+22)]+F.SHARK[2,(i+22)]*lRV[i+22]+F.SHARK[3,(i+22)]*mean(lRV[(18+i):(22+i)])
                          +F.SHARK[4,(i+22)]*mean(lRV[(01+i):(22+i)]))
}

#Residuals of lRV
res.SHARK <- lRV[24:length(RV)]-lRV.hat.SHARK[24:length(RV)]
(MSE.SHARK <- mean(res.SHARK^2))

par(mfrow=c(2,1),mai=c(0.15,0.4,0.1,0.1),mgp=c(0.1,0.5,0.1))
plot(res.SHARK,type="l",xlab="Observations",ylab="Residuals",main="SHARK")
acf(res.SHARK)

#Prediction w/ MGF method
RV.hat.SHARK <-c(0)
for (i in 1:length(23:(length(RV)-1))){
  RV.hat.SHARK[i+23] <- exp(Z%*%(At[,(i+22)])+1/2*Z%*%Pt[[i+22]]%*%t(Z))
}

#Residuals
res.SHARK <- RV[24:length(RV)]-RV.hat.SHARK[24:length(RV)]
(MSE.SHARK <- mean(res.SHARK^2)) 
(QLIKE.SHARK <- mean(RV[24:length(RV)]/RV.hat.SHARK[24:length(RV)]-log(RV[24:length(RV)]/RV.hat.SHARK[24:length(RV)])-1))

#Estimation error
variance <- diag(solve(result_SHARK_BFGS_full$hessian))
(SE.SHARK <- sqrt(variance))

2*pnorm(abs(results/SE.SHARK), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
results-1.96*SE.SHARK; results+1.96*SE.SHARK # 95% confidence interval 

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[24:4096],RV[24:4096]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[24:4096],RV.hat.SHARK[24:4096]); colnames(y2) <- c("date","RV")

ggplot()+
  geom_line(data = y1, aes(date, RV), color = "black")+
  geom_line(data = y2, aes(date, RV), color = "red")+
  labs(y="RV", x="Year")+
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(size = 15))

