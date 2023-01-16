# Description: HAR model
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 27/09/2022
# Edits: 01 Trimming dataset, QLIKE correction.

#Importing data
library(HARModel); library(xts); data(SP500RM); RV <- as.vector(SP500RM$RV)
#Daily RV averaged from past 22 days
RV.m <- c(); for (i in 1:(length(RV)-21)){RV.m[i+21] <- mean(RV[i:(i+21)])}
#Daily RV averaged from past 5 days
RV.w <- c(); for (i in 1:(length(RV)-4)){RV.w[i+4] <- mean(RV[i:(i+4)])}
#OLS, estimation over entire sample size
y <- c(RV[23:length(RV)])
x1 <- c(RV[22:(length(RV)-1)])
x2 <- c(RV.w[22:(length(RV)-1)])
x3 <- c(RV.m[22:(length(RV)-1)])
X <- cbind(1,x1,x2,x3)
Beta.HAR <- solve((t(X)%*%X))%*%t(X)%*%y; results.HAR <- as.vector(Beta.HAR)
#Prediction
RV.hat.HAR <- X%*%Beta.HAR
# #ALT Prediction
# RV.hat.HAR <-c(0)
# for (i in 1:length(22:(length(RV)-1))){
#   RV.hat.HAR[i+22] <- (results.HAR[1]+results.HAR[2]*RV[i+21]+results.HAR[3]*mean(RV[(17+i):(21+i)])
#                        +results.HAR[4]*mean(RV[(i):(21+i)]))
# }
# RV.hat.HAR <- RV.hat.HAR[23:length(RV)]
#Residuals
res.HAR <- RV[23:length(RV)]-RV.hat.HAR
(MSE.HAR <- mean((res.HAR^2)[2:length(res.HAR)])) #MSE and QLIKE are computed here from for RV[24:length(RV)] for consistency in model comparison
(QLIKE.HAR <- mean(RV[24:length(RV)]/RV.hat.HAR[2:length(RV.hat.HAR)]-log(RV[24:length(RV)]/RV.hat.HAR[2:length(RV.hat.HAR)])-1))
#Estimation error standard method
sample.variance.of.error.term <- 1/(length(res.HAR)-4)*(t(res.HAR)%*%res.HAR)
covariance.matrix.of.Beta <- sample.variance.of.error.term[1]*solve(t(X)%*%X)
(SE.HAR <- sqrt(diag(covariance.matrix.of.Beta)))
#Estimation error Newey-West
term1 <- 0
for (t in 1:nrow(X)){term1 <- res.HAR[t]^2*(as.matrix(X[t,])%*%X[t,])+term1}
q <- 5 #Input

outer.sum <- 0
for (v in 1:q){
  inner.sum <- 0
  for (t in (v+1):nrow(X)){
    inner.sum <- as.matrix(X[t,])%*%res.HAR[t]%*%res.HAR[t-v]%*%X[(t-v),]+as.matrix(X[(t-v),])%*%res.HAR[t-v]%*%res.HAR[t]%*%X[t,]+inner.sum
  }
  outer.sum <- (1-v/(q+1))*inner.sum+outer.sum
}
variance <- solve(t(X)%*%(X))%*%(term1+outer.sum)%*%solve(t(X)%*%(X))
(SE.HAR <- sqrt(diag(variance)))
Beta.HAR

2*pnorm(abs(Beta.HAR/SE.HAR), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
Beta.HAR-1.96*SE.HAR; Beta.HAR+1.96*SE.HAR # 95% confidence interval 

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[24:4096],RV[24:4096]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[24:4096],RV.hat.HAR[2:4074]); colnames(y2) <- c("date","RV")

ggplot()+
  geom_line(data = y1, aes(date, RV), color = "black")+
  geom_line(data = y2, aes(date, RV), color = "red")+
  labs(y="RV", x="Year")+
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(size = 15))
################################################################################
