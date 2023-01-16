# Description: HARQ model - (RV)
# Course: SQB7002
# Author: Sarkaaj Singh (S2104193)
# Date: 11/10/2022
# Edits: 

#Importing data
library(HARModel); library(xts); data(SP500RM); RV <- as.vector(SP500RM$RV); RQ <- sqrt(as.vector(SP500RM$RQ))
#Daily RV averaged from past 22 days
RV.m <- c(); for (i in 1:(length(RV)-21)){RV.m[i+21] <- mean(RV[i:(i+21)])}
#Daily RV averaged from past 5 days
RV.w <- c(); for (i in 1:(length(RV)-4)){RV.w[i+4] <- mean(RV[i:(i+4)])}
#OLS, estimation over entire sample size
y <- c(RV[23:length(RV)])
x1 <- c(RV[22:(length(RV)-1)]*RQ[22:(length(RV)-1)])
x2 <- c(RV[22:(length(RV)-1)])
x3 <- c(RV.w[22:(length(RV)-1)])
x4 <- c(RV.m[22:(length(RV)-1)])
X <- cbind(1,x1,x2,x3,x4)
Beta.HARQ <- solve((t(X)%*%X))%*%t(X)%*%y; results.HARQ <- as.vector(Beta.HARQ)

#Prediction
RV.hat.HARQ <- X%*%Beta.HARQ
RV.hat.HARQ[RV.hat.HARQ<0] <- mean(RV.hat.HARQ) #Insanity filter
#Residuals
res.HARQ <- RV[23:length(RV)]-RV.hat.HARQ
(MSE.HARQ <- mean((res.HARQ^2)[2:length(res.HARQ)])) #MSE and QLIKE are computed here from for RV[24:length(RV)] for consistency in model comparison
(QLIKE.HARQ <- mean(RV[24:length(RV)]/RV.hat.HARQ[2:length(RV.hat.HARQ)]-log(RV[24:length(RV)]/RV.hat.HARQ[2:length(RV.hat.HARQ)])-1))

#Estimation error Newey-West
term1 <- 0
for (t in 1:nrow(X)){term1 <- res.HARQ[t]^2*(as.matrix(X[t,])%*%X[t,])+term1}
q <- 5 #Input

outer.sum <- 0
for (v in 1:q){
  inner.sum <- 0
  for (t in (v+1):nrow(X)){
    inner.sum <- as.matrix(X[t,])%*%res.HARQ[t]%*%res.HARQ[t-v]%*%X[(t-v),]+as.matrix(X[(t-v),])%*%res.HARQ[t-v]%*%res.HARQ[t]%*%X[t,]+inner.sum
  }
  outer.sum <- (1-v/(q+1))*inner.sum+outer.sum
}
variance <- solve(t(X)%*%(X))%*%(term1+outer.sum)%*%solve(t(X)%*%(X))
(SE.HARQ <- sqrt(diag(variance)))
Beta.HARQ

2*pnorm(abs(Beta.HARQ/SE.HARQ), 0, 1, lower.tail = FALSE)*100 #p-values of estimates (in percentage)
Beta.HARQ-1.96*SE.HARQ; Beta.HARQ+1.96*SE.HARQ # 95% confidence interval 

#Plots (1200x600)
library(ggplot2)

n <- index(SP500RM)
y1 <- data.frame(n[24:4096],RV[24:4096]); colnames(y1) <- c("date","RV")
y2 <- data.frame(n[24:4096],RV.hat.HARQ[2:4074]); colnames(y2) <- c("date","RV")

ggplot()+
  geom_line(data = y1, aes(date, RV), color = "black")+
  geom_line(data = y2, aes(date, RV), color = "red")+
  labs(y="RV", x="Year")+
  scale_x_date(date_labels = "%Y", date_breaks = "1 year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text =element_text(size = 15))
