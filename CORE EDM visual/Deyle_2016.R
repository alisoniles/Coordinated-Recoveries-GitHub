library(Matrix)
library(quantreg)


params <- list(nu1=0.1,nu2=0.07,lambda1=3.2,lambda2=2.9,C1star=0.5,C2star=0.5,
               mu1=0.15,mu2=0.15,kappa1=2.5,kappa2=2.0,Rstar=0.3,k=1.2)
dR <- function(R,C1,C2,P1,P2) R*(1-R/params$k) -
  params$mu1*params$kappa1*(C1*R)/(R+params$Rstar) -
  params$mu2*params$kappa2*(C2*R)/(R+params$Rstar)
dC1 <- function(R,C1,C2,P1,P2) params$mu1*params$kappa1*(C1*R)/(R+params$Rstar) -
  params$nu1*params$lambda1*(P1*C1)/(C1+params$C1star) - params$mu1*C1
dC2 <- function(R,C1,C2,P1,P2) params$mu2*params$kappa2*(C2*R)/(R+params$Rstar) -
  params$nu2*params$lambda2*(P2*C2)/(C2+params$C2star) - params$mu2*C2
dP1 <- function(R,C1,C2,P1,P2) params$nu1*params$lambda1*(P1*C1)/(C1+params$C1star) -
  params$nu1*P1
dP2 <- function(R,C1,C2,P1,P2) params$nu2*params$lambda2*(P2*C2)/(C2+params$C2star) -
  params$nu2*P2
dF <- function(R,C1,C2,P1,P2) c(dR(R,C1,C2,P1,P2), dC1(R,C1,C2,P1,P2),
                                dC2(R,C1,C2,P1,P2), dP1(R,C1,C2,P1,P2), dP2(R,C1,C2,P1,P2))

tmax <- 10000; tau <- 5; dt <- 1/100; burn <- 200
d <- array(data = 0, dim = c(tmax/tau,5))
colnames(d) <-c('R','C1','C2','P1','P2')

Xi <- c(1,.5,.8,.7,.8)
#>>>RUN BURN<<<
idex <- 0
while(idex< burn*(1/dt)){
  idex <- idex+1;
  Xi <- Xi + dt*do.call(dF,as.list(Xi))
}
#>>>RUN MODEL<<<
d[1,] <- Xi
idex <- 0; tdex <- 2
while(idex<tmax*(1/dt)){
  idex <- idex+1;
  Xi <- Xi + dt*do.call(dF,as.list(Xi))
  if( (idex*dt) %% tau == 0 && tdex < tmax/tau){
    d[tdex,] <- Xi
    tdex <- tdex+1
  }
}

targ_col <- 2


Embedding <- c("R","C1","C2","P1") 
Edim <- length(Embedding)

d <- d[,Embedding]

coeff_names <- sapply(colnames(d),function(x) paste("d", 
                                                    colnames(d)[targ_col], "d", x, sep = ""))
                      
block <- cbind(d[2:dim(d)[1],targ_col],d[1:(dim(d)[1]-1),])
norm_consts <- apply(block, 2, function(x) sd(x))
block <- as.data.frame(apply(block, 2, function(x) (x-mean(x))/sd(x)))  

lib <- 1:dim(block)[1] 
pred <- 1:dim(block)[1]
theta <- 8

coeff <- array(0,dim=c(length(pred),Edim)) 
colnames(coeff) <- coeff_names
coeff <- as.data.frame(coeff)

lm_svdsolve <- function(y, x, ws, subset = seq_along(y)){
  x <- x[subset,]
  y <- y[subset]
  ws <- ws[subset]
# prepended column of 1s for constant term in linear model
  A <- cbind(1, x) * ws 
  
  A_svd <- svd(A)
# >>REMOVE SMALL SINGULAR VALUES<<
  s <- A_svd$d
  s_inv <- matrix(0, nrow = dim(x)[2]+1, ncol = dim(x)[2]+1) 
  for(i in seq_along(s))
  {
    if(s[i] >= max(s) * 1e-5) 
      s_inv[i,i] <- 1/s[i]
  }
  coeff <- A_svd$v %*% s_inv %*% t(A_svd$u) %*% (ws * y)
  coeff <- t(coeff)
  
  colnames(coeff) <- c("const",colnames(x))
  return(coeff)
}              



for (ipred in 1:length(pred)){
  
  #target point is excluded from the fitting procedure
  libs = lib[-pred[ipred]]
  
  # >>CALCULATE WEIGHTS<<
  q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]), 
              ncol=Edim, nrow=length(libs), byrow = T)
  distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2)) 
  dbar <- mean(distances)
  Ws <- exp(-theta*distances/dbar)
  
  # >>REGRESS<<
  svd_fit <- lm_svdsolve(block[libs,1],block[libs,2:dim(block)[2]],Ws) 
  coeff[ipred,] <- svd_fit[-1]
}

coeff <- cbind(pred,coeff) 
colnames(coeff)[1] <- "t"

  
par(mfrow = c(1,1))
coefflims = list(c(0,3),c(0,2),c(-1,0),c(-1.1,0)) 
trange <- 1:500
  
plot(coeff[trange,"t"],coeff[trange,"dC1dR"],type="l",col="blue",xlab="time", 
     ylab=expression(partialdiff*C[1] / partialdiff*R), 
     ylim=coefflims[[1]],xlim=range(trange),lwd=2)
abline(a=0 ,b=0 , lty="dashed", col="black", lwd=.5)

plot(coeff[trange,"t"],coeff[trange,"dC1dC2"],type="l",col="red",xlab="time",
  ylab=expression(partialdiff*C[1] / partialdiff*C[2]), ylim=coefflims[[3]],xlim=range(trange),lwd=2)
abline(a=0 ,b=0 , lty="dashed", col="black", lwd=.5)

plot(coeff[trange,"t"],coeff[trange,"dC1dP1"],type="l",col="green",xlab="time", ylab=expression(partialdiff*C[1] / partialdiff*P[1]), ylim=coefflims[[4]],xlim=range(trange),lwd=2)
abline(a=0 ,b=0 , lty="dashed", col="black", lwd=.5)
  
qL <- rq(coeff$dC1dC2 ~ coeff$dC1dR,.05) 
qU <- rq(coeff$dC1dC2 ~ coeff$dC1dR,.95)

plot(coeff[1001:1999, 'dC1dR'], coeff[1001:1999,'dC1dC2'],
         xlab=expression(partialdiff*C[1] / partialdiff*R), 
         ylab=expression(partialdiff*C[1] / partialdiff*C[2]), 
         xlim=c(0,2),ylim=c(-1,0.2),pch=8)
abline(a=qL$coefficients[1],b=qL$coefficients[2],lty=2,col= 'red', lwd = 2)
abline(a=qU$coefficients[1],b=qU$coefficients[2],lty=2,col='red', lwd = 2) 

plot(partdiffs$npgo, partdiffs$srkw, xlim=c(-0.01,.3), ylim=c(-.2,.2),pch=8)
abline(a=qL$coefficients[1],b=qL$coefficients[2],lty=2,col= 'red', lwd = 2)
abline(a=qU$coefficients[1],b=qU$coefficients[2],lty=2,col='red', lwd = 2) 


plot(coeffs$pdo.spr.4, partdiffs$srkw, xlim=c(-3,3), ylim=c(-.4,.25),pch=8)
abline(a=qL$coefficients[1],b=qL$coefficients[2],lty=2,col= 'red', lwd = 2)
abline(a=qU$coefficients[1],b=qU$coefficients[2],lty=2,col='red', lwd = 2) 

plot(partdiffs$dPdo, partdiffs$dSrkw, 
     xlab=expression(partialdiff*rec4 / partialdiff*pdo), 
     ylab=expression(partialdiff*rec4 / partialdiff*srkw), 
     xlim=c(-.25,.2), ylim=c(-.6,.25),pch=8)
abline(a=qL$coefficients[1],b=qL$coefficients[2],lty=2,col= 'red', lwd = 2)
abline(a=qU$coefficients[1],b=qU$coefficients[2],lty=2,col='red', lwd = 2) 

  