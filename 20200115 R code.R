library(rockchalk)
library(stats)
library(mcmc)
library(dae)

#example 7-1
mu=c(1,2)
sigma=matrix(c(1,0.9,0.9,1),2,2,TRUE)
set.seed(2020)
sam_target=mvrnorm(2500,mu,sigma)
N=5000
start_val=0

alpha <- function(x,y){
  num=exp(-1/2*(y-mu)%*%solve(sigma)%*%(y-mu))
  den=exp(-1/2*(x-mu)%*%solve(sigma)%*%(x-mu))
  min((num/den),1)
}

#Example 7-1 (1)
sam1 <- function(x){
  return(x+c(runif(1,-0.75,0.75),runif(1,-1,1)))
}

result1=matrix(NA,N,2)
result1[1,]=start_val
i=2

for (i in 2:N){
  new1=sam1(result1[i-1,])
  a1=alpha(result1[i-1,],new1)
  u=runif(1)
  if(u<=a1){result1[i,]=new1} 
  else {result1[i,]=result1[i-1,]}
}

dev.new()
plot(sam_target,xlab="x1",ylab="x2",xlim=c(-3,5),ylim=c(-2,6),col="gray")
points(result1,xlab="",ylab="",xlim=c(-3,5),ylim=c(-2,6),col="blue")

#Example 7-1 (2)
sam2 <- function(x) {
  return(x+c(rnorm(1,0,sqrt(0.6)),rnorm(1,0,sqrt(0.4))))
}

result2=matrix(NA,N,2)
result2[1,]=start_val

for (i in 2:N){
  new2=sam2(result2[i-1,])
  a2=alpha(result2[i-1,],new2)
  u=runif(1)
  if(u<=a2){result2[i,]=new2
  } else {result2[i,]=result2[i-1,]}
}

dev.new()
plot(sam_target,xlab="x1",ylab="x2",col="gray")
points(result2,xlab="",ylab="",xlim=c(-3,5),ylim=c(-2,6),col="blue")

#Example 7-1 (3)
D=diag(c(2,2))
c=0.9
f=function(x) 1/(2*pi)*1/sqrt(det(sigma))*exp(-1/2*(x-mu)%*%solve(sigma)%*%(x-mu))
h=function(x) 1/(2*pi)*1/sqrt(det(D))*exp(-1/2*(x-mu)%*%D%*%(x-mu))
result3=matrix(NA,N,2)
result3[1,]=start_val

for (i in 2:N) {
  origin=result3[i-1,]
  accept=FALSE
  while(accept==FALSE) {
    z=c(rmvnorm(c(0,0),D,method="choleski"))
    if (runif(1)<=f(z)/(c*h(z))) {
      new=z
      accept=TRUE
    }
  }
  c1=f(origin)<c*h(origin)
  c2=f(new)<c*h(new)
  if (c1==1 & c2==1) a3=1
  if (c1==0 & c2==1) a3=c*h(origin)/f(origin)
  if (c1==1 & c2==0) a3=1
  if (c1==0 & c2==0) a3=min(((f(new)*h(origin))/(f(origin)*h(new))),1)
  if(runif(1) <= a3){
    result3[i, ] <- new
  } else{result3[i,] <- result3[i-1, ]}
}

dev.new()
plot(sam_target,xlab="x1",ylab="x2",col="gray")
points(result3,xlab="",ylab="",col="blue")

#Example 7-1 (4)
sam4=function(x) {
  return(2*mu-x+c(runif(1,-1,1),runif(1,-1,1)))
}

result4=matrix(NA,N,2)
result4[1,]=start_val

for (i in 2:N) {
  new=sam4(result4[i-1,])
  a=alpha(result4[i-1,],new)
  u=runif(1)
  if (u<=a) {result4[i,]=new
  } else {result4[i,]=result4[i-1,]}
}

dev.new()
plot(sam_target,xlab="x1",ylab="x2",col="gray")
points(result4,xlab="",ylab="",xlim=c(-3,5),ylim=c(-2,6),col="blue")

#Example 7-2
set.seed(1)
n=100
N=5000
AR2=arima.sim(list(order=c(2,0,0),ar=c(1,-0.5)),n,sd=1)
Y=AR2[1:2]
param=matrix(NA,N,3)
y=cbind(AR2[3:n])
wt=cbind(AR2[2:(n-1)],AR2[1:(n-2)])
inv_V_matrix=function(phi) {
  matrix(c(1-phi[2]^2,-phi[1]*(1+phi[2]),-phi[1]*(1+phi[2]),1-phi[2]^2),2,2,TRUE)
}
psi=function(phi,sigma2) {
  inv_V=inv_V_matrix(phi)
   1/sigma2*sqrt(abs(det(inv_V)))*exp(-1/2*sigma2*t(Y)%*%inv_V%*%Y)
}
param[1,]=c(2,0.6,0.3)
i=2

for (i in 2:N) {
  origin_phi=param[i-1,2:3]
  inv_V=inv_V_matrix(origin_phi)
  param[i,1]=1/rgamma(1,n/2,t(Y)%*%inv_V%*%Y+sum((y-wt%*%origin_phi)^2))
  G=t(wt)%*%wt
  phi_mu=solve(G)%*%t(wt)%*%y
  phi_sigma2=param[i,1]*solve(G)
  new_phi=c(rmvnorm(phi_mu,phi_sigma2,method="choleski"))
  alpha=min(psi(new_phi,param[i,1])/psi(origin_phi,param[i,1]),1)
  u=runif(1)
  if (u<=alpha) {param[i,2:3]=new_phi}
  else {param[i,2:3]=origin_phi}
}

param=param[round(0.1*N):N,]
param=cbind(param,sqrt(param[,1]))
result=apply(param,2,function(x) {
  c(mean(x),sd(x),median(x),quantile(x,prob=c(0.025,0.975)),cor(x[2:4501],x[1:4500]))
})
result=t(result)
colnames(result)=c("mean","sd","median","lower","upper","correlation")
result=as.data.frame(result)
result$"Num. SE"=sqrt(c(apply(param,2,function(x) olbm(x,50))))
result=result[-1,]
rownames(result)=c("phi1","phi2","sigma^2")
result=result[,c("mean","Num. SE","sd","median","lower","upper","correlation")]
result=round(result,3)
dev.new()
result


