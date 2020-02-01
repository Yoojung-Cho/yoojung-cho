library(Biobase)
library(rBeta2009)
library(TailRank)
library(ggplot2)
library(graphics)

n=16
alpha=2
beta=4
m=500
k=10
set.seed(19971019)
set.seed(2137)

#Gibbs Sampling
y=matrix(NA,m,k+1)
x=matrix(NA,m,k)
y[,1]=round(runif(m),3)

for (i in 1:m) {
  for (j in 1:k) {
    x[i,j]=rbinom(1,n,y[i,j])
    y[i,j+1]=rbeta(1,(x[i,j]+alpha),(n-x[i,j]+beta))
  }
}
gibbs=x[,k]

#Sampling from Target density
target=rbb(m,n,alpha,beta)

#Data process
result=append(gibbs,target)
group=c(rep("Gibbs",500),rep("Target",500))
data=data.frame(result,group)
colnames(data)=c("sample","group")

#Histogram
color=c(rep(c("black","white"),16))
ggplot(data,aes(x=sample,fill=group))+
  geom_histogram(binwidth=1,position="dodge",color="black")+
  scale_x_continuous(breaks=seq(0,16,1))+
  scale_y_continuous(breaks=seq(0,70,10),limits=c(0,70))+
  scale_fill_manual(values=color)+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  xlab("")+ylab("")

