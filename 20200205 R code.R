library(rstan)
library(brms)

ern=read.csv("C:/data/ERN.csv",stringsAsFactors = FALSE)

model3=brm(ern_mean~anxiety+sex+anxiety*sex,
           data=ern,
           prior=c(set_prior("normal(0,3)",class="b"),
                   set_prior("cauchy(0,2.5)",class="sigma")),
           family=gaussian(),
           chains=4,
           iter=2000,
           warmup=1000,
           seed=20162093)

summary=summary(model3, waic = TRUE)
result1=rbind(as.data.frame(summary$fixed),as.data.frame((summary$spec_pars)))
result1=result1[,c(1,3,4,7,5)]
result1[,c(1,2,3)]=round(result1[,c(1,2,3)],2)
result1[,5]=round(result1[,5],0)
colnames(result1)=c("Post.Mean","95% CI (Lower)","95% CI (Upper)","ESS","R hat")
rownames(result1)=c("Intercept","Anxiety","Sex","Anx*Sex","Res.SD")
waic=WAIC(model3)
waic=as.data.frame(waic$estimates)
loo=LOO(model3)
loo=as.data.frame(loo$estimates)
result2=rbind(waic,loo)
result2=round(result2,2)
colnames(result2)=c("Post.Mean","SE")
rownames(result2)=c("elpd_waic","p_waic","WAIC","elpd_loo","p_loo","LOO")
result2=result2[c(3,6),1,drop=FALSE]

