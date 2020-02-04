library(rstan)
library(brms)

ern=read.csv("C:/data/ERN.csv",stringsAsFactors = FALSE)

#model 1
model1=brm(ern_mean~anxiety,
           data=ern,
           prior=c(set_prior("normal(0,3)",class="b"),
                   set_prior("cauchy(0,2.5)",class="sigma")),
           family=gaussian(),
           chains=4,
           iter=2000,
           warmup=1000,
           seed=20162093)
summary1=summary(model1,waic=TRUE)
result1=rbind(as.data.frame(summary1$fixed),as.data.frame((summary1$spec_pars)))
result1=result1[,c(1,3,4,7,5)]
result1[,c(1,2,3)]=round(result1[,c(1,2,3)],2)
result1[,5]=round(result1[,5],0)
colnames(result1)=c("Post.Mean","95% CI (Lower)","95% CI (Upper)","ESS","R hat")
rownames(result1)=c("Intercept","Anxiety","Res.SD")
waic1=WAIC(model1)
waic1=as.data.frame(waic1$estimates)
loo1=LOO(model1)
loo1=as.data.frame(loo1$estimates)

#model 2
model2=brm(ern_mean~anxiety+sex,
           data=ern,
           prior=c(set_prior("normal(0,3)",class="b"),
                   set_prior("cauchy(0,2.5)",class="sigma")),
           family=gaussian(),
           chains=4,
           iter=2000,
           warmup=1000,
           seed=20162093)
summary2=summary(model2,waic=TRUE)
result2=rbind(as.data.frame(summary2$fixed),as.data.frame((summary2$spec_pars)))
result2=result2[,c(1,3,4,7,5)]
result2[,c(1,2,3)]=round(result2[,c(1,2,3)],2)
result2[,5]=round(result2[,5],0)
colnames(result2)=c("Post.Mean","95% CI (Lower)","95% CI (Upper)","ESS","R hat")
rownames(result2)=c("Intercept","Anxiety","Sex","Res.SD")
waic2=WAIC(model2)
waic2=as.data.frame(waic2$estimates)
loo2=LOO(model2)
loo2=as.data.frame(loo2$estimates)

#model 3
model3=brm(ern_mean~anxiety+sex+anxiety*sex,
           data=ern,
           prior=c(set_prior("normal(0,3)",class="b"),
                   set_prior("cauchy(0,2.5)",class="sigma")),
           family=gaussian(),
           chains=4,
           iter=2000,
           warmup=1000,
           seed=20162093)
summary3=summary(model3,waic=TRUE)
result3=rbind(as.data.frame(summary3$fixed),as.data.frame((summary3$spec_pars)))
result3=result3[,c(1,3,4,7,5)]
result3[,c(1,2,3)]=round(result3[,c(1,2,3)],2)
result3[,5]=round(result3[,5],0)
colnames(result3)=c("Post.Mean","95% CI (Lower)","95% CI (Upper)","ESS","R hat")
rownames(result3)=c("Intercept","Anxiety","Sex","Anx*Sex","Res.SD")
waic3=WAIC(model3)
waic3=as.data.frame(waic3$estimates)
loo3=LOO(model3)
loo3=as.data.frame(loo3$estimates)

#comparing models
model.fit1=rbind(waic1[3,1],loo1[3,1])
model.fit2=rbind(waic2[3,1],loo2[3,1])
model.fit3=rbind(waic3[3,1],loo3[3,1])
model.fit=cbind(model.fit1,model.fit2,model.fit3)
model.fit=round(model.fit,2)
colnames(model.fit)=c("Post.Mean for model1","Post.Mean for model2","Post.Mean for model3")
rownames(model.fit)=c("WAIC","LOO")

#figure 4
model1.samp=posterior_samples(model1)
model1.ll=quantile(model1.samp$b_anxiety,0.025)
model1.ul=quantile(model1.samp$b_anxiety,0.975)
fig4=ggplot(data = model1.samp,aes(x=b_anxiety)) +
  geom_density()
d4=ggplot_build(fig4)$data[[1]]
fig4=fig4+geom_area(data=subset(d4,x>model1.ll&x<model1.ul),
                   aes(x=x,y=y),fill="gray90",color="black") +
     geom_vline(xintercept=mean(model1.samp$b_anxiety),size=1.2) +
     xlim(c(-1, 1)) +
     ylab("Density \n") +
     xlab("\n Anxiety Slope Value") +
     theme_bw() +
     theme(axis.text=element_text(size=10),
           axis.title=element_text(size=12))
fig4

#figure 5
model1.samp.gg=model1.samp[1:25, 1:2]
model1.samp.gg$iternum=seq(1,length(model1.samp.gg$b_Intercept))
fig5=ggplot(data=ern,aes(x=anxiety,y=ern_mean)) +
     geom_point(size=1) +
     geom_abline(data=model1.samp.gg,
                 aes(intercept=b_Intercept,
                     slope=b_anxiety,
                     group=iternum),alpha=0.5) +
     xlab("\n Standardized Trait Anxiety") +
     ylab("Error Related Negativity \n") +
     theme(axis.text=element_text(size=1),
           axis.title=element_text(size=12)) +
     theme_bw()
fig5

#figure 6
expect_val=fitted(model1,summary=TRUE)
colnames(expect_val)=c("Estimate","Est.Error","low","high")
predict_val=predict(model1, summary=TRUE)
colnames(predict_val)=c("Estimate","Est.Error","low","high")
data1=as.data.frame(cbind(ern=standata(model1)$Y,anxiety=standata(model1)$X,expect_val))
data2=as.data.frame(cbind(ern=standata(model1)$Y,anxiety=standata(model1)$X,predict_val))
fig6=ggplot(data2,aes(x=anxiety,y=ern))+
     geom_point()+
     geom_ribbon(data=data2,aes(x=anxiety,ymin=low,ymax=high),fill="gray80",alpha=0.5)+
     geom_ribbon(data=data1,aes(x=anxiety,ymin=low,ymax=high),fill="gray60",alpha=0.5)+
     geom_ribbon(data=data1,aes(x=anxiety,ymin=Estimate-0.01,ymax=Estimate+0.01),fill="black")+
     xlab("\n Standardized Trait Anxiety")+
     ylab("Error Related Negativity \n")+
     theme_bw()+
     theme(axis.text=element_text(size=10),
           axis.title=element_text(size=12))
fig6










