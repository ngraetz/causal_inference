set.seed(1234)

treatment=c(2745.6, 1697.1, 1656.4, 978, 703.4, 489.1, 430, 334.1, 302.8, 274.7, 274.7, 255, 242.5, 200.7, 198.6, 129.6, 119, 118.3, 115.3, 92.4, 40.6, 32.7, 31.4, 17.5, 7.7, 4.1)
control=c(1202.6, 830.1, 372.4, 345.5, 321.2, 244.3, 163, 147.8, 95, 87, 81.2, 68.5, 47.3, 41.1, 36.6, 29, 28.6, 26.3, 26, 24.4, 21.4, 17.3, 11.5, 4.9, 4.9, 1.0)

subsets=function(n,r,v=1:n){
  if(r<=0) NULL else
    if(r>=n) v[1:n] else
      rbind(cbind(v[1],subsets(n-1,r-1,v[-1])),subsets(n-1,r,v[-1]))
}

treat.effect.samplemean.montecarlo.test.func=function(treated.r,control.r,K){
  # Create vectors for r and Z, and find total number in 
  # experiment and number of treated subjects
  r=c(treated.r,control.r);
  Z=c(rep(1,length(treated.r)),rep(0,length(control.r)));
  N=length(r);
  m=length(treated.r);
  
  # Observed test statistic
  obs.test.stat=mean(r[Z==1])-mean(r[Z==0]);
  
  # Monte Carlo simulatoin 
  montecarlo.test.stat=rep(0,K);
  for(i in 1:K){
    treatedgroup=sample(1:N,m);  # Draw random assignment
    controlgroup=(1:N)[-treatedgroup];
    # Compute test statistic for random assignment
    montecarlo.test.stat[i]=mean(r[treatedgroup])-mean(r[controlgroup]);  
    
  }
  # Monte Carlo p-value is proportion of randomly drawn 
  # test statistics that are >= observed test statistic
  pval=sum(montecarlo.test.stat>=obs.test.stat)/K; 
  # 95% CI for true p-value based on Monte Carlo p-value
  lowerci=pval-1.96*sqrt(pval*(1-pval)/K);
  upperci=pval+1.96*sqrt(pval*(1-pval)/K);
  list(pval=pval,lowerci=lowerci,upperci=upperci);
}

## 2. (b)
treat.effect.samplemean.montecarlo.test.func(treatment,control,10000)

## 2. (c) Compare to t-test p-value
t.test(treatment,control,var.equal=TRUE,alternative="greater")

## 2. (d) Test statistics for difference in variance
treat.effect.samplevariance.montecarlo.test.func=function(treated.r,control.r,K){
  # Create vectors for r and Z, and find total number in 
  # experiment and number of treated subjects
  r=c(treated.r,control.r);
  Z=c(rep(1,length(treated.r)),rep(0,length(control.r)));
  N=length(r);
  m=length(treated.r);
  
  # Observed test statistic
  obs.test.stat=var(r[Z==1])-var(r[Z==0]);
  
  # Monte Carlo simulatoin 
  montecarlo.test.stat=rep(0,K);
  for(i in 1:K){
    treatedgroup=sample(1:N,m);  # Draw random assignment
    controlgroup=(1:N)[-treatedgroup];
    # Compute test statistic for random assignment
    montecarlo.test.stat[i]=var(r[treatedgroup])-var(r[controlgroup]);  
    
  }
  # Monte Carlo p-value is proportion of randomly drawn 
  # test statistics that are >= observed test statistic
  pval=sum(montecarlo.test.stat>=obs.test.stat)/K; 
  # 95% CI for true p-value based on Monte Carlo p-value
  lowerci=pval-1.96*sqrt(pval*(1-pval)/K);
  upperci=pval+1.96*sqrt(pval*(1-pval)/K);
  list(pval=pval,lowerci=lowerci,upperci=upperci);
}
set.seed(1234)
treat.effect.samplevariance.montecarlo.test.func(treatment,control,10000)

## 2. (e) Look at distribution of treatment vs. control with boxplots.
plot <- data.table(Rainfall=c(treatment,control),
                   Group=c(rep('Treatment',length(treatment)), rep('Control',length(control))),
                   Type=rep('Rainfall', length(c(treatment,control))))
plot_log <- copy(plot)
plot_log[, Rainfall := log(Rainfall)]
plot_log[, Type := 'Log Rainfall']
plot <- rbind(plot, plot_log)
boxes <- ggplot() + 
  geom_boxplot(data=plot,
               aes(x=Group,
                   y=Rainfall)) +
  facet_wrap(~Type, scales='free_y')
boxes

# Additive treatment effect model implies that the distribution of observed outcomes among treated is the same as among control.
# Our boxplot of rainfall in normal space suggest this is not the case (much more dispersion and extreme outliers among treated). 
# Roughly equal dispersion in log(rainfall). 

## 2. (f) Compare Fisher's sharp null additive model with multiplicative treatment model.
mult_fisher <- t.test(log(treatment),log(control),var.equal=TRUE,alternative="greater")
mult_fisher_effect <- mult_test$estimate[1] - mult_test$estimate
mult_fisher_p <- mult_test$p.value

## 2. (g) Wilcoxon, multiplicative
mult_wilcoxon <- wilcox.test(log(treatment),log(control),conf.int=TRUE)
effect_ci <- mult_wilcoxon$conf.int
mult_wilcoxon

## 2. (h) Calculate multiplicative effect in normal space (most interpretable) by exponentiating additive effect in log space from our multiplicative model test.


