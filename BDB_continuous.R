##################################################################################################
##################################################################################################
# Function for implementing Bayesian Dynamic Borrowing Double Adjustment for a continuous outcome
# with example using simulated data

# Authors: Alfredo E. Farjat (alfredo.farjat@bayer.com) / Ming-Dauh (ming-dauh.wang@bayer.com) 
# Last update: 16.07.2025
##################################################################################################
##################################################################################################

# Remove everything from working space
rm(list=ls())

# Loading R libraries
library(mvtnorm)
library(rjags)       
library(R2jags)
library(coda)
library(MCMCvis)
library(overlapping)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Set working directory (results will be saved in this directory)
#setwd("C:/Users/GNBPU/OneDrive - Bayer/Desktop/Bayer/BIC ECA/Bayesian Dynamic Borrowing/Results")


##################################################################################################
##################################################################################################
# Function

BDB_cont=function(Y, Z, TRT, X, A, a, s, za, K, doubleadj, alpha, beta, tau){
  # BDB_cont() applies Bayesian Dynamic Borrowing double-adjustment method for
  # continuous outcomes. The function is designed for the case when the control
  # arm of a randomized study is augmented with external data.
  # 
  # References: 
  # Wang, C., Li, H., Chen, W.C., Lu, N., Tiwari, R., Xu, Y. and Yue, L.Q., 2019. 
  # Propensity score-integrated power prior approach for incorporating real-world 
  # evidence in single-arm clinical studies. 
  # Journal of biopharmaceutical statistics, 29(5), pp.731-748.
  #
  # Farjat, A., Ji, Y., Kaiser, A., Lu, C., Pap, A., Potts, J., King, B., Wang, M. 
  # A Comprehensive Bayesian Double-Adjustment Approach to Dynamic Borrowing of 
  # External Data. Manuscript in progress.
  #  
  # INPUT
  # Y: n-by-1 vector of continuou outcome  
  # Z: n-by-1 vector indicator for external (Z=0) and current data (Z=1)
  # TRT: n-by-1 active treatment indicator (active treatment: 1, control: 0) 
  # X: n-by-p matrix of confounders
  # A: number of subjects intended to be borrowed from external data source
  # a: tuning parameter for arctangent elastic function (a>=1)
  # s: number of strata of propensity score distribution
  # za: confidence level value
  # K: number of samples from posterior distribution
  # doubleadj: logical value for double-adjustment estimation. Default is TRUE. 
  #            If FALSE, only baseline adjustment is computed.
  # alpha: shape hyper-parameter from gamma prior distribution for the precision  
  # beta : rate hyper-parameter from gamma prior distribution for the precision
  # tau: precision for for the mean prior distribution centered at zero
  #
  # OUTPUT
  # N0_t: number of external subjects after trimming
  # n_s0: s-by-1 vector of number of external subjects on each stratum 
  # n_s1: s-by-1 vector of number of current subjects on each stratum 
  # v: s-by-1 vector of overlapping coefficients
  # r: s-by-1 vector of normalized overlapping coefficients
  # ESS: s-by-1 vector of effective number of subjects borrowed from external sources by stratum
  # A_eff: approximate number of subjects borrowed from external source
  # PPP: s-by-1 vector of posterior predictive probability
  # gPPP: s-by-1 vector of transformed PPP used for further discount based on 
  #       outcome differences
  # gamma: s-by-1 vector of baseline discounting factor
  # theta_s: 3-by-3 matrix of mean difference and CrI by stratum
  # theta.m: vector of Bayesian posterior samples
  # theta_o: 3-by-1 vector of mean difference effect (treatment vs controls) and 
  #        credible interval for the mean difference
  # LCrI: length of credible interval
  # acf_results: ACF results for vector of Bayesian posterior samples (theta.m)
  #
  # Printed Outputs
  # figure1: Forest plot of treatment effect by stratum and overall
  # figure2: Plot the overall distribution of propensity scores by group
  # figure3: Plot the overall distribution of propensity scores by stratum
  # figure4: Trace plot of mean difference
  # figure5: Density plot of posterior samples
  # figure6: figure4 and figure5 printed together
  # figure7: Plot of g(PPP) vs PPP for elastic parameter a
  # figure8: Forest plot of overall outcomes for ICA and ECA
  # figure9: Forest plot of ICA and ECA outcomes by stratum
  
  ##################################################################################################
  # Function begins
  # Bayesian model for the calculation of predictive probability of observing 
  # current mean given external data
  modell.string<-"
  model
  {
    for (i in 1:n){
      y[i] ~ dnorm(m,pre)
  }
  m ~ dnorm(0,tau)          # prior for the mean, centered at zero
  pre ~ dgamma(alpha, beta) # prior for the precision
  # Sampling from current data conditioned on external data
  theta ~ dnorm(m,pre*N1)     
  }
  "
  cat(
    "
  model{
    for(i in 1:n){
      phi[i] <- -(log(pre)/2-(y[i]-mu)^2*pre/2)*eta[i] + 10000 # stabilizing constant
      zeros[i] ~ dpois(phi[i])
  }

  # Prior distributions
    mu ~ dnorm(0, tau)
    pre ~ dgamma(alpha, beta)
  }
  ", file="model_REP.txt")
  
  # Define data object
  data<-data.frame(Y, Z, TRT, X)
 
  ######################################################################################################
  # Step 1: Modeling 
  ######################################################################################################
  # Define exposure mechanism
  mylogit<-glm(Z ~ .-Y - TRT, data = data, family = 'binomial')
  
  ######################################################################################################
  # Step 2: Propensity Score Calculation 
  ######################################################################################################
  prd<-predict(mylogit, type='response')
  data<-data.frame(data, prd)
  prd1<-prd[Z == 1] # propensity scores of current patients
  
  ######################################################################################################
  # Step 3: Trimming 
  ######################################################################################################
  # Define data1 and restrict propensity scores to the range of values of the 
  # current patients, this is just trimming
  data1<-data[between(prd, min(prd1), max(prd1)),]
  # Propensity score of patients in the trial
  prd1<-data1$prd[data1$Z==1]
  # Propensity scores of external patients
  prd0<-data1$prd[data1$Z==0]
  # Overall number of subjects after trimming
  n<-nrow(data1) 
  # Number of current subjects
  N1<-length(prd1)
  # Number of external subjects after trimming
  n0<-n-N1
  
  ######################################################################################################
  # Step 4: Stratification 
  ######################################################################################################
  # s, number of strata
  if (s==1) {q=c(0,1)
  } else { 
    q<-quantile(prd1, seq(1:(s-1))/s)
    q<-c(0, q, 1) 
  }
  # Define matrix strat to identify indices of each stratum 
  strat<-matrix(data=rep(FALSE, length(data1$prd)*s), nrow = s, ncol = length(data1$prd))
  for (i in 1:s){ strat[i,]<-between(data1$prd, q[i], q[i+1]) }
  
  # Define variable s_ind
  data1$stratum<-NA
  for(i in 1:s) {data1$stratum[strat[i,]==1]<-i}
  
  ######################################################################################################
  # Step 5: Overlapping 
  ######################################################################################################
  
  v<-rep(0,s)
  for (i in 1:s){
    # overlapping probability of stratum i
    temp_lst<-list(P0_s=data1$prd[data1$Z==0 & data1$stratum==i], P1_s=data1$prd[data1$Z==1 & data1$stratum==i])
    v[i]<-as.numeric(ovmult(x=temp_lst)) 
  }
  
  ######################################################################################################
  # Step 6; Weighting 
  ######################################################################################################
  # Initialize vector of normalized overlapping coefficients 
  r<-rep(0,s)
  for (i in 1:s){ r[i]<-v[i]/sum(v) }
  
  ######################################################################################################
  # Step 7: Calibration by stratum
  ######################################################################################################
  # Estimate predictive posterior probability
  e.m<-rep(0,s) # observed mean in current study by stratum
  for(i in 1:s){ e.m[i]<-mean(data1$Y[data1$Z==1 & data1$TRT==0 & data1$stratum==i]) }
  
  PPP<-rep(0,s)
  for(i in 1:s){
    # jags model  
    modell<-textConnection(modell.string)
    # Initial values of parameters
    inits<-list(m=e.m[i])
    # Data list as input to JAGS
    n0_s<-sum(data1$Z==0 & data1$stratum==i)
    y0_s<-data1$Y[data1$Z==0 & data1$stratum==i]
    N1_s<-sum(data1$Z==1 & data1$TRT==0 & data1$stratum==i) # compare internal against external controls
    dataList<-list("n"=n0_s,"y"=y0_s,"N1"=N1_s, "alpha"=alpha, "beta"=beta, "tau"=tau) 
    # Run MCMC by JAGS (sampling from external data assuming precision from current data)
    jagsfit <- jags.model(modell, data=dataList, inits=inits, n.chains=1, n.adapt=K) 
    # Draw MCMC samples
    out<-coda.samples(model=jagsfit, variable.names = c("theta"), n.iter=K)
    out<-do.call(rbind.data.frame, out)
    # Calculate posterior predictive probability
    PPP[i]<-mean(out$theta<e.m[i])
  }
  
  ######################################################################################################
  # Step 8: Mapping 
  ######################################################################################################
  # a: value of the tuning parameter in the arc-tangent elastic function
  gPPP<-atan(a*sin(PPP*pi))/atan(a) # further discounting by the elastic function 
  
  # In case of baseline adjustment only
  if(doubleadj==FALSE) {gPPP=rep(1,s)} 
  
  ######################################################################################################
  # Step 9: Discounting 
  ######################################################################################################
  n_s0<-rep(0,s)
  for (i in 1:s){ n_s0[i]<-length(data1$prd[data1$Z==0 & data1$stratum==i]) }
  
  n_s1<-rep(0,s)
  for (i in 1:s){ n_s1[i]<-length(data1$prd[data1$Z==1 & data1$stratum==i]) }
  
  gamma<-rep(0,s)
  for (i in 1:s){ gamma[i]<-min(1,A*r[i]/n_s0[i]) }
  
  ######################################################################################################
  # Steps 10-11 Analysis/Summary 
  ######################################################################################################
  # Initialize post.mean.samples object 
  post.mean.samples<-matrix(data=rep(0, s*K), nrow=s, ncol=K)
  theta_hat_s<-matrix(rep(0, 3*s), nrow=3, ncol=s) 
  # Obtain posterior mean samples by strata
  for (i in 1:s){ 
    # Treatment
    # Data list
    n_<-sum(data1$stratum==i & data1$TRT==1)
    y_<-data1$Y[data1$stratum==i & data1$TRT==1]
    eta_<-data1$Z[data1$stratum==i & data1$TRT==1] 
    zeros_<-rep(0,n_) 
    data_REP<-list(n=n_, y=y_, eta=eta_, zeros=zeros_, "alpha"=alpha, "beta"=beta, "tau"=tau)
    # Initial values
    mu_<-mean(y_)
    pre_<-1/var(y_)
    inits_REP<-function(){ list( mu=mu_, pre=pre_) }
    # Fit model by strata
    fit_Ts <- jags.model(file = "model_REP.txt", data = data_REP, inits = inits_REP, n.chains = 1)
    update(fit_Ts, K)
    fit_Ts <- coda.samples(fit_Ts, variable.names = c("mu"), n.iter = K, thin = 1)
    outT<-do.call(rbind.data.frame, fit_Ts)
    
    # Internal and external controls
    # Data list
    n_<-sum(data1$stratum==i & data1$TRT==0)
    y_<-data1$Y[data1$stratum==i & data1$TRT==0]
    eta_<-data1$Z[data1$stratum==i & data1$TRT==0] + (data1$Z[data1$stratum==i & data1$TRT==0]==0)*gamma[i]*gPPP[i]
    zeros_<-rep(0,n_) 
    data_REP<-list(n=n_, y=y_, eta=eta_, zeros=zeros_, "alpha"=alpha, "beta"=beta, "tau"=tau)
    # Initial values
    mu_<-mean(y_)
    pre_<-1/var(y_)
    inits_REP<-function(){ list( mu=mu_, pre=pre_) }
    # Fit model by strata
    fit_Cs <- jags.model(file = "model_REP.txt", data = data_REP, inits = inits_REP, n.chains = 1)
    update(fit_Cs, K)
    fit_Cs <- coda.samples(fit_Cs, variable.names = c("mu"), n.iter = K, thin = 1)
    outC<-do.call(rbind.data.frame, fit_Cs)
    
    post.mean.samples[i,]<-outT$mu - outC$mu # mean difference
    
    # Posterior mean by stratum and 95% CrI
    theta_hat_s[,i]<-c(mean(post.mean.samples[i,]), as.numeric(quantile(post.mean.samples[i,],c(za/2,1-za/2)))) 
  }
  
  theta.m<-apply(post.mean.samples,2,mean)
  
  Estimate<-mean(theta.m) # posterior mean
  
  CrI<-as.numeric(quantile(theta.m,c(za/2,1-za/2))) # credible interval
  
  LCrI<-CrI[2]-CrI[1] # width of credible interval
  
  # Overall estimate and CrI
  theta_o<-c(Estimate, CrI)
  
  ######################################################################################################
  ######################################################################################################
  # Forest plot of treatment effect by stratum and overall
  temp<-'Stratum 1'
  if(s>1){for (i in 2:s){ temp<-c(temp, paste('Stratum', i)) }}
  temp<-c(temp, 'Overall')
  label<-temp
  mean  <- c(theta_hat_s[1,], theta_o[1]) 
  lower <- c(theta_hat_s[2,], theta_o[2]) 
  upper <- c(theta_hat_s[3,], theta_o[3]) 
  
  df <- data.frame(label, mean, lower, upper)
  
  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label))
  
  fp1 <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    geom_pointrange(size=1, lwd=1) + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    geom_text(x=1.2, y=theta_o[1], label=paste(round(theta_o[1],2), ' [', round(theta_o[2],2),', ', round(theta_o[3],2),']', sep = '')) +
    xlab('') + 
    ylab('Posterior Mean (95% CrI)') +
    ggtitle('Treatment Effect by Stratum and Overall') +
    theme_bw()  # use a white background
  
  X11()
  print(fp1)
  ######################################################################################################
  ######################################################################################################
  # Plot the overall distribution of propensity scores by treatment group
  df<-data.frame(x = data1$prd, 
                 Control = c(rep('Current', sum(data1$Z)),rep('External',sum(data1$Z==0))) )
  
  fp2<-ggplot(df, aes(x = x, fill = Control)) +
    geom_density(alpha = 0.7, bw=0.05) +
    labs(title = 'Propensity Score Distribution by Group', 
         x='Propensity Score', y='Density') + theme_minimal()
  
  X11()
  print(fp2)
  ######################################################################################################
  ######################################################################################################
  # Plot the overall distribution of propensity scores by stratum
  df<-data.frame(prd = data1$prd, stratum=data1$stratum,
                 Control = c(rep('Current', sum(data1$Z)),rep('External',sum(data1$Z==0))) )
  fp3<-ggplot(df, aes(x = prd, fill = Control)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ stratum, scales='free') +
    labs(title = "Propensity Score Distribution by Group and Stratum",
         x = "Propensity Score",
         y = "Density",
         fill = "Group") +
    theme_minimal()
  
  X11(width = 20, height=8)
  print(fp3)
  ######################################################################################################
  ######################################################################################################
  # Trace plot of mean difference
  df<-data.frame(iteration=1:K, theta=theta.m)
  fp4<-ggplot(df, aes(x=iteration, y=theta)) +
       geom_line( color="black", size=0.5, alpha=0.9, linetype=1) +
       xlab('Iteration') + 
       ylab('theta') +
       ggtitle('Trace Plot')
  
  X11(width = 20, height=8)
  print(fp4)
  ######################################################################################################
  ######################################################################################################
  # Density of posterior samples
  df<-data.frame(x=theta.m)
  fp5<-ggplot(df, aes(x) ) +
    geom_density(alpha = 0.7) +
    labs(title = 'Density of posterior samples ', 
         x='theta', y='Density')
  X11()
  print(fp5)
  #####################################################################################################
  ######################################################################################################
  
  fp6 <- ggarrange(fp4, fp5, nrow = 2, ncol = 1)
  X11()
  print(fp6)
  
  ######################################################################################################
  ######################################################################################################
  # Plot - Elastic parameter 
  ######################################################################################################
  ######################################################################################################
  tuning_par=function(alpha, x){atan(a*sin(x*pi))/atan(a) }
  x<-seq(from=0, to=1, by=0.001)
  y<-tuning_par(a, x)
  
  # Define data frame
  df<-data.frame(x = x, y = y)
  
  fp7<-ggplot(df, aes(x = x, y = y)) + 
    geom_line() +
    labs(title = paste('Elastic Function, a =',a), 
         x='Posterior Predictive Probability (PPP)', y='g(PPP)')
  X11()
  print(fp7)
  ######################################################################################################
  ######################################################################################################
  # Plot - Outcome Intervals for ICA and ECA
  ######################################################################################################
  ######################################################################################################
  label <- c('ICA', 'ECA')
  #calculate mean, standard deviation, and count
  mean_ICA <- mean(data1$Y[which((data1$Z==1) & (data1$TRT==0))])
  mean_ECA <- mean(data1$Y[which((data1$Z==0) & (data1$TRT==0))])
  sd_ICA <- sd(data1$Y[which((data1$Z==1) & (data1$TRT==0))])
  sd_ECA <- sd(data1$Y[which((data1$Z==0) & (data1$TRT==0))])
  n_ICA <- length(data1$Y[which((data1$Z==1) & (data1$TRT==0))])
  n_ECA <- length(data1$Y[which((data1$Z==0) & (data1$TRT==0))])
  
  #calculate confidence intervals
  t_alpha <- 0.05
  tval_ICA <- qt(1 - t_alpha/2, df = n_ICA - 1)
  tval_ECA <- qt(1 - t_alpha/2, df = n_ECA - 1)
  error_margin_ICA <- tval_ICA * (sd_ICA / sqrt(n_ICA))
  error_margin_ECA <- tval_ECA * (sd_ECA / sqrt(n_ECA))
  
  #get 95% CI bounds
  CI_lower_ICA <- mean_ICA - error_margin_ICA
  CI_upper_ICA <- mean_ICA + error_margin_ICA
  CI_lower_ECA <- mean_ECA - error_margin_ECA
  CI_upper_ECA <- mean_ECA + error_margin_ECA
  
  mean  <- c(mean_ICA, mean_ECA) 
  lower <- c(CI_lower_ICA, CI_lower_ECA) 
  upper <- c(CI_upper_ICA, CI_upper_ECA)
  
  df <- data.frame(label, mean, lower, upper)
  
  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label))
  
  fp8 <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    geom_pointrange(size=1, lwd=1) + 
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab('') + 
    ylab('Posterior Mean (95% CrI)') +
    ggtitle('Overall Outcomes for ICA and ECA') +
    theme_bw()  # use a white background
  
  X11()
  print(fp8)
  ######################################################################################################
  ######################################################################################################
  # Plot - Forest plot of ICA and ECA outcomes by stratum
  ######################################################################################################
  ######################################################################################################
  temp<-'Stratum 1'
  if(s>1){for (i in 2:s){ temp<-c(temp, paste('Stratum', i)) }}
  temp<-c(temp, temp)
  label<-temp
  group<-c(rep('ECA',s), rep('ICA',s))
  ICA_mean <- c()
  ECA_mean <- c()
  ICA_lower <- c()
  ICA_upper <- c()
  ECA_lower <- c()
  ECA_upper <- c()
  t_alpha <- 0.05
  for (strat in 1:s) {
    #subset by stratum
    df_s <- data1[data1$stratum==strat,]
    #calculate mean
    ICA_mean[strat] <- mean(df_s$Y[which((df_s$Z==1) & (df_s$TRT==0))])
    ECA_mean[strat] <- mean(df_s$Y[which((df_s$Z==0) & (df_s$TRT==0))])
    #calculate SD
    ICA_sd <- sd(df_s$Y[which((df_s$Z==1) & (df_s$TRT==0))])
    ECA_sd <- sd(df_s$Y[which((df_s$Z==0) & (df_s$TRT==0))])
    #get counts
    ICA_n <- length(df_s$Y[which((df_s$Z==1) & (df_s$TRT==0))])
    ECA_n <- length(df_s$Y[which((df_s$Z==0) & (df_s$TRT==0))])
    #calculate confidence intervals
    tval_ICA <- qt(1 - t_alpha/2, df = ICA_n - 1)
    tval_ECA <- qt(1 - t_alpha/2, df = ECA_n - 1)
    error_margin_ICA <- tval_ICA * (ICA_sd / sqrt(ICA_n))
    error_margin_ECA <- tval_ECA * (ECA_sd / sqrt(ECA_n))
    #get 95% CI bounds
    ICA_lower[strat] <- ICA_mean[strat] - error_margin_ICA
    ICA_upper[strat] <- ICA_mean[strat] + error_margin_ICA
    ECA_lower[strat] <- ECA_mean[strat] - error_margin_ECA
    ECA_upper[strat] <- ECA_mean[strat] + error_margin_ECA
  }
  
  mean_outcome <- c(ECA_mean, ICA_mean)
  lower <- c(ECA_lower, ICA_lower)
  upper <- c(ECA_upper, ICA_upper)
  
  df <- data.frame(group, label, mean_outcome, lower, upper)
  
  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label)[1:s])
  df$group<-factor(df$group, levels=c('ICA', 'ECA'))
  
  fp9 <- ggplot(data=df, aes(x=label, y=mean_outcome, ymin=lower, ymax=upper, fill=group)) +
    geom_linerange(size=5, position=position_dodge(width = 0.25), colour="lightgrey") +
    geom_point(size=3, shape=21, stroke = 0.25, position=position_dodge(width = 0.25)) +
    scale_x_discrete(name=" ") +
    scale_y_continuous(name="Mean Outcome [95% CrI]") +
    coord_flip() +
    ggtitle('Outcomes for ICA and ECA by Stratum') +
    theme_bw()
  
  X11()
  print(fp9)
  ######################################################################################################
  ######################################################################################################
  # Autocorrelation function for trace plot
  fp10<-acf(theta.m, plot=F)
  
  ######################################################################################################
  ######################################################################################################
  
output<-list(N0_t=n0, n_s0=n_s0, n_s1=n_s1, 
             v=v, r=r, ESS=A*r*gPPP, A_eff=A*sum(r*gPPP), 
             PPP=PPP, gPPP=gPPP, gamma=gamma,
             theta_s = theta_hat_s, theta.m = theta.m,
             theta_o = theta_o, LCrI = LCrI,
             acf_results = fp10)
             
output
}


######################################################################################################
######################################################################################################
# Example: Simulate Data
######################################################################################################
######################################################################################################
# Parameters
N0<-3000 # number of external subjects
N1<-300 # number of current subjects

p<-10 # number of variables
rho<-0.1 # correlation between the p variables 
mu0<-rep(1.1, p) # mean of variables from external data source 
mu1<-rep(1, p) # mean of variables from current study
var0<-0.12  # variance of variables from external data source
var1<-0.1 # variance of variables from current study

# Covariance matrices
cov0<-rep(sqrt(var0*rho),p)%*%t(rep(sqrt(var0*rho),p))+diag(p)*var0*(1-rho)
cov1<-rep(sqrt(var1*rho),p)%*%t(rep(sqrt(var1*rho),p))+diag(p)*var1*(1-rho)

# Set pseudo-random generator seed
#set.seed(3/4*pi)

# Generate covariates
x0<-rmvnorm(N0,mu0,cov0)
x1<-rmvnorm(N1,mu1,cov1)

# Convert first four variables into binary variables
x0[,1:4]<-(x0[,1:4]>0)*1
x1[,1:4]<-(x1[,1:4]>0)*1

# True coefficients 
beta<-c(0.2,0.4,0.1,0.3,1,1,1,1,1,1)

# True treatment effect
delta_treat1=-7
# True mean difference between internal and external controls
delta_treat0=2

TRT0<-rep(0, N0)
TRT1<-rbinom(N1, size=1, prob=2/3)

# Outcome
y0 <- 100 + TRT0 + delta_treat0 + x0%*%beta  + rnorm(N0, mean = 0, sd=8)
y1 <- 100 + TRT1*delta_treat1 + x1%*%beta  + rnorm(N1, mean = 0, sd=5)

# Create data object including covariates (X), data source indicator (Z), and outcome (Y)
Y<-c(y0,y1)
TRT<-c(TRT0, TRT1)
Z<-c(rep(0,N0),rep(1,N1)) # data source indicator 0: external, 1: current
X<-rbind(x0,x1)

# Define data object
data<-data.frame(X, Z, TRT, Y)


######################################################################################################
######################################################################################################
# Run Bayesian Dynamic Borrowing Double Adjustment approach

# Estimate treatment effect from a trial with N1=300 subjects randomized according to a 2:1 ratio 
# treatment vs control and borrow information equivalent to A=100 subjects from a external source 
# of N0=3000 patients

output <- BDB_cont(Y=data$Y, 
         Z=data$Z, 
         TRT=data$TRT, 
         X=data[,c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10')], 
         A=100, 
         a=5, 
         s=3, 
         za=0.05, 
         K=1000,
         doubleadj=TRUE,
         alpha=0.01,
         beta=0.01,
         tau=0.00001)
  