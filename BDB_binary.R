######################################################################################################
######################################################################################################
######################################################################################################
# Brief description: Function for implementing Bayesian Dynamic Borrowing approach when the response
# is binary
#
# Author: Alfredo E. Farjat (alfredo.farjat@bayer.com)
#
# Last update: 31.07.2025
# 
######################################################################################################
######################################################################################################
######################################################################################################
# Removes everything from memory
rm(list = ls())

# Libraries
library(overlapping)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(ggpubr)


######################################################################################################
######################################################################################################
# BDB_binary

BDB_binary=function(Y, Z, TRT, X, A, a, s, alpha, beta, za, K, doubleadj){
  # BDB_binary() applies Bayesian Dynamic Borrowing method 
  # for binary outcomes 
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
  # Y: n-by-1 vector of binary outcomes (0s and 1s) 
  # Z: n-by-1 vector indicator for external (Z=0) and current data (Z=1)
  # TRT: n-by-1 active treatment indicator (active treatment: 1, control: 0) 
  # X: n-by-p matrix of confounders
  # A: number of subjects intended to be borrowed from external data source
  # s: number of strata of propensity score distribution
  # a: tuning parameter for arctangent elastic function
  # alpha: shape1 hyper-parameter from prior beta distribution of binomial distribution success parameter p 
  # beta : shape2 hyper-parameter from prior beta distribution of binomial distribution success parameter p
  # za: confidence level value
  # K: number of samples from posterior distribution
  # doubleadj: logical value for double-adjustment estimation. Default is TRUE. 
  #            If FALSE, only baseline adjustment is computed.
  #
  # OUTPUT
  # n_s1: s-by-1 vector with number of current subjects on each stratum 
  # n_s0: s-by-1 vector with number of external subjects on each stratum 
  # n_s11: s-by-1 vector with number of current subjects in treatment arm on each stratum
  # n_s10: s-by-1 vector with number of internal control subjects on each stratum 
  # ne_s11: s-by-1 vector with number of events in current treatment arm by stratum 
  # ne_s10: s-by-1 vector with number of events in current control arm by stratum 
  # ne_s0: s-by-1 vector with number of events in external data source by stratum 
  # r_s11: s-by-3 matrix with event rate of current treatment arm and confidence interval by stratum 
  # r_s10: s-by-3 matrix with event rate of internal control arm and confidence interval by stratum
  # r_s0: s-by-3 matrix with event rate of external controls and confidence interval by stratum
  # r_ica: overall event rate and 95% CrI of internal control arm
  # r_eca: overall event rate and 95% CrI of external control arm
  # v: s-by-1 vector of overlapping coefficients
  # r: s-by-1 vector of normalized overlapping coefficients 
  # ESS: s-by-1 vector of effective number of subjects borrowed from external sources by stratum
  # A_eff: Overall effective number of subjects borrowed from external sources
  # PPP: posterior predictive probability
  # gPPP: transformed posterior predictive probability used for further discounting
  # gamma: s-by-1 vector of power discounting parameter
  # theta_s: s-by-3 matrix of treatment effect estimates by stratum
  # theta_o: 1-by-3 vector of overall treatment effect and corresponding 95% CrI. 
  # The treatment effect is defined as the difference of event rates between treatment 
  # and control arms 
  # LCrI: length of credible interval
  # theta.m: vector of Bayesian posterior samples
  #
  # Printed Outputs
  # figure1: Forest plot of treatment effect by stratum and overall
  # figure2: Plot the overall distribution of propensity scores by group
  # figure3: Plot the overall distribution of propensity scores by stratum
  # figure4: Trace plot of mean difference (not displayed)
  # figure5: Density plot of posterior samples (not displayed)
  # figure6: figure4 and figure5 printed together
  # figure7: Plot of g(PPP) vs PPP for elastic parameter a
  # figure8: Forest plot of overall outcome for ICA and ECA
  # figure9: Forest plot of ICA and ECA outcomes by stratum
  
  ######################################################################################################
  # Function begins
  # Define dataZX, dataset of Z and X
  dataZX<-data.frame(Z, X)
  # Define data object
  data<-data.frame(Y, Z, TRT, X)
  
  # Define exposure mechanism
  mylogit <- glm(Z ~ ., family = binomial, data = dataZX)
  
  ######################################################################################################
  ######################################################################################################
  # Step 2: Propensity Score Calculation 
  ######################################################################################################
  ######################################################################################################
  # Propensity scores
  prd<-predict(mylogit, type='response')
  # Add propensity scores to data object
  data<-data.frame(data, prd)
  # Propensity scores of current patients (patients in the trial, RCT: Z=1, ECA: Z=0)
  prd1<-prd[Z == 1] 
  
  ######################################################################################################
  ######################################################################################################
  # Step 3: Trimming 
  ######################################################################################################
  ######################################################################################################
  # Define data1 and restrict propensity scores to the range of values of the 
  # current patients, this is just trimming
  data1<-data[between(prd, min(prd1), max(prd1)),]
  # Sort data 
  data1<-data1[order(data1$Z, decreasing = TRUE),] 
  # Trimmed Propensity scores
  # Propensity score of patients in the trial
  prd1<-data1$prd[data1$Z==1]
  # Propensity scores of external patients
  prd0<-data1$prd[data1$Z==0]
  #prd0<-prd[Z == 0][between(prd[Z == 0], min(prd1), max(prd1))]  # equivalently
  # Define prd10, as the propensity score for subjects in the control arm of the trial
  prd10<-data1$prd[data1$Z==1 & data1$TRT==0]
  # Define prd11, as the propensity score for the subjects in the treatment arm of the trial
  prd11<-data1$prd[data1$Z==1 & data1$TRT==1]
  
  ######################################################################################################
  ######################################################################################################
  # Step 4: Stratification 
  ######################################################################################################
  ######################################################################################################
  # s,  number of strata
  q<-quantile(prd1, seq(1:(s-1))/s)
  q<-c(0, q, 1)
  n_s1<-rep(0, s)
  n_s0<-rep(0, s)
  n_s10<-rep(0, s)
  n_s11<-rep(0, s)
  s_Z1<-matrix(data=rep(FALSE, length(prd1)*s), nrow=s, ncol=length(prd1))
  s_Z0<-matrix(data=rep(FALSE, length(prd0)*s), nrow=s, ncol=length(prd0))
  s_Z10<-matrix(data=rep(FALSE, length(prd10)*s), nrow=s, ncol=length(prd10))
  s_Z11<-matrix(data=rep(FALSE, length(prd11)*s), nrow=s, ncol=length(prd11))
  for (i in 1:s){
    s_Z1[i,]<-between(prd1, q[i], q[i+1])
    s_Z0[i,]<-between(prd0, q[i], q[i+1])
    s_Z10[i,]<-between(prd10, q[i], q[i+1])
    s_Z11[i,]<-between(prd11, q[i], q[i+1])
    n_s1[i]<-sum(s_Z1[i,]) # number of current subjects on each stratum (both treatment and control)
    n_s0[i]<-sum(s_Z0[i,]) # number of external subjects on each stratum
    n_s10[i]<-sum(s_Z10[i,]) # number of internal control subjects on each stratum
    n_s11[i]<-sum(s_Z11[i,]) # number of current subjects on active treatment on each stratum
  }
  
  # define matrix strat to indentify indices of each stratum 
  strat<-matrix(data=rep(FALSE, length(data1$prd)*s), nrow = s, ncol = length(data1$prd))
  for (i in 1:s){
    strat[i,]<-between(data1$prd, q[i], q[i+1])
  }
  
  # Define variable stratum 
  data1$stratum<-NA
  for(i in 1:s) {data1$stratum[strat[i,]==1]<-i}
  
  ######################################################################################################
  ######################################################################################################
  # Step 5: Overlapping coefficient 
  ######################################################################################################
  ######################################################################################################
  # Calculate vector of overlapping coefficients
  v<-rep(0,s)
  
  for (i in 1:s){
    # Overlapping probability of stratum i
    temp_lst=list(P1_s=prd1[s_Z1[i,]], P0_s=prd0[s_Z0[i,]])
    v[i]<-as.numeric(ovmult(x=temp_lst)) 
  }
  
  ######################################################################################################
  ######################################################################################################
  # Step 6: Weighting 
  ######################################################################################################
  ######################################################################################################
  # Initiate vector of normalized overlapping coefficients 
  r<-rep(0,s)
  for (i in 1:s){
    r[i]<-v[i]/sum(v)
  }
  
  ######################################################################################################
  ######################################################################################################
  # Step 7: Calibration 
  ######################################################################################################
  ######################################################################################################
  
  # number of events by stratum for current treatment 
  ne_s11<-rep(0,s)
  for(i in 1:s){ ne_s11[i]<-sum(data1$Y[data1$Z==1 & data1$TRT==1 & data1$stratum==i]) }
  
  # Estimate predictive posterior probability
  # number of events by stratum for internal control
  ne_s10<-rep(0,s)
  for(i in 1:s){ ne_s10[i]<-sum(data1$Y[data1$Z==1 & data1$TRT==0 & data1$stratum==i]) }
  
  e.m<-ne_s10/n_s10 # observed mean of current control
  
  # number of event by stratum for external data
  ne_s0<-rep(0,s)
  for(i in 1:s){ ne_s0[i]<- sum(data1$Y[data1$Z==0 & data1$stratum==i]) }
  
  # Number of random samples
  N<-10000
  # beta-binomial posterior for external data by stratum
  prb<-matrix(0, nrow = s, ncol=N)
  for(i in 1:s){ prb[i,] <- rbeta(N, shape1 = alpha + ne_s0[i], shape2 = beta + n_s0[i] - ne_s0[i]) }
  
  # Posterior predictive probability (sampling from current study with probabilities from external data and 
  # compare against current study)
  PPP<-rep(0, s)
  for(i in 1:s){ PPP[i]<-mean(rbinom(N, n_s10[i], prb[i,])/n_s10[i] > e.m[i]) }
  
  ######################################################################################################
  ######################################################################################################
  # Step 8: Mapping 
  ######################################################################################################
  ######################################################################################################
  # a: value of the tuning parameter in the arctangent elastic function
  
  # Further discounting by the elastic function 
  gPPP<-atan(a*sin(PPP*pi))/atan(a) # mapped predictive posterior probability
  
  # In case of baseline adjustment only
  if(doubleadj==FALSE) {gPPP=rep(1,s)} 
  
  ######################################################################################################
  ######################################################################################################
  # Step 9: Discounting 
  ######################################################################################################
  ######################################################################################################
  # A: number of subjects intended to be borrowed from external source
  # Define discounting factor
  gamma<-rep(0,s)
  for (i in 1:s){ gamma[i]<-min(1,A*r[i]/n_s0[i]) }
  
  # Number of subjects actually borrowed from external source by stratum
  A_borrow<-rep(0,s)
  for (i in 1:s){ A_borrow[i]<-A*r[i] }
  
  ######################################################################################################
  ######################################################################################################
  # Steps 10-11: Analysis/Summary 
  ######################################################################################################
  ######################################################################################################
  
  ######################################################################################################
  ######################################################################################################
  # Both baseline and outcome adjustment 
  ######################################################################################################
  ######################################################################################################
  # Initialize matrices for later use
  thetaC<-matrix(rep(0, s*K), nrow=s, ncol=K) 
  thetaC_s<-matrix(rep(0, 3*s), nrow=3, ncol=s) 
  thetaT<-matrix(rep(0, s*K), nrow=s, ncol=K) 
  thetaT_s<-matrix(rep(0, 3*s), nrow=3, ncol=s) 
  theta<-matrix(rep(0, s*K), nrow=s, ncol=K) 
  theta_s<-matrix(rep(0, 3*s), nrow=3, ncol=s) 
  for (i in 1:s){
    # Treatment arm
    trial.treat_s<-data1[(data1$Z==1 & data1$TRT==1),][s_Z11[i,],]
    n11e<-sum(trial.treat_s$Y) # number of subjects in trial treatment with events
    N11<-nrow(trial.treat_s) # number of subjects in trial treatment
    # Trial control
    trial.cont_s<-data1[(data1$Z == 1 & data1$TRT==0),][s_Z10[i,],]
    n10e<-sum(trial.cont_s$Y) # number of subjects in trial control with events
    N10<-nrow(trial.cont_s) # number of subjects in trial control
    # External control
    ext.cont_s<-data1[(data1$Z == 0),][s_Z0[i,],]
    n0e<-sum(ext.cont_s$Y) # number of external subjects with events 
    N0<-nrow(ext.cont_s)  # number of external subjects
    # Treatment
    # Sampling from beta-binomial distribution with parameters shape parameters alpha and beta
    thetaT[i,]<-rbeta(K, alpha + n11e, beta + N11-n11e)
    # Mean rate by stratum
    thetaT_s[,i]<-c(mean(thetaT[i,]), as.numeric(quantile(thetaT[i,],c(za/2,1-za/2))))
    # Controls
    # Sampling from beta-binomial distribution with parameters shape parameters alpha and beta
    thetaC[i,]<-rbeta(K, alpha + n10e + n0e*gamma[i]*gPPP[i], beta + (N10-n10e) + (N0-n0e)*gamma[i]*gPPP[i])
    # Mean rate by stratum
    thetaC_s[,i]<-c(mean(thetaC[i,]), as.numeric(quantile(thetaC[i,],c(za/2,1-za/2))))
    # Treatment effect (Treatment - Controls)
    theta[i,]<-thetaT[i,] - thetaC[i,]
    # Mean rate by stratum
    theta_s[,i]<-c(mean(theta[i,]), as.numeric(quantile(theta[i,],c(za/2,1-za/2))))
  }
  # Calculate mean across strata (posterior distribution)
   theta.m<-colMeans(theta)
  # Credible Interval
   CrI<-as.numeric(quantile(theta.m, c(za/2,1-za/2)))
  # Length of credible interval 
   LCrI<-CrI[2] - CrI[1]
   # Overall posterior mean with credible interval
   theta_o<-c(mean(theta.m), CrI)
   
   ######################################################################################################
   ######################################################################################################
   # Current treatment event rate 
   ######################################################################################################
   ######################################################################################################
   # Event rate by stratum and CI (Clopper-Person)
   r_s11<-t(rbind(ne_s11/n_s11 , qbeta(za/2, ne_s11, n_s11 - ne_s11 +1) , qbeta(1-za/2, ne_s11 + 1, n_s11 - ne_s11)))
   
   ######################################################################################################
   ######################################################################################################
   # ICA event rate 
   ######################################################################################################
   ######################################################################################################
   # Event rate by stratum and CI (Clopper-Person)
   r_s10<-t(rbind(ne_s10/n_s10 , qbeta(za/2, ne_s10, n_s10 - ne_s10 +1) , qbeta(1-za/2, ne_s10 + 1, n_s10 - ne_s10)))
   
   # n10e, number of subjects in the trial in control with events
   n10e<-sum(data1$Y[data1$Z==1 & data1$TRT==0])
   # N10, number of subjects in the trial in control
   N10<-length(data1$Y[data1$Z==1 & data1$TRT==0])
   
   # Overall event rate from ICA and confidence interval from beta distribution (Clopper-Person interval)
   R_ICA<-c(n10e/N10 , qbeta(za/2, n10e, N10 - n10e +1) , qbeta(1-za/2, n10e+ 1, N10- n10e) )
   
   ######################################################################################################
   ######################################################################################################
   # ECA event rate 
   ######################################################################################################
   ######################################################################################################
   # Event rate by stratum and CrI
   r_s0<-t(rbind(ne_s0/n_s0 , qbeta(za/2, ne_s0, n_s0 - ne_s0 +1) , qbeta(1-za/2, ne_s0 + 1, n_s0 - ne_s0)))
   
   # n0e, number of subjects in the external in control with events
   n0e<-sum(data1$Y[data1$Z==0 & data1$TRT==0])
   # N0, number of subjects in the external control
   N0<-length(data1$Y[data1$Z==0 & data1$TRT==0])
   
   # Overall event rate from ECA and confidence interval from beta distribution (Clopper-Person interval)
   R_ECA<-c(n0e/N0, qbeta(za/2, n0e, N0 - n0e +1) , qbeta(1-za/2, n0e+ 1, N0- n0e) )
   
  ######################################################################################################
  ######################################################################################################
  # Plot - Forest plot of treatment effect by stratum and overall
  ######################################################################################################
  ######################################################################################################
  temp<-'Stratum 1'
  if(s>1){for (i in 2:s){ temp<-c(temp, paste('Stratum', i)) }}
  temp<-c(temp, 'Overall')
  label<-temp
  mean  <- c(theta_s[1,], theta_o[1]) 
  lower <- c(theta_s[2,], theta_o[2]) 
  upper <- c(theta_s[3,], theta_o[3]) 
  
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
  # Plot - Propensity Score Distributions by Group
  ######################################################################################################
  ######################################################################################################
  # Define data frame 
  df<-data.frame(x = data1$prd, 
                 Control = c(rep('Current', sum(data1$Z)),rep('External',sum(data1$Z==0))) )
  
  fp2<-ggplot(df, aes(x = x, fill = Control)) +
    geom_density(alpha = 0.7, bw=0.05) +
    labs(title = 'Propensity Score Distribution by Group', 
         x='Propensity Score', y='Density')
  X11()
  print(fp2)
  
  
  ######################################################################################################
  ######################################################################################################
  # Plot - Propensity Score Distributions by Stratum
  ######################################################################################################
  ######################################################################################################
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
  
  ######################################################################################################
  ######################################################################################################
  # Density of posterior samples
  df<-data.frame(x=theta.m)
  fp5<-ggplot(df, aes(x) ) +
    geom_density(alpha = 0.7, bw=0.05) +
    labs(title = 'Density of posterior samples ', 
         x='theta', y='Density')

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
  # Plot - Incidence proportion for ICA and ECA
  ######################################################################################################
  ######################################################################################################
  label <- c('ICA', 'ECA')
  mean  <- c(R_ICA[1], R_ECA[1]) 
  lower <- c(R_ICA[2], R_ECA[2]) 
  upper <- c(R_ICA[3], R_ECA[3]) 
  
  df <- data.frame(label, mean, lower, upper)
  
  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label))
  
  fp8 <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    geom_pointrange(size=1, lwd=1) + 
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab('') + 
    ylab('Posterior Mean (95% CrI)') +
    ggtitle('Overall event rate for ICA and ECA') +
    theme_bw()  # use a white background
  
  X11()
  print(fp8)
  
  ######################################################################################################
  ######################################################################################################
  # Plot - Forest plot of ICA and ECA event rates by stratum
  ######################################################################################################
  ######################################################################################################
  temp<-'Stratum 1'
  if(s>1){for (i in 2:s){ temp<-c(temp, paste('Stratum', i)) }}
  temp<-c(temp, temp)
  label<-temp
  group<-c(rep('ECA',s), rep('ICA',s))
  rate  <- c(r_s0[,1], r_s10[,1]) 
  lower <- c(r_s0[,2], r_s10[,2]) 
  upper <- c(r_s0[,3], r_s10[,3]) 
  
  df <- data.frame(group, label, rate, lower, upper)
  
  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels=rev(df$label)[1:s])
  df$group<-factor(df$group, levels=c('ICA', 'ECA'))
  
  fp9 <- ggplot(data=df, aes(x=label, y=rate, ymin=lower, ymax=upper, fill=group)) +
    geom_linerange(size=5, position=position_dodge(width = 0.25), colour="lightgrey") +
    geom_point(size=3, shape=21, stroke = 0.25, position=position_dodge(width = 0.25)) +
    scale_x_discrete(name=" ") +
    scale_y_continuous(name="Event rate [95% CrI]", limits = c(0, 1)) +
    coord_flip() +
    ggtitle('Event rate for ICA and ECA by stratum') +
    theme_bw()  
  X11()
  print(fp9)
  
  ######################################################################################################
  ######################################################################################################
  # Function output
  output<-list(n_s1 = n_s1, n_s0 = n_s0, n_s11 = n_s11, n_s10 = n_s10, 
               ne_s11=ne_s11, ne_s10=ne_s10, ne_s0=ne_s0,
               r_s11 = r_s11, r_s10 = r_s10, r_s0 = r_s0,
               r_ica=R_ICA, r_eca=R_ECA,
               v=v, r=r, ESS=A*r*gPPP, A_eff=sum(A*r*gPPP),
               PPP=PPP, gPPP=gPPP, gamma=gamma, 
               theta_s = theta_s, theta.m = theta.m,
               theta_o = theta_o, LCrI=LCrI)
  
  return(output)
  # function ends
}


######################################################################################################
######################################################################################################
# Simulate data
######################################################################################################
######################################################################################################
# Set pseudo-random generator seed
set.seed(pi)

# Parameters
N0<-3000 # number of external subjects
N1<-750 # number of current subjects

p<-10 # number of variables
rho<-0.1 # correlation between the p variables 
mu0<-rep(1.1, p) # mean of variables from external data source 
mu1<-rep(1, p) # mean of variables from current study
var0<-0.12  # variance of variables from external data source
var1<-0.1 # variance of variables from current study

# Covariance matrices
cov0<-rep(sqrt(var0*rho),p)%*%t(rep(sqrt(var0*rho),p))+diag(p)*var0*(1-rho)
cov1<-rep(sqrt(var1*rho),p)%*%t(rep(sqrt(var1*rho),p))+diag(p)*var1*(1-rho)

# Generate covariates
x0<-rmvnorm(N0,mu0,cov0)
x1<-rmvnorm(N1,mu1,cov1)

# Convert first four variables into binary variables
x0[,1:4]<-(x0[,1:4]>1)*1
x1[,1:4]<-(x1[,1:4]>1)*1

# True coefficients 
beta0<- -10.479 # intercept of linear predictor
beta<-c(0.2,0.4,0.1,0.3,1,1,1,1,1,1) 

# True treatment effect
delta_treat1=4.8
delta_treat10=-0.7
# True mean difference between internal and external controls
delta_treat0=4.5

TRT0<-rep(0, N0)
TRT1<-rbinom(N1, size=1, prob=2/3)

# Outcome
y0<-rbinom(N0, size = 1, prob = plogis(beta0 + TRT0 + delta_treat0 + x0%*%beta + rnorm(N0, mean = 0, sd=0.1) ) )
y1<-rbinom(N1, size = 1, prob = plogis(beta0 + TRT1*delta_treat10 + delta_treat1  + x1%*%beta + rnorm(N1, mean = 0, sd=0.1) ) )

# Create data object including covariates (X), data source indicator (Z), and outcome (Y)
Y<-c(y0,y1)
TRT<-c(TRT0, TRT1)
Z<-c(rep(0,N0),rep(1,N1)) # data source indicator 0: external, 1: current
X<-rbind(x0,x1)

# Define data object
data<-data.frame(X, Z, TRT, Y)

#####################################################################################################
######################################################################################################
# Example with simulated data
######################################################################################################
######################################################################################################
# Variable names for the adjustment
var.names<-c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10') 
             
output<-BDB_binary(Y=data$Y, 
                   Z=data$Z, 
                   TRT=data$TRT, 
                   X = data[,var.names], 
                   A = 255, 
                   s=3,
                   a=1,
                   alpha = 0.5, 
                   beta = 0.5,
                   za=0.05,
                   K=1e5,
                   doubleadj=TRUE)
