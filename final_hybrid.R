library(mvtnorm)
library(tidyverse)
library(survey)
library(doParallel)
library(foreach)
library(doSNOW)
library(parallel)
library(progress)

plogit = function(X){
  return(1/(1+exp(-X)))
}

logit = function(x){
  return(log(x)-log(1-x))
}

logit_inv = function(x){
  return(exp(x)/(1+exp(x)))
}

population = readRDS("population.rds")
SAMPLE = readRDS("SAMPLE.rds")

Hybrid = function(sub_dat, niter = 10000, burnin = 5000, pop_X1_mean = sum(population$X1), pop_X2_mean = sum(population$X2)){
  
  n_unit0 = sum(sub_dat$U==0) # unit respondents
  n_unit1 = sum(sub_dat$U==1) # unit nonrespondents
  N = nrow(population)
  
  # adjust unit nonrespondents weights
  sub_dat$origin_W = sub_dat$W
  sub_dat$W = sub_dat$W*10
  sub_dat$pi = sub_dat$pi/10
  sub_dat[sub_dat$U==1, c('X1', 'X2', 'Y1', 'Y2', 'Y3', 'Y4', 'W', 'pi','Rx1','Ry1','Ry2','Ry4')] = NA
  sub_dat[which(sub_dat$U == 1),]$W = (N - sum(sub_dat$W[which(sub_dat$U == 0)]))/sum(sub_dat$U) 
  sub_dat$pi = 1/sub_dat$W
  
  sub_dat[sub_dat$U==0 & sub_dat$Rx1==1,]$X1 = NA
  sub_dat[sub_dat$U==0 & sub_dat$Ry1==1,]$Y1 = NA
  sub_dat[sub_dat$U==0 & sub_dat$Ry2==1,]$Y2 = NA
  sub_dat[sub_dat$U==0 & sub_dat$Ry4==1,]$Y4 = NA
  
  impdat = mice(sub_dat[,c("X1","X2","Y1","Y2","Y3","Y4",'Rx1','Ry1','Ry2','Ry4')],
                m=1, maxit = 5, seed=1)
  adj_dat = complete(impdat, 1) %>% mutate(U=sub_dat$U, W=sub_dat$W, pi=sub_dat$pi, origin_W=sub_dat$origin_W)
  
  # HT estimator
  pop_X1_sd_HT = sqrt(sum((adj_dat$X1/adj_dat$pi)^2*(1-adj_dat$pi))) 
  pop_X2_sd_HT = sqrt(sum((adj_dat$X2/adj_dat$pi)^2*(1-adj_dat$pi)))
  
  ALPHA = matrix(NA,niter-burnin,2)
  BETA = matrix(NA,niter-burnin,2)
  ETA = matrix(NA, niter-burnin, 3)
  TAU = matrix(NA, niter-burnin, 4)
  PHI = matrix(NA, niter-burnin, 5)
  THETA = matrix(NA, niter-burnin, 6)
  OMEGA.Rx1 = OMEGA.Ry1 = OMEGA.Ry2 = matrix(NA, niter-burnin, 5)
  OMEGA.Ry4 = matrix(NA, niter-burnin, 6)
  X1_impute = matrix(NA,niter-burnin,length(adj_dat$X1)) 
  X2_impute = matrix(NA,niter-burnin,length(adj_dat$X2)) 
  Y1_impute = matrix(NA,niter-burnin,length(adj_dat$Y1)) 
  Y2_impute = matrix(NA,niter-burnin,length(adj_dat$Y2)) 
  Y3_impute = matrix(NA,niter-burnin,length(adj_dat$Y3)) 
  Y4_impute = matrix(NA,niter-burnin,length(adj_dat$Y4))
  
  # iteration
  for (i in 1:niter){
    
    ###### Impute X2 ######
    # IM1-2
    T_X2_hat = rnorm(1, pop_X2_mean, pop_X2_sd_HT) 
    n_2U = floor( (T_X2_hat-sum((adj_dat$X2*adj_dat$W)[which(adj_dat$U == 0)])) / mean(adj_dat$W[which(adj_dat$U==1)]) ) 
    
    # IM3: sample alpha
    m_X2 = glm(X2~W, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    X2_mat = model.matrix(m_X2)
    alpha_X2 = alpha_X2_mod = rmvnorm(1, mean = coef(m_X2), sigma = vcov(m_X2))
    
    # IM4-IM5: impute X2 for U=1
    m_X2_1 = glm(X2~W, data = adj_dat[which(adj_dat$U == 1),], family=binomial())
    X2_mat1 = model.matrix(m_X2_1)
    
    if(n_2U/n_unit1 < 1 && n_2U/n_unit1 > 0){
      u_X2 = logit(n_2U/n_unit1) - mean(X2_mat1 %*% t(alpha_X2))
      alpha_X2_mod[1] = alpha_X2_mod[1] + u_X2
      X2_prob = logit_inv(X2_mat1%*%t(alpha_X2_mod))
      adj_dat$X2[which(adj_dat$U == 1)] = rbinom(n_unit1, 1, prob = X2_prob)
    } else{
      if(n_2U > n_unit1){
        adj_dat$X2[which(adj_dat$U == 1)] = rep(1,n_unit1)
      }
      
      if(n_2U <= 0){
        adj_dat$X2[which(adj_dat$U == 1)] = rep(0,n_unit1)
      }
    }
    
    ###### unit nonresponse for X1 ######
    # IM6-IM7
    T_X1_hat = rnorm(1, pop_X1_mean, pop_X2_sd_HT) 
    n_1U = floor( (T_X1_hat-sum((adj_dat$X1*adj_dat$W)[which(adj_dat$U == 0)])) / mean(adj_dat$W[which(adj_dat$U==1)]) ) 
    
    # IM8: sample beta 
    m_X1 = glm(X1~X2, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    X1_mat = model.matrix(m_X1)
    beta_X1 = beta_X1_mod = rmvnorm(1, mean = coef(m_X1), sigma = vcov(m_X1))
    
    # IM9-IM10: impute X1 for U=1
    m_X1_1 = glm(X1~X2, data = adj_dat[which(adj_dat$U == 1),], family=binomial())
    X1_mat1 = model.matrix(m_X1_1)
    
    if(n_1U/n_unit1 < 1 && n_1U/n_unit1 > 0){
      u_X1 = logit(n_1U/n_unit1) - mean(X1_mat1%*%t(beta_X1))
      beta_X1_mod[1] = beta_X1_mod[1] + u_X1
      X1_prob = logit_inv(X1_mat1%*%t(beta_X1_mod))
      adj_dat$X1[which(adj_dat$U == 1)] = rbinom(n_unit1, 1, prob = X1_prob)
    } else{
      if(n_1U > n_unit1){
        adj_dat$X1[which(adj_dat$U == 1)] = rep(1,n_unit1)
      }
      
      if(n_1U <= 0){
        adj_dat$X1[which(adj_dat$U == 1)] = rep(0,n_unit1)
      }
    }
    
    ###### Item nonresponse for X1 ######
    sub_X = adj_dat[which(adj_dat$Rx1 == 1 & adj_dat$U == 0),]
    n_mis_X = nrow(sub_X)
    pr_X_miss = matrix(0, ncol=2, nrow=n_mis_X)
    colnames(pr_X_miss) = c("0","1")
    
    # P(X1|beta)
    subX_mat_mis = model.matrix(glm(X1~X2, family = binomial(), data = sub_X))
    pi_X = pr_X_miss
    pi_X[,1] = 1-logit_inv(subX_mat_mis%*%t(beta_X1))
    pi_X[,2] = logit_inv(subX_mat_mis%*%t(beta_X1))
    
    # P(Y1|X1...)
    m_Y1 = glm(Y1~X1+X2, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Y1_mat = model.matrix(m_Y1)
    eta = rmvnorm(1,mean = c(t(coef(m_Y1))),sigma = vcov(m_Y1))
    
    pi_Y1 = pr_X_miss
    Y1_mat_mis = model.matrix(glm(Y1~X1+X2, family = binomial(), data = sub_X))
    Y1_mat_mis[,"X1"] = 0
    pi_Y1[,1] = ifelse(sub_X$Y1==1, logit_inv(Y1_mat_mis%*%t(eta)), 1-logit_inv(Y1_mat_mis%*%t(eta)))
    Y1_mat_mis[,"X1"] = 1
    pi_Y1[,2] = ifelse(sub_X$Y1==1, logit_inv(Y1_mat_mis%*%t(eta)), 1-logit_inv(Y1_mat_mis%*%t(eta)))
    
    # P(Y2|X1...)
    m_Y2 = glm(Y2~X1+X2+Y1, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Y2_mat = model.matrix(m_Y2)
    tau = rmvnorm(1,mean = c(t(coef(m_Y2))),sigma = vcov(m_Y2))
    
    pi_Y2 = pr_X_miss
    Y2_mat_mis = model.matrix(glm(Y2~X1+X2+Y1, family = binomial(), data = sub_X))
    Y2_mat_mis[, "X1"] = 0
    pi_Y2[,1] = ifelse(sub_X$Y2==1, logit_inv(Y2_mat_mis%*%t(tau)), 1-logit_inv(Y2_mat_mis%*%t(tau)))
    Y2_mat_mis[, "X1"] = 1
    pi_Y2[,2] = ifelse(sub_X$Y2==1, logit_inv(Y2_mat_mis%*%t(tau)), 1-logit_inv(Y2_mat_mis%*%t(tau)))
    
    ### P(Y3|X1...)
    m_Y3 = lm(Y3~X1+X2+Y1+Y2, data = adj_dat[which(adj_dat$U == 0),])
    Y3_mat = model.matrix(m_Y3)
    phi = as.matrix(as.numeric(coefficients(m_Y3)))
    s2.Y3 = var(adj_dat[which(adj_dat$U == 0),]$Y3)
    
    pi_Y3 = pr_X_miss
    Y3_mat_mis = model.matrix(lm(Y3~X1+X2+Y1+Y2, data = sub_X))
    Y3_mat_mis[, "X1"] = 0
    pi_Y3[,1] = dnorm(sub_X$Y3, mean=Y3_mat_mis%*%phi, sd=sqrt(s2.Y3)) 
    Y3_mat_mis[, "X1"] = 1
    pi_Y3[,2] = dnorm(sub_X$Y3, mean=Y3_mat_mis%*%phi, sd=sqrt(s2.Y3))
    
    ### P(Y4|X1...)
    m_Y4 = lm(Y4~X1+X2+Y1+Y2+Y3, data = adj_dat[which(adj_dat$U == 0),])
    Y4_mat = model.matrix(m_Y4)
    theta = as.matrix(as.numeric(coefficients(m_Y4)))
    s2.Y4 = var(adj_dat[which(adj_dat$U == 0),]$Y4)
    
    pi_Y4 = pr_X_miss
    Y4_mat_mis = model.matrix(lm(Y4~X1+X2+Y1+Y2+Y3, data = sub_X))
    Y4_mat_mis[, "X1"] = 0
    pi_Y4[,1] = dnorm(sub_X$Y4, mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4)) 
    Y4_mat_mis[, "X1"] = 1
    pi_Y4[,2] = dnorm(sub_X$Y4, mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4))
    
    ### P(Ry1|X1...)
    m_Ry1 = glm(Ry1~X1+X2+Y2+Y3, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Ry1_mat = model.matrix(m_Ry1)
    omega.Ry1 = rmvnorm(1,mean = c(t(coef(m_Ry1))),sigma = vcov(m_Ry1))
    
    pi_Ry1 = pr_X_miss
    Ry1_mat_mis = model.matrix(glm(Ry1~X1+X2+Y2+Y3, family = binomial(), data = sub_X))
    Ry1_mat_mis[, "X1"] = 0
    pi_Ry1[,1] = ifelse(sub_X$Ry1==1, logit_inv(Ry1_mat_mis%*%t(omega.Ry1)), 1-logit_inv(Ry1_mat_mis%*%t(omega.Ry1)))
    Ry1_mat_mis[, "X1"] = 1
    pi_Ry1[,2] = ifelse(sub_X$Ry1==1, logit_inv(Ry1_mat_mis%*%t(omega.Ry1)), 1-logit_inv(Ry1_mat_mis%*%t(omega.Ry1)))
    
    ### P(Ry2|X1...)
    m_Ry2 = glm(Ry2~X1+X2+Y1+Y3, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Ry2_mat = model.matrix(m_Ry2)
    omega.Ry2 = rmvnorm(1,mean = c(t(coef(m_Ry2))),sigma = vcov(m_Ry2))
    
    pi_Ry2 = pr_X_miss
    Ry2_mat_mis = model.matrix(glm(Ry2~X1+X2+Y1+Y3, family = binomial(), data = sub_X))
    Ry2_mat_mis[, "X1"] = 0
    pi_Ry2[,1] = ifelse(sub_X$Ry2==1, logit_inv(Ry2_mat_mis%*%t(omega.Ry2)), 1-logit_inv(Ry2_mat_mis%*%t(omega.Ry2)))
    Ry2_mat_mis[, "X1"] = 1
    pi_Ry2[,2] = ifelse(sub_X$Ry2==1, logit_inv(Ry2_mat_mis%*%t(omega.Ry2)), 1-logit_inv(Ry2_mat_mis%*%t(omega.Ry2)))
    
    ### P(Ry4|X1...)
    m_Ry4 = glm(Ry4~X1+X2+Y1+Y2+Y3, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Ry4_mat = model.matrix(m_Ry4)
    omega.Ry4 = rmvnorm(1,mean = c(t(coef(m_Ry4))),sigma = vcov(m_Ry4))
    
    pi_Ry4 = pr_X_miss
    Ry4_mat_mis = model.matrix(glm(Ry4~X1+X2+Y1+Y2+Y3, family = binomial(), data = sub_X))
    Ry4_mat_mis[, "X1"] = 0
    pi_Ry4[,1] = ifelse(sub_X$Ry4==1, logit_inv(Ry4_mat_mis%*%t(omega.Ry4)), 1-logit_inv(Ry4_mat_mis%*%t(omega.Ry4)))
    Ry2_mat_mis[, "X1"] = 1
    pi_Ry4[,2] = ifelse(sub_X$Ry4==1, logit_inv(Ry4_mat_mis%*%t(omega.Ry4)), 1-logit_inv(Ry4_mat_mis%*%t(omega.Ry4)))
    
    # combine & calculate prob
    pr_X_miss = pi_X*pi_Y1*pi_Y2*pi_Y3*pi_Y4*pi_Ry1*pi_Ry2*pi_Ry4
    pr_X_miss_prob = ifelse(is.na(pr_X_miss/rowSums(pr_X_miss)),0,pr_X_miss/rowSums(pr_X_miss))
    adj_dat[which(adj_dat$Rx == 1 & adj_dat$U == 0),]$X1 = rbinom(n_mis_X, 1, pr_X_miss_prob[,2])
    
    ###### unit nonresponse for Y1 ###### 
    m_Y1 = glm(Y1~X1+X2, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Y1_mat = model.matrix(m_Y1)
    eta_Y1 = eta_Y1_mod = rmvnorm(1, mean = coef(m_Y1), sigma = vcov(m_Y1))
    
    m_Y1_1 = glm(Y1~X1+X2, data = adj_dat[which(adj_dat$U == 1),], family=binomial())
    Y1_mat1 = model.matrix(m_Y1_1)
    Y1_prob = logit_inv(Y1_mat1%*%t(eta_Y1_mod))
    adj_dat$Y1[which(adj_dat$U == 1)] = rbinom(n_unit1, 1, prob = Y1_prob)
    
    ###### item nonresponse for Y1 ###### 
    sub_Y1 = adj_dat[which(adj_dat$Ry1 == 1 & adj_dat$U == 0),]
    n_mis_Y1 = nrow(sub_Y1)
    pr_Y1_miss = matrix(0, ncol=2, nrow=n_mis_Y1)
    colnames(pr_Y1_miss) = c("1","2")
    
    ### P(Y1|...)
    subY1_mat_mis = model.matrix(glm(Y1~X1+X2, family = binomial(), data = sub_Y1))
    pi_Y1 = pr_Y1_miss
    pi_Y1[,1] = 1-logit_inv(subY1_mat_mis%*%t(eta))
    pi_Y1[,2] = logit_inv(subY1_mat_mis%*%t(eta))
    
    ### P(Y2|Y1...)
    pi_Y2 = pr_Y1_miss
    Y2_mat_mis = model.matrix(glm(Y2~X1+X2+Y1, family = binomial(), data = sub_Y1))
    Y2_mat_mis[,"Y1"] = 0
    pi_Y2[,1] = ifelse(sub_Y1$Y2==1, logit_inv(Y2_mat_mis%*%t(tau)), 1-logit_inv(Y2_mat_mis%*%t(tau)))
    Y2_mat_mis[,"Y1"] = 1
    pi_Y2[,2] = ifelse(sub_Y1$Y2==1, logit_inv(Y2_mat_mis%*%t(tau)), 1-logit_inv(Y2_mat_mis%*%t(tau)))
    
    ### P(Y3|Y1...)
    pi_Y3 = pr_Y1_miss
    Y3_mat_mis = model.matrix(lm(Y3~X1+X2+Y1+Y2, data = sub_Y1))
    Y3_mat_mis[, "Y1"] = 0
    pi_Y3[,1] = dnorm(sub_Y1$Y3, mean=Y3_mat_mis%*%phi, sd=sqrt(s2.Y3)) 
    Y3_mat_mis[, "Y1"] = 1
    pi_Y3[,2] = dnorm(sub_Y1$Y3, mean=Y3_mat_mis%*%phi, sd=sqrt(s2.Y3))
    
    ### P(Y4|Y1...)
    pi_Y4 = pr_Y1_miss
    Y4_mat_mis = model.matrix(lm(Y4~X1+X2+Y1+Y2+Y3, data = sub_Y1))
    Y4_mat_mis[, "Y1"] = 0
    pi_Y4[,1] = dnorm(sub_Y1$Y4, mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4)) 
    Y4_mat_mis[, "Y1"] = 1
    pi_Y4[,2] = dnorm(sub_Y1$Y4, mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4))
    
    ### P(Rx1|Y1...)
    m_Rx1 = glm(Rx1~X2+Y1+Y2+Y3, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Rx1_mat = model.matrix(m_Rx1)
    omega.Rx1 = rmvnorm(1,mean = c(t(coef(m_Rx1))),sigma = vcov(m_Rx1))
    
    pi_Rx1 = pr_Y1_miss
    Rx1_mat_mis = model.matrix(glm(Rx1~X2+Y1+Y2+Y3, family = binomial(), data = sub_Y1))
    Rx1_mat_mis[, "Y1"] = 0
    pi_Rx1[,1] = ifelse(sub_Y1$Rx1==1, logit_inv(Rx1_mat_mis%*%t(omega.Rx1)), 1-logit_inv(Rx1_mat_mis%*%t(omega.Rx1)))
    Rx1_mat_mis[, "Y1"] = 1
    pi_Rx1[,2] = ifelse(sub_Y1$Rx1==1, logit_inv(Rx1_mat_mis%*%t(omega.Rx1)), 1-logit_inv(Rx1_mat_mis%*%t(omega.Rx1)))
    
    ### P(Ry2|Y1...)
    pi_Ry2 = pr_Y1_miss
    Ry2_mat_mis = model.matrix(glm(Ry2~X1+X2+Y1+Y3, family = binomial(), data = sub_Y1))
    Ry2_mat_mis[, "Y1"] = 0
    pi_Ry2[,1] = ifelse(sub_Y1$Ry2==1, logit_inv(Ry2_mat_mis%*%t(omega.Ry2)), 1-logit_inv(Ry2_mat_mis%*%t(omega.Ry2)))
    Ry2_mat_mis[, "Y1"] = 1
    pi_Ry2[,2] = ifelse(sub_Y1$Ry2==1, logit_inv(Ry2_mat_mis%*%t(omega.Ry2)), 1-logit_inv(Ry2_mat_mis%*%t(omega.Ry2)))
    
    ### P(Ry4|X1...)
    pi_Ry4 = pr_Y1_miss
    Ry4_mat_mis = model.matrix(glm(Ry4~X1+X2+Y1+Y2+Y3, family = binomial(), data = sub_Y1))
    Ry4_mat_mis[, "Y1"] = 0
    pi_Ry4[,1] = ifelse(sub_Y1$Ry4==1, logit_inv(Ry4_mat_mis%*%t(omega.Ry4)), 1-logit_inv(Ry4_mat_mis%*%t(omega.Ry4)))
    Ry4_mat_mis[, "Y1"] = 1
    pi_Ry4[,2] = ifelse(sub_Y1$Ry4==1, logit_inv(Ry4_mat_mis%*%t(omega.Ry4)), 1-logit_inv(Ry4_mat_mis%*%t(omega.Ry4)))
    
    # combine & calculate prob
    pr_Y1_miss = pi_Y1*pi_Y2*pi_Y3*pi_Y4*pi_Rx1*pi_Ry2*pi_Ry4
    pr_Y1_miss_prob = ifelse(is.na(pr_Y1_miss/rowSums(pr_Y1_miss)),0,pr_Y1_miss/rowSums(pr_Y1_miss))
    adj_dat[which(adj_dat$Ry1 == 1 & adj_dat$U == 0),]$Y1 = rbinom(n_mis_Y1, 1, pr_Y1_miss_prob[,2])
    
    ###### unit nonresponse for Y2 ###### 
    m_Y2 = glm(Y2~X1+X2+Y1, data = adj_dat[which(adj_dat$U == 0),], family=binomial())
    Y2_mat = model.matrix(m_Y2)
    tau_Y2 = tau_Y2_mod = rmvnorm(1, mean = coef(m_Y2), sigma = vcov(m_Y2))
    
    m_Y2_1 = glm(Y2~X1+X2+Y1, data = adj_dat[which(adj_dat$U == 1),], family=binomial())
    Y2_mat1 = model.matrix(m_Y2_1)
    Y2_prob = logit_inv(Y2_mat1%*%t(tau_Y2_mod))
    adj_dat$Y2[which(adj_dat$U == 1)] = rbinom(n_unit1, 1, prob = Y2_prob)
    
    ###### Item nonresponse for Y2 ######
    sub_Y2 = adj_dat[which(adj_dat$Ry2 == 1 & adj_dat$U == 0),]
    n_mis_Y2 = nrow(sub_Y2)
    pr_Y2_miss = matrix(0, ncol=2, nrow=n_mis_Y2)
    colnames(pr_Y2_miss) = c("1","2")
    
    ### P(Y2|...)
    subY2_mat_mis = model.matrix(glm(Y2~X1+X2+Y1, family = binomial(), data = sub_Y2))
    pi_Y2 = pr_Y2_miss
    pi_Y2[,1] = 1-logit_inv(subY2_mat_mis%*%t(tau))
    pi_Y2[,2] = logit_inv(subY2_mat_mis%*%t(tau))
    
    ### P(Y3|Y2...)
    pi_Y3 = pr_Y2_miss
    Y3_mat_mis = model.matrix(lm(Y3~X1+X2+Y1+Y2, data = sub_Y2))
    Y3_mat_mis[, "Y2"] = 0
    pi_Y3[,1] = dnorm(sub_Y2$Y3, mean=Y3_mat_mis%*%phi, sd=sqrt(s2.Y3)) 
    Y3_mat_mis[, "Y2"] = 1
    pi_Y3[,2] = dnorm(sub_Y2$Y3, mean=Y3_mat_mis%*%phi, sd=sqrt(s2.Y3))
    
    ### P(Y4|Y2...)
    pi_Y4 = pr_Y2_miss
    Y4_mat_mis = model.matrix(lm(Y4~X1+X2+Y1+Y2+Y3, data = sub_Y2))
    Y4_mat_mis[, "Y2"] = 0
    pi_Y4[,1] = dnorm(sub_Y2$Y4, mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4)) 
    Y4_mat_mis[, "Y2"] = 1
    pi_Y4[,2] = dnorm(sub_Y2$Y4, mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4))
    
    ### P(Rx1|Y2...)
    pi_Rx1 = pr_Y2_miss
    Rx1_mat_mis = model.matrix(glm(Rx1~X2+Y1+Y2+Y3, family = binomial(), data = sub_Y2))
    Rx1_mat_mis[, "Y2"] = 0
    pi_Rx1[,1] = ifelse(sub_Y2$Rx1==1, logit_inv(Rx1_mat_mis%*%t(omega.Rx1)), 1-logit_inv(Rx1_mat_mis%*%t(omega.Rx1)))
    Rx1_mat_mis[, "Y2"] = 1
    pi_Rx1[,2] = ifelse(sub_Y2$Rx1==1, logit_inv(Rx1_mat_mis%*%t(omega.Rx1)), 1-logit_inv(Rx1_mat_mis%*%t(omega.Rx1)))
    
    ### P(Ry1|Y2...)
    pi_Ry1 = pr_Y2_miss
    Ry1_mat_mis = model.matrix(glm(Ry1~X1+X2+Y2+Y3, family = binomial(), data = sub_Y2))
    Ry1_mat_mis[, "Y2"] = 0
    pi_Ry1[,1] = ifelse(sub_Y2$Ry1==1, logit_inv(Ry1_mat_mis%*%t(omega.Ry1)), 1-logit_inv(Ry1_mat_mis%*%t(omega.Ry1)))
    Ry1_mat_mis[, "Y2"] = 1
    pi_Ry1[,2] = ifelse(sub_Y2$Ry1==1, logit_inv(Ry1_mat_mis%*%t(omega.Ry1)), 1-logit_inv(Ry1_mat_mis%*%t(omega.Ry1)))
    
    ### P(Ry4|Y2...)
    pi_Ry4 = pr_Y2_miss
    Ry4_mat_mis = model.matrix(glm(Ry4~X1+X2+Y1+Y2+Y3, family = binomial(), data = sub_Y2))
    Ry4_mat_mis[, "Y2"] = 0
    pi_Ry4[,1] = ifelse(sub_Y2$Ry4==1, logit_inv(Ry4_mat_mis%*%t(omega.Ry4)), 1-logit_inv(Ry4_mat_mis%*%t(omega.Ry4)))
    Ry4_mat_mis[, "Y2"] = 1
    pi_Ry4[,2] = ifelse(sub_Y2$Ry4==1, logit_inv(Ry4_mat_mis%*%t(omega.Ry4)), 1-logit_inv(Ry4_mat_mis%*%t(omega.Ry4)))
    
    # combine & calculate prob
    pr_Y2_miss = pi_Y2*pi_Y3*pi_Y4*pi_Rx1*pi_Ry1*pi_Ry4
    pr_Y2_miss_prob = ifelse(is.na(pr_Y2_miss/rowSums(pr_Y2_miss)),0,pr_Y2_miss/rowSums(pr_Y2_miss))
    adj_dat[which(adj_dat$Ry2 == 1 & adj_dat$U == 0),]$Y2 = rbinom(n_mis_Y2, 1, pr_Y2_miss_prob[,2])
    
    ###### unit nonresponse for Y3 ###### 
    m_Y3 = lm(Y3~X1+X2+Y1+Y2, data = adj_dat[which(adj_dat$U == 0),])
    Y3_mat = model.matrix(m_Y3)
    phi_Y3 = rmvnorm(1, mean = coef(m_Y3), sigma = vcov(m_Y3))
    
    m_Y3_1 = lm(Y3~X1+X2+Y1+Y2, data = adj_dat[which(adj_dat$U == 1),])
    Y3_mat1 = model.matrix(m_Y3_1)
    adj_dat$Y3[which(adj_dat$U == 1)] = rnorm(n_unit1, mean=Y3_mat1%*%t(phi_Y3), sd=sqrt(s2.Y3))
    
    ###### Impute Y4 ######
    m_Y4 = lm(Y4~X1+X2+Y1+Y2+Y3, data = adj_dat[which(adj_dat$U == 0),])
    Y4_mat = model.matrix(m_Y4)
    theta_Y4 = rmvnorm(1, mean = coef(m_Y4), sigma = vcov(m_Y4))
    
    m_Y4_1 = lm(Y4~X1+X2+Y1+Y2+Y3, data = adj_dat[which(adj_dat$U == 1),])
    Y4_mat1 = model.matrix(m_Y4_1)
    adj_dat$Y4[which(adj_dat$U == 1)] = rnorm(n_unit1, mean=Y4_mat1%*%t(theta_Y4), sd=sqrt(s2.Y4))

    sub_Y4 = adj_dat[which(adj_dat$Ry4 == 1 & adj_dat$U == 0),]
    Y4_mat_mis = model.matrix(lm(Y4~X1+X2+Y1+Y2+Y3, data = sub_Y4))
    adj_dat[which(adj_dat$Ry4 == 1 & adj_dat$U == 0),]$Y4 = rnorm(nrow(sub_Y4), mean=Y4_mat_mis%*%theta, sd=sqrt(s2.Y4))
    
    if (i>burnin){
      ALPHA[i-burnin,] = alpha_X2
      BETA[i-burnin,] = beta_X1
      ETA[i-burnin,] = eta
      TAU[i-burnin,] = tau
      PHI[i-burnin,] = t(phi)
      THETA[i-burnin,] = t(theta)
      OMEGA.Rx1[i-burnin,] = omega.Rx1
      OMEGA.Ry1[i-burnin,] = omega.Ry1
      OMEGA.Ry2[i-burnin,] = omega.Ry2
      OMEGA.Ry4[i-burnin,] = omega.Ry4
      X1_impute[i-burnin,] = adj_dat$X1
      X2_impute[i-burnin,] = adj_dat$X2
      Y1_impute[i-burnin,] = adj_dat$Y1
      Y2_impute[i-burnin,] = adj_dat$Y2
      Y3_impute[i-burnin,] = adj_dat$Y3
      Y4_impute[i-burnin,] = adj_dat$Y4
    }
    
  }
  
  # Create L=50 multiple imputation datasets for every 100 posterior samples
  X1L = X1_impute[seq(1,(niter-burnin),100), ] 
  X2L = X2_impute[seq(1,(niter-burnin),100), ] 
  Y1L = Y1_impute[seq(1,(niter-burnin),100), ] 
  Y2L = Y2_impute[seq(1,(niter-burnin),100), ] 
  Y3L = Y3_impute[seq(1,(niter-burnin),100), ] 
  Y4L = Y4_impute[seq(1,(niter-burnin),100), ]
  
  # obtain results
  return (result=list(MI_X1 = X1L, MI_X2 = X2L, MI_Y1=Y1L, MI_Y2=Y2L, MI_Y3=Y3L, MI_Y4=Y4L, 
                      alpha=ALPHA, beta=BETA, eta=ETA, tau=TAU, phi=PHI, theta=THETA,
                      omega.Rx1=OMEGA.Rx1, omega.Ry1=OMEGA.Ry1, omega.Ry2=OMEGA.Ry2, omega.Ry4=OMEGA.Ry4))
}

numCores = 40
cl = makeCluster(numCores)
registerDoSNOW(cl)
iterations = length(SAMPLE)
pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

RESULTS = foreach(i = 1:iterations, .options.snow = opts) %dopar% {
  Hybrid(sub_dat=SAMPLE[[i]])
}
close(pb)
stopCluster(cl)
saveRDS(RESULTS,file = "MI_hybrid.rds")

################# Results ###################
n_sim = length(SAMPLE)
T_est = T_premiss = T_SRS = matrix(NA, nrow=n_sim, ncol=6) 
cond_prob = cond_prob_premiss = matrix(NA, nrow=n_sim, ncol=8) 
joint_prob = joint_prob_premiss = matrix(NA, nrow=n_sim, ncol=4)
cond_prob_given12 = cond_prob_given12_premiss = matrix(NA, nrow=n_sim, ncol=8)

# variance
V_B_T = matrix(NA, nrow=n_sim, ncol=6) 
V_B_prob = matrix(NA, nrow=n_sim, ncol=8) 
V_B_given12 = matrix(NA, nrow=n_sim, ncol=8) 
V_B_joint = matrix(NA, nrow=n_sim, ncol=4) 

V_W_T = matrix(NA, nrow=n_sim, ncol=6) 
V_W_cond = matrix(NA, n_sim, 8) 
V_W_given12 = matrix(NA, n_sim, 8) 
V_W_joint = matrix(NA, n_sim, 4) 

V_T = matrix(NA, nrow=n_sim, ncol=6)
N = nrow(population)
L = 50 

for(i in 1:n_sim){
  sub_dat = SAMPLE[[i]]
  n = nrow(sub_dat)
  sub_dat$origin_W = sub_dat$W
  sub_dat$W = sub_dat$W*10
  sub_dat$pi = sub_dat$pi/10
  
  W = ifelse(sub_dat$U==1, (N-sum(sub_dat$W*(1-sub_dat$U)))/sum(sub_dat$U), sub_dat$W)
  pi = 1/W
  
  X1L = RESULTS[[i]]$MI_X1
  X2L = RESULTS[[i]]$MI_X2
  Y1L = RESULTS[[i]]$MI_Y1
  Y2L = RESULTS[[i]]$MI_Y2
  Y3L = RESULTS[[i]]$MI_Y3
  Y4L = RESULTS[[i]]$MI_Y4
  
  # pre-missing
  T_premiss[i,] = sapply(sub_dat[c("X1", "X2", "Y1", "Y2", "Y3", "Y4")], function(X) sum(X * sub_dat$W)) 
  T_SRS[i,] = sapply(sub_dat[c("X1", "X2", "Y1", "Y2", "Y3", "Y4")], function(X) sum(X * N/n))
  
  premiss_dat = sub_dat %>% select(X1,X2,Y1,Y2,U) %>% 
    mutate(X1=as.factor(X1), X2=as.factor(X2), Y1=as.factor(Y1), Y2=as.factor(Y2)) %>%
    mutate(W = as.numeric(W), pi = as.numeric(pi), fpc = as.numeric(nrow(population)))
  mydesign_premiss = svydesign(id=~1, data = premiss_dat, weight = ~sub_dat$W, fpc=~fpc)
  
  cond_prob_premiss[i,1:2] = svyby(~as.factor(X2), ~as.factor(X1), mydesign_premiss, svymean)[1:2,2] 
  cond_prob_premiss[i,3:4] = svyby(~as.factor(X1), ~as.factor(X2), mydesign_premiss, svymean)[1:2,2] 
  cond_prob_premiss[i,5:6] = svyby(~as.factor(Y2), ~as.factor(Y1), mydesign_premiss, svymean)[1:2,2]  
  cond_prob_premiss[i,7:8] = svyby(~as.factor(Y1), ~as.factor(Y2), mydesign_premiss, svymean)[1:2,2]  
  
  joint_prob_premiss[i,1] = as.numeric(coef(svymean(~I(X1==0 & Y1==0), design=mydesign_premiss))[2]) 
  joint_prob_premiss[i,2] = as.numeric(coef(svymean(~I(X1==1 & Y1==0), design=mydesign_premiss))[2]) 
  joint_prob_premiss[i,3] = as.numeric(coef(svymean(~I(X1==0 & Y1==1), design=mydesign_premiss))[2]) 
  joint_prob_premiss[i,4] = as.numeric(coef(svymean(~I(X1==1 & Y1==1), design=mydesign_premiss))[2]) 
  
  cond_prob_given12_premiss[i, 1:4] = svyby(~as.factor(Y1), ~X1+X2, design=mydesign_premiss, svymean)[1:4,3]
  cond_prob_given12_premiss[i, 5:8] = svyby(~as.factor(Y2), ~X1+X2, design=mydesign_premiss, svymean)[1:4,3]
  
  # store results
  T_EST = matrix(NA, nrow=L, ncol=6)  
  VW_T = matrix(NA, nrow=L, ncol=6)  
  COND_PROB = matrix(NA, nrow=L, ncol=8)
  VW_COND_PROB = matrix(NA, nrow=L, ncol=8)
  JOINT_PROB = matrix(NA, nrow=L, ncol=4)
  VW_JOINT_PROB = matrix(NA, nrow=L, ncol=4)
  GIVEN12 = matrix(NA, L, 8)
  VW_GIVEN12 = matrix(NA, L, 8)
  
  for(j in 1:L){
    IMP_dat = data.frame(id = as.numeric(rownames(sub_dat)),
                         X1 = as.factor(X1L[j,]),
                         X2 = as.factor(X2L[j,]),
                         Y1 = as.factor(Y1L[j,]),
                         Y2 = as.factor(Y2L[j,]),
                         Y3 = as.numeric(Y3L[j,]),
                         Y4 = as.numeric(Y4L[j,]),
                         U = as.factor(sub_dat$U),
                         W = as.numeric(W),
                         pi = as.numeric(pi),
                         fpc = as.numeric(nrow(population)))
    
    mydesign_imp = svydesign(id=~1, data = IMP_dat, weight = ~W) 
    T_EST[j,] = svytotal(~X1+X2+Y1+Y2+Y3+Y4, mydesign_imp)[c(2,4,6,8,9,10)]
    VW_T[j, ] = (SE(svytotal(~X1+X2+Y1+Y2+Y3+Y4, mydesign_imp))^2)[c(2,4,6,8,9,10)]
    
    # cond_prob
    COND_PROB[j,1:2] = svyby(~as.factor(X2), ~as.factor(X1), mydesign_imp, svymean)[1:2,2]  
    COND_PROB[j,3:4] = svyby(~as.factor(X1), ~as.factor(X2), mydesign_imp, svymean)[1:2,2]  
    COND_PROB[j,5:6] = svyby(~as.factor(Y2), ~as.factor(Y1), mydesign_imp, svymean)[1:2,2]  
    COND_PROB[j,7:8] = svyby(~as.factor(Y1), ~as.factor(Y2), mydesign_imp, svymean)[1:2,2]  
    VW_COND_PROB[j,1:2] = (svyby(~as.factor(X2), ~as.factor(X1), mydesign_imp, svymean)[1:2,4])^2 
    VW_COND_PROB[j,3:4] = (svyby(~as.factor(X1), ~as.factor(X2), mydesign_imp, svymean)[1:2,4])^2 
    VW_COND_PROB[j,5:6] = (svyby(~as.factor(Y2), ~as.factor(Y1), mydesign_imp, svymean)[1:2,4])^2 
    VW_COND_PROB[j,7:8] = (svyby(~as.factor(Y1), ~as.factor(Y2), mydesign_imp, svymean)[1:2,4])^2 
    
    # joint_prob
    JOINT_PROB[j,1] = as.numeric(coef(svymean(~I(X1==0 & Y1==0), design=mydesign_imp))[2]) 
    JOINT_PROB[j,2] = as.numeric(coef(svymean(~I(X1==1 & Y1==0), design=mydesign_imp))[2]) 
    JOINT_PROB[j,3] = as.numeric(coef(svymean(~I(X1==0 & Y1==1), design=mydesign_imp))[2]) 
    JOINT_PROB[j,4] = as.numeric(coef(svymean(~I(X1==1 & Y1==1), design=mydesign_imp))[2])
    VW_JOINT_PROB[j,1] = vcov(svymean(~I(X1==0 & Y1==0), design=mydesign_imp))[1]
    VW_JOINT_PROB[j,2] = vcov(svymean(~I(X1==1 & Y1==0), design=mydesign_imp))[1]
    VW_JOINT_PROB[j,3] = vcov(svymean(~I(X1==0 & Y1==1), design=mydesign_imp))[1]
    VW_JOINT_PROB[j,4] = vcov(svymean(~I(X1==1 & Y1==1), design=mydesign_imp))[1]
    
    # cond prob given X1,X2
    GIVEN12[j,1:4] = svyby(~as.factor(Y1), ~X1+X2, design=mydesign_imp, svymean)[1:4,3]
    GIVEN12[j,5:8] = svyby(~as.factor(Y2), ~X1+X2, design=mydesign_imp, svymean)[1:4,3]
    VW_GIVEN12[j,1:4] = (svyby(~as.factor(Y1), ~X1+X2, design=mydesign_imp, svymean)[1:4, 5])^2
    VW_GIVEN12[j,5:8] = (svyby(~as.factor(Y2), ~X1+X2, design=mydesign_imp, svymean)[1:4, 5])^2
  }
  
  T_est[i,] = colMeans(T_EST)
  V_B_T[i,] = apply(T_EST, 2, var)
  V_W_T[i,] = colMeans(VW_T)
  
  cond_prob[i,] = colMeans(COND_PROB)
  V_B_prob[i,] = apply(COND_PROB, 2, var)
  V_W_cond[i,] = colMeans(VW_COND_PROB)
  
  joint_prob[i,] = colMeans(JOINT_PROB)
  V_B_joint[i,] = apply(JOINT_PROB, 2, var)
  V_W_joint[i,] = colMeans(VW_JOINT_PROB)
  
  cond_prob_given12[i,] = colMeans(GIVEN12)
  V_B_given12[i,] = apply(GIVEN12, 2, var)
  V_W_given12[i,] = colMeans(VW_GIVEN12)
}

colMeans(T_est)
colMeans(cond_prob)
colMeans(joint_prob)
colMeans(cond_prob_given12)

cond_prob_true = c(nrow(population %>% filter(X1==0,X2==0))/nrow(population %>% filter(X1==0)),
                   nrow(population %>% filter(X1==1,X2==0))/nrow(population %>% filter(X1==1)), 
                   nrow(population %>% filter(X1==0,X2==0))/nrow(population %>% filter(X2==0)), 
                   nrow(population %>% filter(X1==0,X2==1))/nrow(population %>% filter(X2==1)),
                   nrow(population %>% filter(Y1==0,Y2==0))/nrow(population %>% filter(Y1==0)), 
                   nrow(population %>% filter(Y1==1,Y2==0))/nrow(population %>% filter(Y1==1)), 
                   nrow(population %>% filter(Y1==0,Y2==0))/nrow(population %>% filter(Y2==0)), 
                   nrow(population %>% filter(Y1==0,Y2==1))/nrow(population %>% filter(Y2==1)) 
) 

joint_prob_true = c(nrow(population %>% filter(X1==0,Y1==0))/nrow(population), 
                    nrow(population %>% filter(X1==1,Y1==0))/nrow(population), 
                    nrow(population %>% filter(X1==0,Y1==1))/nrow(population), 
                    nrow(population %>% filter(X1==1,Y1==1))/nrow(population))

cond_prob_given12_true =
  c(nrow(population %>% filter(Y1==0,X1==0,X2==0))/nrow(population %>% filter(X1==0,X2==0)), 
    nrow(population %>% filter(Y1==0,X1==1,X2==0))/nrow(population %>% filter(X1==1,X2==0)), 
    nrow(population %>% filter(Y1==0,X1==0,X2==1))/nrow(population %>% filter(X1==0,X2==1)), 
    nrow(population %>% filter(Y1==0,X1==1,X2==1))/nrow(population %>% filter(X1==1,X2==1)), 
    nrow(population %>% filter(Y2==0,X1==0,X2==0))/nrow(population %>% filter(X1==0,X2==0)), 
    nrow(population %>% filter(Y2==0,X1==1,X2==0))/nrow(population %>% filter(X1==1,X2==0)), 
    nrow(population %>% filter(Y2==0,X1==0,X2==1))/nrow(population %>% filter(X1==0,X2==1)), 
    nrow(population %>% filter(Y2==0,X1==1,X2==1))/nrow(population %>% filter(X1==1,X2==1))  
  )

V_T = V_W_T + (1+1/50)*V_B_T
V_cond = V_W_cond + (1+1/50)*V_B_prob
V_given12 = V_W_given12 + (1+1/50)*V_B_given12
V_joint = V_W_joint+ (1+1/50)*V_B_joint

CI.T1 = 1 - sum(T_est[,1] - 2.009*sqrt(colMeans(V_T)[1]) >  sum(population$X1) | T_est[,1] + 2.009*sqrt((colMeans(V_T)[1]))  <  sum(population$X1))/n_sim
CI.T2 = 1 - sum(T_est[,2] - 2.009*sqrt(colMeans(V_T)[2]) >  sum(population$X2) | T_est[,2] + 2.009*sqrt(colMeans(V_T)[2]) <  sum(population$X2))/n_sim 
CI.T3 = 1 - sum(T_est[,3] - 2.009*sqrt(colMeans(V_T)[3]) >  sum(population$Y1) | T_est[,3] + 2.009*sqrt(colMeans(V_T)[3]) <  sum(population$Y1))/n_sim 
CI.T4 = 1 - sum(T_est[,4] - 2.009*sqrt(colMeans(V_T)[4]) >  sum(population$Y2) | T_est[,4] + 2.009*sqrt(colMeans(V_T)[4]) <  sum(population$Y2))/n_sim 
CI.T5 = 1 - sum(T_est[,5] - 2.009*sqrt(colMeans(V_T)[5]) >  sum(population$Y3) | T_est[,5] + 2.009*sqrt(colMeans(V_T)[5]) <  sum(population$Y3))/n_sim 
CI.T6 = 1 - sum(T_est[,6] - 2.009*sqrt(colMeans(V_T)[6]) >  sum(population$Y4) | T_est[,6] + 2.009*sqrt(colMeans(V_T)[6]) <  sum(population$Y4))/n_sim 

cond_prob_true1 = nrow(population %>% filter(X2==0, X1==0))/nrow(population %>% filter(X1==0))
CI.cond1 = 1 - sum(cond_prob[,1] - 2.009*sqrt(colMeans(V_cond)[1]) >  cond_prob_true1 | cond_prob[,1] + 2.009*sqrt(colMeans(V_cond)[1]) <  cond_prob_true1)/n_sim 

cond_prob_true2 = nrow(population %>% filter(X2==0, X1==1))/nrow(population %>% filter(X1==1))
CI.cond2 = 1 - sum(cond_prob[,2] - 2.009*sqrt(colMeans(V_cond)[2]) >  cond_prob_true2 | cond_prob[,2] + 2.009*sqrt(colMeans(V_cond)[2]) <  cond_prob_true2)/n_sim 

cond_prob_true3 = nrow(population %>% filter(X2==0, X1==0))/nrow(population %>% filter(X2==0))
CI.cond3 = 1 - sum(cond_prob[,3] - 2.009*sqrt(colMeans(V_cond)[3]) >  cond_prob_true3 | cond_prob[,3] + 2.009*sqrt(colMeans(V_cond)[3]) <  cond_prob_true3)/n_sim

cond_prob_true4 = nrow(population %>% filter(X2==1, X1==0))/nrow(population %>% filter(X2==1))
CI.cond4 = 1 - sum(cond_prob[,4] - 2.009*sqrt(colMeans(V_cond)[4]) >  cond_prob_true4 | cond_prob[,4] + 2.009*sqrt(colMeans(V_cond)[4]) <  cond_prob_true4)/n_sim 

cond_prob_true5 = nrow(population %>% filter(Y1==0, Y2==0))/nrow(population %>% filter(Y1==0))
CI.cond5 = 1 - sum(cond_prob[,5] - 2.009*sqrt(colMeans(V_cond)[5]) >  cond_prob_true5 | cond_prob[,5] + 2.009*sqrt(colMeans(V_cond)[5]) <  cond_prob_true5)/n_sim 

cond_prob_true6 = nrow(population %>% filter(Y1==1, Y2==0))/nrow(population %>% filter(Y1==1))
CI.cond6 = 1 - sum(cond_prob[,6] - 2.009*sqrt(colMeans(V_cond)[6]) >  cond_prob_true6 | cond_prob[,6] + 2.009*sqrt(colMeans(V_cond)[6]) <  cond_prob_true6)/n_sim 

cond_prob_true7 = nrow(population %>% filter(Y1==0, Y2==0))/nrow(population %>% filter(Y2==0))
CI.cond7 = 1 - sum(cond_prob[,7] - 2.009*sqrt(colMeans(V_cond)[7]) >  cond_prob_true7 | cond_prob[,7] + 2.009*sqrt(colMeans(V_cond)[7]) <  cond_prob_true7)/n_sim 

cond_prob_true8 = nrow(population %>% filter(Y1==0, Y2==1))/nrow(population %>% filter(Y2==1))
CI.cond8 = 1 - sum(cond_prob[,8] - 2.009*sqrt(colMeans(V_cond)[8]) >  cond_prob_true8 | cond_prob[,8] + 2.009*sqrt(colMeans(V_cond)[8]) <  cond_prob_true8)/n_sim 

CI_given12 = c()
for (k in 1:8){
  p_true = cond_prob_given12_true[k]
  var = colMeans(V_given12)[k]
  CI_given12[k] = 1 - sum(cond_prob_given12[,k] - 2.009*sqrt(var) > p_true | cond_prob_given12[,k] + 2.009*sqrt(var) < p_true)/n_sim
}

CI_joint = c()
joint_true = c(nrow(population%>%filter(X1==0, Y1==0))/N,
               nrow(population%>%filter(X1==1, Y1==0))/N,
               nrow(population%>%filter(X1==0, Y1==1))/N,
               nrow(population%>%filter(X1==1, Y1==1))/N)

for (k in 1:4){
  p_true = joint_true[k]
  var = V_joint[k]
  CI_joint[k] = 1 - sum(joint_prob[,k] - 2.009*sqrt(var) > p_true | joint_prob[,k] + 2.009*sqrt(var) < p_true)/n_sim
}

CI.list = list(CI.T = c(CI.T1, CI.T2, CI.T3, CI.T4, CI.T5, CI.T6), 
               CI.cond = c(CI.cond1, CI.cond2, CI.cond3, CI.cond4, CI.cond5, CI.cond6, CI.cond7, CI.cond8),
               CI.joint = CI_joint, CI.given12 = CI_given12)

VarMDAM.list = list(Premissvar.T = apply(T_premiss, 2, var), 
                    Premissvar.cond = apply(cond_prob_premiss, 2, var),
                    Premissvar.joint = apply(joint_prob_premiss, 2, var),
                    Premissvar.given12 = apply(cond_prob_given12_premiss, 2, var),
                    var.T = apply(T_est, 2, var), var.cond = apply(cond_prob, 2, var),
                    var.joint = apply(joint_prob, 2, var), var.given12 = apply(cond_prob_given12, 2, var),
                    AvgEstVar.T = colMeans(V_T), AvgEstVar.cond = colMeans(V_cond),
                    AvgEstVar.joint = colMeans(V_joint), AvgEstVar.given12 = colMeans(V_given12))

