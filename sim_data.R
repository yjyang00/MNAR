library(tidyverse)
library(survey)

plogit = function(X){
  return(1/(1+exp(-X)))
}

logit = function(x){
  return(log(x)-log(1-x))
}

logit_inv = function(x){
  return(exp(x)/(1+exp(x)))
}

micro = read.csv("ACSPUMS1Y2022_2024-01-17T121210.csv")

# generate population
gen_pop = function(dat=micro, nu = -1.2, alpha = c(0.06, -2.0, -0.00023792),
                   beta=c(0.2, 0.4, -2.0), eta=c(0.2, 0.3, 0.1),
                   tau = c(0.2, 0.4, 0.4, 0.1),  phi = c(0.4, 1.2, -0.9, 0.1, 0.2),
                   theta=c(0.4, 1.2, -0.9, 0.1, -0.1, 0.1)){

  N = nrow(dat)
  W = dat$PWGTP
  pi = 1/W

  pi_U = plogit(nu)
  U = rbinom(N, 1, p=pi_U)

  X2 = U
  pi_2 = plogit(model.matrix(X2~U+W)%*%alpha)
  X2 = rbinom(N, 1, p=pi_2)
  
  X1 = X2
  pi_1 = plogit(model.matrix(X1~X2+U)%*%beta)
  X1 = rbinom(N, 1, p=pi_1)

  Y1 = X2
  pi_3 = plogit(model.matrix(Y1~X1+X2)%*%eta)
  Y1 = rbinom(N, 1, p=pi_3)

  Y2 = Y1
  pi_4 = plogit(model.matrix(Y2~X1+X2+Y1)%*%tau)
  Y2 = rbinom(N, 1, p=pi_4)

  Y3 = Y2
  Y3 = model.matrix(Y3~X1+X2+Y1+Y2)%*%phi+rnorm(N, mean=0, sd=0.5)

  Y4 = Y3
  Y4 = model.matrix(Y4~X1+X2+Y1+Y2+Y3)%*%theta+rnorm(N, mean=0, sd=0.5)

  return(data.frame(X1=X1, X2=X2, Y1=Y1, Y2=Y2, Y3=Y3, Y4=Y4, U=U, W=W, pi=pi))
}

population = gen_pop(dat=micro)
#population = gen_pop(dat=micro, alpha=c(0.06, -0.5, -0.00023792)) # for -0.5 scenario

prob_table = population %>% group_by(X1, X2) %>% summarize(prob_U = mean(U), .groups = 'drop')
population = population %>% left_join(prob_table, by = c("X1", "X2"))

getSample = function(population = population){
  N = dim(population)[1]
  pi = population$pi/10
  sampleID = rbinom(N, 1, pi)
  sub_dat = population[sampleID==1,]

  # unit nonresponse indicator
  sub_dat = sub_dat %>% rowwise() %>% mutate(U = rbinom(1, 1, prob_U)) %>% ungroup()

  # item nonresponse indicator
  Rx1 = sub_dat$X1
  pi_Rx1 = plogit(model.matrix(Rx1~X2+Y1+Y2+Y3, data=sub_dat)%*%c(-1.4, 0.1, 0.1, 0.1, 0.1))
  Rx1 = rbinom(nrow(sub_dat), 1, p=pi_Rx1)

  Ry1 = sub_dat$Y1
  pi_Ry1 = plogit(model.matrix(Ry1~X1+X2+Y2+Y3, data=sub_dat)%*%c(-1.4, 0.1, 0.1, 0.1, 0.1))
  Ry1 = rbinom(nrow(sub_dat), 1, p=pi_Ry1)

  Ry2 = sub_dat$Y2
  pi_Ry2 = plogit(model.matrix(Ry2~X1+X2+Y1+Y3, data=sub_dat)%*%c(-1.4, 0.1, 0.1, 0.1, 0.1))
  Ry2 = rbinom(nrow(sub_dat), 1, p=pi_Ry2)

  Ry4 = sub_dat$Y4
  pi_Ry4 = plogit(model.matrix(Ry4~X1+X2+Y1+Y2+Y3, data=sub_dat)%*%c(-1.4, 0.1, 0.1, 0.1, 0.1, 0.1))
  Ry4 = rbinom(nrow(sub_dat), 1, p=pi_Ry4)

  sub_dat = sub_dat %>% mutate(Rx1=Rx1, Ry1=Ry1, Ry2=Ry2, Ry4=Ry4) %>% select(-prob_U)

  return(sub_dat)
}

# generate 500 sampled datasets
SAMPLE = list()
n_sample = 500
for (i in 1:n_sample){
  SAMPLE[[i]] = getSample(population)
}

population = saveRDS("population.rds")
SAMPLE = saveRDS("SAMPLE.rds")
