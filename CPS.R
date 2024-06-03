library(mvtnorm)
library(survey)
library(mice)
library(tidyverse)

logit = function(x){
  return(log(x)-log(1-x))
}

logit_inv = function(x){
  return(exp(x)/(1+exp(x)))
}

cps18_test_data = read.csv("cps18_more_var.csv")
nc_dat = cps18_test_data[,-1]
summary(nc_dat) 

# adjust weights
nc_adj = nc_dat
total_W = sum(nc_dat$weight, na.rm = T)
nc_adj[which(nc_adj$unit == 0),]$weight = nc_adj[which(nc_adj$unit == 0),]$weight*(1-sum(nc_adj$unit)/nrow(nc_adj))
nc_adj[which(nc_adj$unit == 1),]$weight = total_W/nrow(nc_adj)
nc_adj$pi = 1/nc_adj$weight

# add R
nc_adj$Rv = ifelse(is.na(nc_adj$voted), 1, 0)
nc_adj$Rs = ifelse(is.na(nc_adj$sex), 1, 0)
nc_adj$Ra = ifelse(is.na(nc_adj$age), 1, 0)
nc_adj$Re = ifelse(is.na(nc_adj$race), 1, 0)
nc_adj$Rc = ifelse(is.na(nc_adj$educ), 1, 0)
nc_adj$Rf = ifelse(is.na(nc_adj$faminc), 1, 0)
nc_adj$Rm = ifelse(is.na(nc_adj$marst), 1, 0)
nc_adj$Ro = ifelse(is.na(nc_adj$voteres), 1, 0)
nc_adj$Rp = ifelse(is.na(nc_adj$voteresp), 1, 0)

# reorder columns based on missing rate
columns_ordered = c("diff", "empstat", "sex", "race", "marst", "educ", "age", 
                     "voteresp", "voteres", "voted", "faminc", "weight", "pi", 
                     "unit", "Rs", "Re", "Rm", "Rc", "Ra", "Rp", "Ro", "Rv", "Rf")
nc_adj = nc_adj[, columns_ordered]
nc_adj[, -which(names(nc_adj) %in% c("age", "weight", "pi"))] = lapply(nc_adj[, -which(names(nc_adj) %in% c("age", "weight", "pi"))], factor)
nc_adj[which(nc_adj$unit == 1), c("Rv","Rs","Ra","Re","Rc","Rf","Rm","Ro","Rp")] = NA

# generate 1 completed dataset for item nonresponse only
Imp = mice(nc_adj[1:2013, c(1:11, 15:23)], 
           method = c("", "", "logreg", "polyreg", "polyreg", "polr", "pmm", "logreg", "logreg", "logreg", "polr", "", "", "", "", "", "", "", "", ""),
           m=1, maxit = 5, seed=1) 

datmice.completed = complete(Imp, 1)
datmice.completed$id = rownames(datmice.completed)
datmice.completed$oriweight = cps18_test_data$weight[1:2013]
mydesign = svydesign(id = ~id, data = datmice.completed, weight = ~oriweight)

# Get HT estimators for T_sex (reference: https://github.com/jiuruitang/missing-data/blob/main/4.2_CPS_analysis.R)
T_sex = total_W*0.524
Sd_sex = data.frame(svytotal(~sex, design = mydesign))[2,2]

# HT estimators for race
T_raceW = total_W*0.699
T_raceB = total_W*0.218
T_raceH = total_W*0.039
T_raceR = total_W*0.044
Sd_raceW = data.frame(svytotal(~factor(race), design = mydesign))[1,2]
Sd_raceB = data.frame(svytotal(~factor(race), design = mydesign))[2,2] 
Sd_raceH = data.frame(svytotal(~factor(race), design = mydesign))[3,2] 
Sd_raceR = data.frame(svytotal(~factor(race), design = mydesign))[4,2] 

# HT estimators for vote
T_vote = total_W*0.49
Sd_vote = data.frame(svytotal(~as.factor(voted), design = mydesign))[2,2] 

# use MICE to impute item nonrespondents (generate 50 datasets)
itemimpdat = mice(nc_adj[1:2013, c(1:11, 15:23)], 
                  method=c("", "", "logreg", "polyreg", "polyreg","polr", "pmm", "logreg", "logreg", "logreg", "polr","", "", "", "", "", "", "", "", ""),
                  m=50, maxit = 5, seed=1) # include R
IMP_DAT = list()
n_unit0 = sum(nc_adj$unit == 0) # unit respondents
n_unit1 = sum(nc_adj$unit == 1) # unit nonrespondents

for (i in 1:50) {
  tempdat = complete(itemimpdat, i)
  adj_dat = cbind(rbind(tempdat, tail(nc_adj[, c(1:11, 15:23)], 913)), 
                  nc_adj %>% select(unit, weight, pi))
  
  adj_dat$sex[which(adj_dat$unit == 1)] = rbinom(sum(as.numeric(as.character(adj_dat$unit))),1,0.524)
  adj_dat$race[which(adj_dat$unit == 1)] = sample(c(1:4),sum(as.numeric(as.character(adj_dat$unit))),prob = c(0.699,0.218,0.039,0.044),replace = TRUE)
  adj_dat$voted[which(adj_dat$unit == 1)] = rbinom(sum(as.numeric(as.character(adj_dat$unit))),1,0.49)
  
  #### impute sex ####
  T_sex_hat = rnorm(1, T_sex, Sd_sex)
  total_one_sex = floor( (T_sex_hat-sum((as.numeric(as.character(adj_dat$sex))*adj_dat$weight)[which(adj_dat$unit == 0)])) / 
                           mean(adj_dat$weight[which(adj_dat$unit==1)]) )
  
  m_sex = glm(sex~1, data = adj_dat[which(adj_dat$unit == 0),],family=binomial())
  sex_mat = model.matrix(m_sex)
  alpha_g_origin = alpha_g = rmvnorm(1,mean = coef(m_sex),sigma = vcov(m_sex))
  
  m_sex1 = glm(sex~1,data = adj_dat[which(adj_dat$unit == 1),], family=binomial())
  sex_mat1 = model.matrix(m_sex1)
  
  u_sex = logit(total_one_sex/n_unit1) - mean(sex_mat1%*%t(alpha_g))
  alpha_g[1] = alpha_g[1] + u_sex
  sex_prob = logit_inv(sex_mat1 %*% t(alpha_g))
  adj_dat$sex[which(adj_dat$unit == 1)] = rbinom(n_unit1, 1, prob = sex_prob)
  
  #### Impute race #### 
  T_black_hat = rnorm(1,T_raceB, Sd_raceB)
  T_hist_hat = rnorm(1,T_raceH,Sd_raceH)
  T_rest_hat = rnorm(1,T_raceR,Sd_raceR)
  
  imp_one_rest = max(floor((T_rest_hat - sum(adj_dat$weight[which(adj_dat$unit == 0 & adj_dat$race == 4)]))/mean(adj_dat$weight[which(adj_dat$unit==1)])),0)
  imp_one_black = max(floor((T_black_hat - sum(adj_dat$weight[which(adj_dat$unit == 0 & adj_dat$race == 2)]))/mean(adj_dat$weight[which(adj_dat$unit==1)])),0)
  imp_one_hist = max(floor((T_hist_hat - sum(adj_dat$weight[which(adj_dat$unit == 0 & adj_dat$race == 3)]))/mean(adj_dat$weight[which(adj_dat$unit==1)])),0)
  
  m_eth = multinom(race ~ as.factor(sex), data = adj_dat[which(adj_dat$unit == 0),]) # base level is white(1)
  alpha_e = rmvnorm(1, mean = c(t(coef(m_eth))), sigma = vcov(m_eth))
  eth_mat = model.matrix(m_eth)
  alpha_e_mat = matrix(alpha_e, nrow = 3, byrow = T)
  
  m_eth1 = multinom(race ~ as.factor(sex), 
                    data = adj_dat[which(adj_dat$unit == 1),])
  eth_mat1 = model.matrix(m_eth1)
  
  imp_rest_prob = imp_one_rest/n_unit1
  imp_black_prob = imp_one_black/n_unit1
  imp_hist_prob = imp_one_hist/n_unit1
  imp_white_prob = (n_unit1 - imp_one_rest - imp_one_black - imp_one_hist)/n_unit1
  
  imp_black_prob_overW = imp_black_prob/imp_white_prob
  imp_hist_prob_overW = imp_hist_prob/imp_white_prob
  imp_rest_prob_overW = imp_rest_prob/imp_white_prob
  
  black_adj = logit(imp_black_prob_overW/(1+imp_black_prob_overW)) - mean(eth_mat1 %*%(alpha_e_mat[1,]))
  hist_adj = logit(imp_hist_prob_overW/(1+imp_hist_prob_overW)) - mean(eth_mat1 %*%(alpha_e_mat[2,]))
  rest_adj = logit(imp_rest_prob_overW/(1+imp_rest_prob_overW)) - mean(eth_mat1 %*%(alpha_e_mat[3,]))
  
  # change intercepts
  alpha_e_mat[1] = alpha_e_mat[1] + black_adj
  alpha_e_mat[2] = alpha_e_mat[2] + hist_adj
  alpha_e_mat[3] = alpha_e_mat[3] + rest_adj
  
  black_prob = logit_inv(eth_mat1 %*%(alpha_e_mat[1,]))
  hist_prob = logit_inv(eth_mat1 %*%(alpha_e_mat[2,]))
  rest_prob = logit_inv(eth_mat1 %*%(alpha_e_mat[3,]))
  white_prob = rep(1,n_unit1)
  
  black_prob = black_prob/(1-black_prob)
  hist_prob = hist_prob/(1-hist_prob)
  rest_prob = rest_prob/(1-rest_prob)
  
  black_prob_adj = black_prob/(black_prob+hist_prob+rest_prob+white_prob)
  hist_prob_adj = hist_prob/(black_prob+hist_prob+rest_prob+white_prob)
  rest_prob_adj = rest_prob/(black_prob+hist_prob+rest_prob+white_prob)
  white_prob_adj = white_prob/(black_prob+hist_prob+rest_prob+white_prob)
  
  race_prob = matrix(NA, n_unit1, 4)
  race_prob[,1] = white_prob_adj
  race_prob[,2] = black_prob_adj
  race_prob[,3] = hist_prob_adj
  race_prob[,4] = rest_prob_adj
  
  adj_dat$race[which(adj_dat$unit == 1)] = sapply(1:n_unit1, function(x)sample(c(1,2,3,4),size = 1, replace = FALSE, prob = race_prob[x,]))
  
  #### Impute vote
  T_vote_hat = rnorm(1, T_vote, Sd_vote)
  total_one_vote = floor((T_vote_hat-sum((as.numeric(as.character(adj_dat$voted))*adj_dat$weight)[which(adj_dat$unit == 0)]))/mean(adj_dat$weight[which(adj_dat$unit==1)]))
  
  m_vote = glm(voted ~ as.factor(sex) + as.factor(race) + as.factor(sex)*as.factor(race),
               data = adj_dat[which(adj_dat$unit == 0),],family=binomial())
  vote_mat = model.matrix(m_vote)
  alpha_v_origin = alpha_v = rmvnorm(1, mean = c(t(coef(m_vote))), sigma = vcov(m_vote))
  
  vote_mat1 = vote_mat[1:n_unit1, ]
  vote_mat1[,1] = rep(1,n_unit1) # intercept column
  vote_mat1[,2] = adj_dat$sex[which(adj_dat$unit == 1)] # sex column
  vote_mat1[,3] = ifelse(adj_dat$race[which(adj_dat$unit == 1)] == 2, 1, 0) # race columns
  vote_mat1[,4] = ifelse(adj_dat$race[which(adj_dat$unit == 1)] == 3, 1, 0)
  vote_mat1[,5] = ifelse(adj_dat$race[which(adj_dat$unit == 1)] == 4,1,0)
  vote_mat1[,6] = ifelse(adj_dat$sex[which(adj_dat$unit == 1)] == 1 &
                           adj_dat$race[which(adj_dat$unit == 1)] == 2,1,0) # interactions: sex*race
  vote_mat1[,7] = ifelse(adj_dat$sex[which(adj_dat$unit == 1)] == 1 &
                           adj_dat$race[which(adj_dat$unit == 1)] == 3,1,0)
  vote_mat1[,8] = ifelse(adj_dat$sex[which(adj_dat$unit == 1)] == 1 &
                           adj_dat$race[which(adj_dat$unit == 1)] == 4,1,0)
  # change intercepts
  u_adj = logit(total_one_vote/n_unit1) - mean(vote_mat1%*%t(alpha_v))
  alpha_v[1] = alpha_v[1] + u_adj
  vote_prob = logit_inv(vote_mat1%*%t(alpha_v))
  adj_dat$voted[which(adj_dat$unit == 1)] = rbinom(n_unit1, 1, prob = vote_prob)
  
  #### hot deck ####
  cell_010 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==1 & adj_dat$voted==0)
  cell_011 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==1 & adj_dat$voted==1)
  cell_020 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==2 & adj_dat$voted==0)
  cell_021 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==2 & adj_dat$voted==1)
  cell_030 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==3 & adj_dat$voted==0)
  cell_031 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==3 & adj_dat$voted==1)
  cell_040 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==4 & adj_dat$voted==0)
  cell_041 = which(adj_dat$unit==0 & adj_dat$sex==0 & adj_dat$race==4 & adj_dat$voted==1)
  cell_110 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==1 & adj_dat$voted==0)
  cell_111 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==1 & adj_dat$voted==1)
  cell_120 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==2 & adj_dat$voted==0)
  cell_121 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==2 & adj_dat$voted==1)
  cell_130 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==3 & adj_dat$voted==0)
  cell_131 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==3 & adj_dat$voted==1)
  cell_140 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==4 & adj_dat$voted==0)
  cell_141 = which(adj_dat$unit==0 & adj_dat$sex==1 & adj_dat$race==4 & adj_dat$voted==1)
  
  imp_lookup = list("010" = cell_010, "011" = cell_011, "020" = cell_020, "021" = cell_021,
                    "030" = cell_030, "031" = cell_031, "040" = cell_040, "041" = cell_041,
                    "110" = cell_110, "111" = cell_111, "120" = cell_120, "121" = cell_121,
                    "130" = cell_130, "131" = cell_131, "140" = cell_140, "141" = cell_141)

  rows_to_impute = which(adj_dat$unit==1)
  for (k in rows_to_impute) {
    imp_key = paste0(adj_dat[k, "sex"], adj_dat[k, "race"], adj_dat[k, "voted"])
    adj_dat[k, c("diff", "empstat", "marst", "educ", "age", "voteresp", "voteres", "faminc")] = adj_dat[sample(imp_lookup[[imp_key]], 1), 
                                                                                                        c("diff", "empstat", "marst", "educ", "age", "voteresp", "voteres", "faminc")]
  }
  IMP_DAT[[i]] = adj_dat
}

#### marginal distributions ####
ans_S = matrix(NA,50,2)
ul_S = matrix(NA,50,2)
ans_R = matrix(NA,50,4)
ul_R = matrix(NA,50,4)

for (j in 1:50){
  IMP_dat = cbind(id=as.numeric(rownames(nc_adj)), IMP_DAT[[j]][,c(1:11, 22:23)])
  mydesign = svydesign(id = ~id, data = IMP_dat, weight = ~weight)
  
  wgt_df_R = data.frame(svymean(~as.factor(race), design=mydesign))
  ans_R[j,] = wgt_df_R[,1]
  ul_R[j,] = wgt_df_R[,2]
  
  wgt_df_S = data.frame(svymean(~as.factor(sex), design=mydesign))
  ans_S[j,] = wgt_df_S[,1]
  ul_S[j,] = wgt_df_S[,2]
}

colMeans(ans_S)
u_L_S = colMeans(ul_S^2)
b_L_S = apply((scale(ans_S,scale=FALSE))^2/(50-1),2,sum)
T_L_S = (1+1/50)*b_L_S + u_L_S
sqrt(T_L_S) 

colMeans(ans_R) 
u_L_R = colMeans(ul_R^2)
b_L_R = apply((scale(ans_R,scale=FALSE))^2/(50-1),2,sum)
T_L_R = (1+1/50)*b_L_R + u_L_R
sqrt(T_L_R) 

#### marginal comparison: MICE ####
ans_S_mice = matrix(NA,50,2)
ul_S_mice = matrix(NA,50,2)
ans_R_mice = matrix(NA,50,4)
ul_R_mice = matrix(NA,50,4)

itemimpdat = mice(nc_adj[1:2013, c(1:11, 15:23)], 
                  method=c("", "", "logreg", "polyreg", "polyreg","polr", "pmm", "logreg", "logreg", "logreg", "polr","", "", "", "", "", "", "", "", ""),
                  m=50, maxit = 5, seed=1)

completemice = list()
for (i in 1:50) {
  completemice[[i]] = complete(itemimpdat, i)
  completemice[[i]]$id = rownames(nc_dat[nc_dat$unit==0,])
  completemice[[i]]$oriweight = cps18_test_data[cps18_test_data$unit==0,]$weight
  mydesign = svydesign(id = ~id, data = completemice[[i]], weight = ~oriweight)
  
  wgt_df_R = data.frame(svymean(~as.factor(race), design=mydesign))
  ans_R_mice[i,] = wgt_df_R[,1]
  ul_R_mice[i,] = wgt_df_R[,2]
  
  wgt_df_S = data.frame(svymean(~as.factor(sex), design=mydesign))
  ans_S_mice[i,] = wgt_df_S[,1]
  ul_S_mice[i,] = wgt_df_S[,2]
}

colMeans(ans_S_mice)
u_L_S_mice = colMeans(ul_S_mice^2)
b_L_S_mice = apply((scale(ans_S_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_S_mice = (1+1/50)*b_L_S_mice + u_L_S_mice
sqrt(T_L_S_mice) 

colMeans(ans_R_mice)
u_L_R_mice = colMeans(ul_R_mice^2)
b_L_R_mice = apply((scale(ans_R_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_R_mice = (1+1/50)*b_L_R_mice + u_L_R_mice 
sqrt(T_L_R_mice)

#### voter turnout across subgroups ####
ans_S = matrix(NA,50,2)
ul_S = matrix(NA,50,2)
ans_R = matrix(NA,50,4)
ul_R = matrix(NA,50,4)
ans_C = matrix(NA,50,3)
ul_C = matrix(NA,50,3)
ans = ul = matrix(NA,50,2)
ans_F = ul_F = matrix(NA,50,3)
ans_T = ul_T = matrix(NA,50,3)
ans_M = ul_M = matrix(NA,50,3)
ans_D = ul_D = matrix(NA,50,2)
ans_O = ul_O = matrix(NA,50,2)
ans_P = ul_P = matrix(NA,50,2)

for (j in 1:50){
  IMP_dat = cbind(id=as.numeric(rownames(nc_adj)), IMP_DAT[[j]][,c(1:11, 22:23)])
  mydesign = svydesign(id = ~id, data = IMP_dat, weight = ~weight)
  
  wgt_tblR = svyby(~as.factor(voted), ~as.factor(race), mydesign, svymean)
  wgt_dfR = data.frame(wgt_tblR)
  ans_R[j,] = wgt_dfR[,3]
  ul_R[j,] = wgt_dfR[,5]
  
  wgt_tblS = svyby(~as.factor(voted),~as.factor(sex),mydesign,svymean)
  wgt_dfS = data.frame(wgt_tblS)
  ans_S[j,] = wgt_dfS[,3]
  ul_S[j,] = wgt_dfS[,5]
  
  wgt_tblC = svyby(~as.factor(voted),~as.factor(educ),mydesign,svymean)
  wgt_dfC = data.frame(wgt_tblC)
  ans_C[j,] = wgt_dfC[,3]
  ul_C[j,] = wgt_dfC[,5]
  
  wgt_tblF = svyby(~as.factor(voted),~as.factor(faminc),mydesign,svymean)
  wgt_dfF = data.frame(wgt_tblF)
  ans_F[j,] = wgt_dfF[,3]
  ul_F[j,] = wgt_dfF[,5]
  
  wgt_tblT = svyby(~as.factor(voted),~as.factor(empstat),mydesign,svymean)
  wgt_dfT = data.frame(wgt_tblT)
  ans_T[j,] = wgt_dfT[,3]
  ul_T[j,] = wgt_dfT[,5]
  
  wgt_tblM = svyby(~as.factor(voted),~as.factor(marst),mydesign,svymean)
  wgt_dfM = data.frame(wgt_tblM)
  ans_M[j,] = wgt_dfM[,3]
  ul_M[j,] = wgt_dfM[,5]
  
  wgt_tblD = svyby(~as.factor(voted),~as.factor(diff),mydesign,svymean)
  wgt_dfD = data.frame(wgt_tblD)
  ans_D[j,] = wgt_dfD[,3]
  ul_D[j,] = wgt_dfD[,5]
  
  wgt_tblO = svyby(~as.factor(voted),~as.factor(voteres),mydesign,svymean)
  wgt_dfO = data.frame(wgt_tblO)
  ans_O[j,] = wgt_dfO[,3]
  ul_O[j,] = wgt_dfO[,5]
  
  wgt_tblP = svyby(~as.factor(voted),~as.factor(voteresp),mydesign,svymean)
  wgt_dfP = data.frame(wgt_tblP)
  ans_P[j,] = wgt_dfP[,3]
  ul_P[j,] = wgt_dfP[,5]
  
  wgt_df = data.frame(svymean(~voted, design=mydesign))
  ans[j,] = wgt_df[,1]
  ul[j,] = wgt_df[,2]
}

colMeans(ans) 
u_L = mean(ul^2)
b_L = apply((scale(ans,scale=FALSE))^2/(50-1),2,sum)
T_L = (1+1/50)*b_L + u_L

colMeans(ans_S) 
u_L_S = colMeans(ul_S^2)
b_L_S = apply((scale(ans_S,scale=FALSE))^2/(50-1),2,sum)
T_L_S = (1+1/50)*b_L_S + u_L_S

colMeans(ans_R) 
u_L_R = colMeans(ul_R^2)
b_L_R = apply((scale(ans_R,scale=FALSE))^2/(50-1),2,sum)
T_L_R = (1+1/50)*b_L_R + u_L_R 

colMeans(ans_C) 
u_L_C = colMeans(ul_C^2)
b_L_C = apply((scale(ans_C,scale=FALSE))^2/(50-1),2,sum)
T_L_C = (1+1/50)*b_L_C + u_L_C 

colMeans(ans_F) 
u_L_F = colMeans(ul_F^2)
b_L_F = apply((scale(ans_F,scale=FALSE))^2/(50-1),2,sum)
T_L_F = (1+1/50)*b_L_F + u_L_F 

colMeans(ans_T) 
u_L_T = colMeans(ul_T^2)
b_L_T = apply((scale(ans_T,scale=FALSE))^2/(50-1),2,sum)
T_L_T = (1+1/50)*b_L_T + u_L_T

colMeans(ans_M) 
u_L_M = colMeans(ul_M^2)
b_L_M = apply((scale(ans_M,scale=FALSE))^2/(50-1),2,sum)
T_L_M = (1+1/50)*b_L_M + u_L_M
T_L_M

colMeans(ans_D) 
u_L_D = colMeans(ul_D^2)
b_L_D = apply((scale(ans_D,scale=FALSE))^2/(50-1),2,sum)
T_L_D = (1+1/50)*b_L_D + u_L_D
T_L_D

colMeans(ans_O) 
u_L_O = colMeans(ul_O^2)
b_L_O = apply((scale(ans_O,scale=FALSE))^2/(50-1),2,sum)
T_L_O = (1+1/50)*b_L_O + u_L_O
T_L_O

colMeans(ans_P)
u_L_P = colMeans(ul_P^2)
b_L_P = apply((scale(ans_P,scale=FALSE))^2/(50-1),2,sum)
T_L_P = (1+1/50)*b_L_P + u_L_P
T_L_P

#### comparison: MICE ####
ans_S_mice = matrix(NA,50,2)
ul_S_mice = matrix(NA,50,2)
ans_R_mice = matrix(NA,50,4)
ul_R_mice = matrix(NA,50,4)
ans_C_mice = matrix(NA,50,3)
ul_C_mice = matrix(NA,50,3)
ans_mice = ul_mice = matrix(NA,50, 2)
ans_F_mice = ul_F_mice = matrix(NA,50,3)
ans_T_mice = ul_T_mice = matrix(NA,50,3)
ans_M_mice = ul_M_mice = matrix(NA,50,3)
ans_D_mice = ul_D_mice = matrix(NA,50,2)
ans_O_mice = ul_O_mice = matrix(NA,50,2)
ans_P_mice = ul_P_mice = matrix(NA,50,2)

completemice = list()
for (i in 1:50) {
  completemice[[i]] = complete(itemimpdat, i)
  completemice[[i]]$id = rownames(nc_dat[nc_dat$unit==0,])
  completemice[[i]]$oriweight = cps18_test_data[cps18_test_data$unit==0,]$weight
  mydesign = svydesign(id = ~id, data = completemice[[i]], weight = ~oriweight)
  
  wgt_tblR = svyby(~as.factor(voted), ~as.factor(race), mydesign, svymean)
  wgt_dfR = data.frame(wgt_tblR)
  ans_R_mice[i,] = wgt_dfR[,3]
  ul_R_mice[i,] = wgt_dfR[,5]
  
  wgt_tblS = svyby(~as.factor(voted),~as.factor(sex),mydesign,svymean)
  wgt_dfS = data.frame(wgt_tblS)
  ans_S_mice[i,] = wgt_dfS[,3]
  ul_S_mice[i,] = wgt_dfS[,5]
  
  wgt_tblC = svyby(~as.factor(voted),~as.factor(educ),mydesign,svymean)
  wgt_dfC = data.frame(wgt_tblC)
  ans_C_mice[i,] = wgt_dfC[,3]
  ul_C_mice[i,] = wgt_dfC[,5]
  
  wgt_tblF = svyby(~as.factor(voted),~as.factor(faminc),mydesign,svymean)
  wgt_dfF = data.frame(wgt_tblF)
  ans_F_mice[i,] = wgt_dfF[,3]
  ul_F_mice[i,] = wgt_dfF[,5]
  
  wgt_tblT = svyby(~as.factor(voted),~as.factor(empstat),mydesign,svymean)
  wgt_dfT = data.frame(wgt_tblT)
  ans_T_mice[i,] = wgt_dfT[,3]
  ul_T_mice[i,] = wgt_dfT[,5]
  
  wgt_tblM = svyby(~as.factor(voted),~as.factor(marst),mydesign,svymean)
  wgt_dfM = data.frame(wgt_tblM)
  ans_M_mice[i,] = wgt_dfM[,3]
  ul_M_mice[i,] = wgt_dfM[,5]
  
  wgt_tblD = svyby(~as.factor(voted),~as.factor(diff),mydesign,svymean)
  wgt_dfD = data.frame(wgt_tblD)
  ans_D_mice[i,] = wgt_dfD[,3]
  ul_D_mice[i,] = wgt_dfD[,5]
  
  wgt_tblO = svyby(~as.factor(voted),~as.factor(voteres),mydesign,svymean)
  wgt_dfO = data.frame(wgt_tblO)
  ans_O_mice[i,] = wgt_dfO[,3]
  ul_O_mice[i,] = wgt_dfO[,5]
  
  wgt_tblP = svyby(~as.factor(voted),~as.factor(voteresp),mydesign,svymean)
  wgt_dfP = data.frame(wgt_tblP)
  ans_P_mice[i,] = wgt_dfP[,3]
  ul_P_mice[i,] = wgt_dfP[,5]
  
  wgt_df = data.frame(svymean(~voted, design=mydesign))
  ans_mice[i,] = wgt_df[,1]
  ul_mice[i,] = wgt_df[,2]
}

colMeans(ans_mice) 
u_L_mice = colMeans(ul_mice^2)
b_L_mice = apply((scale(ans_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_mice = (1+1/50)*b_L_mice + u_L_mice

colMeans(ans_S_mice) 
u_L_S_mice = colMeans(ul_S_mice^2)
b_L_S_mice = apply((scale(ans_S_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_S_mice = (1+1/50)*b_L_S_mice + u_L_S_mice

colMeans(ans_R_mice) 
u_L_R_mice = colMeans(ul_R_mice^2)
b_L_R_mice = apply((scale(ans_R_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_R_mice = (1+1/50)*b_L_R_mice + u_L_R_mice 

colMeans(ans_C_mice) 
u_L_C_mice = colMeans(ul_C_mice^2)
b_L_C_mice = apply((scale(ans_C_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_C_mice = (1+1/50)*b_L_C_mice + u_L_C_mice 

colMeans(ans_F_mice) 
u_L_F_mice = colMeans(ul_F_mice^2)
b_L_F_mice = apply((scale(ans_F_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_F_mice = (1+1/50)*b_L_F_mice + u_L_F_mice 

colMeans(ans_T_mice)
u_L_T_mice = colMeans(ul_T_mice^2)
b_L_T_mice = apply((scale(ans_T_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_T_mice = (1+1/50)*b_L_T_mice + u_L_T_mice 

colMeans(ans_M_mice)
u_L_M_mice = colMeans(ul_M_mice^2)
b_L_M_mice = apply((scale(ans_M_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_M_mice = (1+1/50)*b_L_M_mice + u_L_M_mice 

colMeans(ans_D_mice)
u_L_D_mice = colMeans(ul_D_mice^2)
b_L_D_mice = apply((scale(ans_D_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_D_mice = (1+1/50)*b_L_D_mice + u_L_D_mice 

colMeans(ans_O_mice)
u_L_O_mice = colMeans(ul_O_mice^2)
b_L_O_mice = apply((scale(ans_O_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_O_mice = (1+1/50)*b_L_O_mice + u_L_O_mice 

colMeans(ans_P_mice)
u_L_P_mice = colMeans(ul_P_mice^2)
b_L_P_mice = apply((scale(ans_P_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_P_mice = (1+1/50)*b_L_P_mice + u_L_P_mice 

#### voter turnout by educ and diff #### 
ans_CD = ans_CD_mice = matrix(NA,50,6)
ul_CD = ul_CD_mice = matrix(NA,50,6)

for (j in 1:50){
  IMP_dat = cbind(id=as.numeric(rownames(nc_adj)), IMP_DAT[[j]][,c(1:11, 22:23)])
  mydesign = svydesign(id = ~id, data = IMP_dat, weight = ~weight)
  wgt_df_CD = data.frame(svyby(~as.factor(voted),~as.factor(diff)+~as.factor(educ),mydesign,svymean))
  ans_CD[j,] = wgt_df_CD[,4]
  ul_CD[j,] = wgt_df_CD[,6] 
}

colMeans(ans_CD) 
u_L_CD = colMeans(ul_CD^2)
b_L_CD = apply((scale(ans_CD,scale=FALSE))^2/(50-1),2,sum)
T_L_CD = (1+1/50)*b_L_CD + u_L_CD 

completemice = list()
for (i in 1:50) {
  completemice[[i]] = complete(itemimpdat, i)
  completemice[[i]]$id = rownames(nc_dat[nc_dat$unit==0,])
  completemice[[i]]$oriweight = cps18_test_data[cps18_test_data$unit==0,]$weight
  mydesign = svydesign(id = ~id, data = completemice[[i]], weight = ~oriweight)
  
  wgt_df_CD = data.frame(svyby(~as.factor(voted),~as.factor(diff)+~as.factor(educ), mydesign, svymean))
  ans_CD_mice[i,] = wgt_df_CD[,4]
  ul_CD_mice[i,] = wgt_df_CD[,6]
}

colMeans(ans_CD_mice)
u_L_CD_mice = colMeans(ul_CD_mice^2)
b_L_CD_mice = apply((scale(ans_CD_mice,scale=FALSE))^2/(50-1),2,sum)
T_L_CD_mice = (1+1/50)*b_L_CD_mice + u_L_CD_mice