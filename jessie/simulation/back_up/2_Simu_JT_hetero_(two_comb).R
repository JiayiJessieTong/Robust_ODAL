library(mvtnorm)
library(ggplot2)
library(matlib)

##################################################
###Distributed algorithm for logistic regression##
### homogenrous
### setting one (K = 10, each n = seq(100, 1000, by=100), total N = 1000), 
##################################################
setwd("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation")
# column names of result data.frame
col_names_list = c("iter",
                   "odal1 mean beta0", "odal1 mean beta1", "odal1 mean beta2", 
                   "odal1 median beta0", "odal1 median beta1", "odal1 median beta2", 
                   "odal2 beta0 mean mean", "odal2 beta1 mean mean", "odal2 beta2 mean mean", 
                   "odal2 beta0 mean median", "odal2 beta1 mean median", "odal2 beta2 mean median",  
                   "odal2 beta0 median mean", "odal2 beta1 median mean", "odal2 beta2 median mean", 
                   "odal2 beta0 med med", "odal2 beta1 med med", "odal2 beta2 med med",
                   "odal1 mean var0", "odal1 mean var1", "odal1 mean var2", 
                   "odal1 median var0", "odal1 median var1", "odal1 median var2",
                   "odal2 var0 mean mean", "odal2 var1 mean mean", "odal2 var2 mean mean", 
                   "odal2 var0 mean med", "odal2 var1 mean med", "odal2 var2 mean med", 
                   "odal2 var0 med mean", "odal2 var1 med mean", "odal2 var2 med mean",
                   "odal2 var0 med med", "odal2 var1 med med", "odal2 var2 med med",
                   "local beta0", "local beta1", "local beta2",
                   "pooled beta0", "pooled beta1", "pooled beta2", 
                   "local var0", "local var1", "local var2", 
                   "pooled var0", "pooled var1", "pooled var2")


##########functions############################
expit = function(x){exp(x)/(1+exp(x))}

#likelihood function for logistic regression, the input X is a n*p matrix where 
#each patient has p covariates stored in each row.
Lik = function(beta,X,Y){
  design = cbind(1,X)
  sum(Y*(design%*%t(t(beta)))-log(1+exp(design%*%t(t(beta)))))/length(Y)
}


#Y is the binary vector with disease status per patient (1 indicates case and 0 indicates control)

#first order gradient
#beta:p*1 X:n*p Y:n*1
# mean
Lgradient = function(beta,X,Y){
  design = cbind(1,X)
  t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
}

# median
Lgradient_median = function(beta, X, Y, k){
  l = seq(1:length(1:dim(X)[1]))
  # ind = sample(rep(1:k, each = length(l)/k))
  # sp_ind = split(l,ind)
  
  num = length(l)/k
  mat = list()
  for (i in 1:k){
    # id = sp_ind[[i]]
    id = seq(1:num)+(i-1)*num
    X_single = X[id,]
    design = cbind(1,X_single)
    mat[[i]] =  t(Y[id]-expit(design%*%t(t(beta))))%*%design/length(Y[id])
  }
  final_mat = apply(simplify2array(mat), 2, median)
  final_mat
}

# second order gradient
# mean
Lgradient2 =function(beta,X){
  design = cbind(1,X)
  Z=expit(design%*%beta)
  t(c(-1*Z*(1-Z))*design)%*%design/nrow(X)
}

# k is the number of hospitals
Lgradient2_median =function(beta,X,k){
  l = seq(1:length(1:dim(X)[1]))
  # ind = sample(rep(1:k, each = length(l)/k))
  # sp_ind = split(l,ind)

  num = length(l)/k
  mat = list()
  for (i in 1:k){
    # id = sp_ind[[i]]
    id = seq(1:num)+(i-1)*num
    X_single = X[id,]
    design = cbind(1,X_single)
    Z=expit(design%*%beta)
    mat[[i]] = t(c(-1*Z*(1-Z))*design)%*%design/nrow(X_single)
  }
  final_mat = apply(simplify2array(mat), 1:2, median)
  final_mat
}

# gradient for variance
Lgradient_meat = function(beta,X,Y){
  design = cbind(1,X)
  Z = expit(design%*%beta)
  t(c(Y-Z)*design)%*%(c(Y-Z)*design)/nrow(X)
}

# sandwich variance
Sandwich = function(beta,X,Y,N){  
  # vec_L1 = Lgradient(beta,X,Y)
  # mat_L1 = t(length(Y)/N *vec_L1)%*%t(t(vec_L1))
  mat_L1 = Lgradient_meat(beta,X,Y)
  #mat_L2 = diag(1,5)
  mat_L2 = Lgradient2(beta,X)
  
  inv_L2 = solve.default(mat_L2)

  out = inv_L2%*%mat_L1%*%inv_L2/N

  return(out)
}

#generate the Data to Simulation, nn is the size of each hospital (N/k)
Generate = function(N,beta,nn,x1_mu,x1_var,x2_var,b0,b1){
  
  X1 = rnorm(N-nn)         #generate N numbers with N(0,1)
  X2 = rbinom(N-nn,1,0.3)  #generate N numbers with B(1,0.3)
  
  X = cbind(1,X1,X2)
  meanY = expit(X%*%beta)
  Y = rbinom(N-nn,1,meanY)
  
  allMinisHetero = cbind(X1,X2,Y)
  
  # the heterogeneous site
  X1.hetero = rnorm(nn,x1_mu,x1_var)
  X2.hetero = rbinom(nn,1,x2_var)
  # X1.hetero = rnorm(2*nn,1,1)
  # X2.hetero = rbinom(2*nn,1,0.5)
   
  beta.hetero = c(b0, b1, -1)
  # beta.hetero = c(1.5, 0.8, -1)
  X.hetero = cbind(1,X1.hetero,X2.hetero)
  meanY.hetero = expit(X.hetero%*%beta.hetero)
  Y.hetero = rbinom(nn,1,meanY.hetero)
  
  Hetero = cbind(X1.hetero, X2.hetero, Y.hetero)

  
  # return(rbind(allMinisHetero))
  return(rbind(allMinisHetero, Hetero))
}

#simulation
simu = function(beta,N,p,nn,x1_mu,x1_var,x2_var,b0,b1){
  # ODAL 1: mean of gradient of each site (only homo site)
  SL_homo = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L_homo%*%beta
  }
  
  # ODAL 1: median of gradient of each site (only homo site)
  SL_median_homo = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L_median_homo%*%beta
  }
  
  # ODAL 1: mean of gradient of each site
  SL = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L%*%beta
  }
  
  # ODAL 1: median of gradient of each site
  SL_median = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L_median%*%beta
  }
  
  # ODAL 2: mean of first gradient, mean of second gradient
  SL2_mean_mean = function(beta){
    beta_temp = beta - beta0
    #mat_temp = t(t(beta_temp))%*%t(beta_temp)
    -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2%*%t(t(beta_temp))
  }
  
  # ODAL 2: mean of first gradient, median of second gradient
  SL2_mean_median = function(beta){
    beta_temp = beta - beta0
    #mat_temp = t(t(beta_temp))%*%t(beta_temp)
    -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2_median%*%t(t(beta_temp))
  }
  
  # ODAL 2: median of first gradient, mean of second gradient
  SL2_median_mean = function(beta){
    beta_temp = beta - beta0
    #mat_temp = t(t(beta_temp))%*%t(beta_temp)
    -Lik(beta,Xlocal,Ylocal) - L_median%*%beta - 0.5*t(beta_temp)%*%L2%*%t(t(beta_temp))
  }
  
  # ODAL 2: median of first gradient, median of second gradient
  SL2_median_median = function(beta){
    beta_temp = beta - beta0
    #mat_temp = t(t(beta_temp))%*%t(beta_temp)
    -Lik(beta,Xlocal,Ylocal) - L_median%*%beta - 0.5*t(beta_temp)%*%L2_median%*%t(t(beta_temp))
  }
  
  
  #generate the data
  Data = Generate(N,beta,nn,x1_mu,x1_var,x2_var,b0,b1)
  # nlocal = round(p*N,1)   #round to 1e-1
  X = Data[,1:2]
  Y = Data[,3]
  Xlocal = Data[1:nn,1:2]
  Ylocal = Data[1:nn,3]
  X.homo = Data[1:(N-nn),1:2]
  Y.homo = Data[1:(N-nn),3]
  
  #Local beta and var
  E_local = summary(glm(Ylocal~Xlocal, family = "binomial"(link = "logit")))
  beta0 = E_local$coefficients[,1]
  std0 = E_local$coefficients[,2]
  var0 = std0^2
  
  #pooled beta and var
  fit = summary(glm(Y~X, family = "binomial"(link = "logit")))
  beta_all = fit$coefficients[,1]
  std_all = fit$coefficients[,2]
  var_all =std_all^2
  
  # meta method 
  l = seq(1:dim(X)[1])
  num = length(l)*p
  beta_meta_list = c()
  se_meta_list = c()
  for (i in c(1:(1/p))){
    id = seq(1:num)+(i-1)*num
    X_each = X[id,]
    Y_each = Y[id]
    fit_each = summary(glm(Y_each~X_each, family = "binomial"(link = "logit")))
    beta_meta_list[i] = fit_each$coefficients[3,1]
    se_meta_list[i] = fit_each$coefficients[3,2]
  }
  beta_meta <-metagen(beta_meta_list, se_meta_list,comb.fixed = TRUE,comb.random = FALSE,prediction=TRUE,sm="OR")$TE.fixed
  
  iter = fit$iter+1
  
  #estimatior beta and var
  L = Lgradient(beta0,X,Y)
  L_median = Lgradient_median(beta0,X,Y,1/p)
  
  L_homo = Lgradient(beta0,X.homo,Y.homo)
  L_median_homo = Lgradient_median(beta0,X.homo,Y.homo,1/p)
  
  L2 = Lgradient2(beta0,X) - Lgradient2(beta0,Xlocal)
  L2_median = Lgradient2_median(beta0,X,1/p) - Lgradient2(beta0,Xlocal)
  
  #mean of the estimator
  beta_tilde1_mean_homo = optim(beta0,SL_homo,control = list(maxit = 10000,reltol = 1e-10))$par
  
  beta_tilde1_median_homo = optim(beta0,SL_median_homo,control = list(maxit = 10000,reltol = 1e-10))$par
  
  beta_tilde1_mean = optim(beta0,SL,control = list(maxit = 10000,reltol = 1e-10))$par
  
  beta_tilde1_median = optim(beta0,SL_median,control = list(maxit = 10000,reltol = 1e-10))$par
  
  beta_tilde2_mean_mean = optim(beta_tilde1_mean,SL2_mean_mean,control = list(maxit = 10000,reltol = 1e-10))$par

  beta_tilde2_mean_median = optim(beta_tilde1_mean,SL2_mean_median,control = list(maxit = 10000,reltol = 1e-10))$par

  beta_tilde2_median_mean = optim(beta_tilde1_mean,SL2_median_mean,control = list(maxit = 10000,reltol = 1e-10))$par

  beta_tilde2_median_median = optim(beta_tilde1_mean,SL2_median_median,control = list(maxit = 10000,reltol = 1e-10))$par

  #var of the estimator
  # var_tilde1_mean_mat = Sandwich(beta_tilde1_mean,Xlocal,Ylocal,N)
  # var_tilde1_median_mat = Sandwich(beta_tilde1_median,Xlocal,Ylocal,N)
  # 
  # var_tilde2_mean_mean_mat = Sandwich(beta_tilde2_mean_mean,Xlocal,Ylocal,N)
  # var_tilde2_mean_median_mat = Sandwich(beta_tilde2_mean_median,Xlocal,Ylocal,N)
  # var_tilde2_median_mean_mat = Sandwich(beta_tilde2_median_mean,Xlocal,Ylocal,N)
  # var_tilde2_median_median_mat = Sandwich(beta_tilde2_median_median,Xlocal,Ylocal,N)
  # 
  # var_tilde1_mean = diag(var_tilde1_mean_mat)
  # var_tilde1_median = diag(var_tilde1_median_mat)
  # 
  # var_tilde2_mean_mean = diag(var_tilde2_mean_mean_mat)
  # var_tilde2_mean_median = diag(var_tilde2_mean_median_mat)
  # var_tilde2_median_mean = diag(var_tilde2_median_mean_mat)
  # var_tilde2_median_median = diag(var_tilde2_median_median_mat)
  
  #output the outcome
  total = c(iter,
            beta_tilde1_mean_homo, beta_tilde1_median_homo,
            beta_tilde1_mean, beta_tilde1_median,
            beta_tilde2_mean_mean, beta_tilde2_mean_median, beta_tilde2_median_mean, beta_tilde2_median_median,
            beta0,beta_all,beta_meta)
  # total = c(iter, beta_tilde1_mean, beta_tilde1_median, beta_tilde2_mean_mean, beta_tilde2_mean_median, beta_tilde2_median_mean, beta_tilde2_median_median,
  #           var_tilde1_mean, var_tilde1_median, var_tilde2_mean_mean, var_tilde2_mean_median, var_tilde2_median_mean, var_tilde2_median_median, 
  #           beta0,beta_all,var0,var_all)
  # total = c(iter,beta0,beta_all,var0,var_all,
  #           beta_tilde1,beta_tilde2,beta_tilde2_median,
  #           var_tilde1,var_tilde2,var_tilde2_median)
  return(total)
}



####Simulation Setting######
###Overall Setting##########
# True beta#
set.seed(1111)
beta_true = c(1,1,-1)
# setting 1
K1 = 10
nn1 = seq(100,1000,by=100)   #number in a site
NN1 = nn1*K1
p1 = 0.1
Nsim = 50
Result_setting1 = matrix(0,nrow =length(NN1) ,ncol = 32)
for(i in 1:length(NN1)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN1[i],p1,
                            nn = nn1[i],x1_mu = 1,x1_var = 1,x2_var = 0.5,b0= 1.5,b1 = 0.8))  #run func with Nsim times
  set_out1 = apply(out,1,mean)
  Result_setting1[i,] = set_out1
}

# setting 2
set.seed(1111)
Result_setting2 = matrix(0,nrow =length(NN1) ,ncol = 32)
for(i in 1:length(NN1)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN1[i],p1,
                            nn = 2*nn1[i],x1_mu = 1,x1_var = 1,x2_var = 0.5,beta0= 1.5,beta1 = 0.8))  #run func with Nsim times
  set_out2 = apply(out,1,mean)
  Result_setting2[i,] = set_out2
}

# setting 3
set.seed(1111)
Result_setting3 = matrix(0,nrow =length(NN1) ,ncol = 32)
for(i in 1:length(NN1)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN1[i],p1,
                            nn = nn1[i],x1_mu = 5,x1_var = 5,x2_var = 0.9,beta0= 3,beta1 = 2))  #run func with Nsim times
  set_out3 = apply(out,1,mean)
  Result_setting3[i,] = set_out3
}

# setting 4
set.seed(1111)
Result_setting4 = matrix(0,nrow =length(NN1) ,ncol = 32)
for(i in 1:length(NN1)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN1[i],p1,
                            nn = 2*nn1[i],x1_mu = 5,x1_var = 5,x2_var = 0.9,beta0= 3,beta1 = 2))  #run func with Nsim times
  set_out4 = apply(out,1,mean)
  Result_setting4[i,] = set_out4
}


write.csv(Result_setting1, "two_var/Result_setting1.csv")
write.csv(Result_setting2, "two_var/Result_setting2.csv")
write.csv(Result_setting3, "two_var/Result_setting3.csv")
write.csv(Result_setting4, "two_var/Result_setting4.csv")


# colnames(Result_setting1) = col_names_list
# 
# 
# #Setting 2 Fix n, increase K
# KK2 = seq(2,92,by = 10) # different number of hospitals
# n = 1000 # each site has 1000 patients
# NN2 = n*KK2
# p2 = n/NN2
# Nsim = 500
# Result_setting2 = matrix(0,nrow =length(NN2) ,ncol =51 )
# for(i in 1:length(NN2)){
#   print(i)
#   out = replicate(Nsim,simu(beta_true,NN2[i],p2[i]))
#   set_out2 = apply(out,1,mean)
#   Result_setting2[i,] = set_out2
# }
# 
# colnames(Result_setting2) = col_names_list
# 
# #Setting 3 Fix N, increase K
# N3 = 10000
# KK3 = seq(2,92,by = 10)
# p3 = 1/KK3
# NN3 = rep(N3,10)
# Nsim = 500
# Result_setting3 = matrix(0,nrow =length(NN3) ,ncol =51 )
# for(i in 1:length(NN3)){
#   print(i)
#   out = replicate(Nsim,simu(beta_true,NN3[i],p3[i]))
#   set_out3 = apply(out,1,mean)
#   Result_setting3[i,] = set_out3
# }
# 
# colnames(Result_setting3) = col_names_list
# 
# #Setting 4 Fix N, fix K, increase p
# N4 = 10000
# #KK = rep(10,100,by = 10)
# p4 = seq(0.05,0.99,length=10)
# nn4 = p4*N4
# NN4 = rep(N4,10)
# Nsim = 500
# Result_setting4 = matrix(0,nrow =length(NN4) ,ncol =51 )
# for(i in 1:length(NN4)){
#   print(i)
#   out = replicate(Nsim,simu(beta_true,NN4[i],p4[i]))
#   set_out4 = apply(out,1,mean)
#   Result_setting4[i,] = set_out4
# }
# 
# colnames(Result_setting4) = col_names_list
# 
# 
# # save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/2_hete_homo.Rdata")
# # save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/2_hete_diff1.Rdata")
# # save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/2_hete_diff2.Rdata")
# save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/2_hete_diff3.Rdata")


