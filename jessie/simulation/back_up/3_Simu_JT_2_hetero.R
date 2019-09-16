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
col_names_list = c("iter", "local beta0", "local beta1", "local beta2", "local beta3","local beta4",
                   "all beta0", "all beta1", "all beta2", "all beta3", "all beta4", 
                   "local var0", "local var1", "local var2", "local var3", "local var4",
                   "all var0", "all var1", "all var2", "all var3", "all var4", 
                   "odal1 beta0", "odal1 beta1", "odal1 beta2", "odal1 beta3", "odal1 beta4", 
                   "odal2 beta0 mean", "odal2 beta1 mean", "odal2 beta2 mean", "odal2 beta3 mean", "odal2 beta4 mean", 
                   "odal2 beta0 med", "odal2 beta1 med", "odal2 beta2 med", "odal2 beta3 med", "odal2 beta3 med", 
                   "odal1 var0", "odal1 var1", "odal1 var2", "odal1 var3", "odal1 var4", 
                   "odal2 var0 mean", "odal2 var1 mean", "odal2 var2 mean", "odal2 var3 mean", "odal2 var4 mean", 
                   "odal2 var0 med", "odal2 var1 med", "odal2 var2 med", "odal2 var3 med", "odal2 var4 med")


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
Lgradient = function(beta,X,Y){
  design = cbind(1,X)
  t(Y-expit(design%*%t(t(beta))))%*%design/length(Y)
}

Lgradient2 =function(beta,X){
  design = cbind(1,X)
  Z=expit(design%*%beta)
  t(c(-1*Z*(1-Z))*design)%*%design/nrow(X)
}

# k is the number of hospitals
Lgradient2_median =function(beta,X,k){
  l = seq(1:length(1:dim(X)[1]))
  ind = sample(rep(1:k, each = length(l)/k))
  sp_ind = split(l,ind)

  mat = list()
  for (i in 1:k){
    id = sp_ind[[i]]
    X_single = X[id,]
    design = cbind(1,X_single)
    Z=expit(design%*%beta)
    mat[[i]] = t(c(-1*Z*(1-Z))*design)%*%design/nrow(X_single)
  }
  final_mat = apply(simplify2array(mat), 1:2, median)
  final_mat
}

Lgradient_meat = function(beta,X,Y){
  design = cbind(1,X)
  Z = expit(design%*%beta)
  t(c(Y-Z)*design)%*%(c(Y-Z)*design)/nrow(X)
}

#surogate likelihood, suppose the local data are stored in Xlocal,Ylocal
SL = function(beta){
  -Lik(beta,Xlocal,Ylocal) - L%*%beta
}

SL2 = function(beta){
  beta_temp = beta - beta0
  #mat_temp = t(t(beta_temp))%*%t(beta_temp)
  -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2%*%t(t(beta_temp))
}

SL2_median = function(beta){
  beta_temp = beta - beta0
  #mat_temp = t(t(beta_temp))%*%t(beta_temp)
  -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2_median%*%t(t(beta_temp))
}

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
Generate = function(N,beta,nn){
  
  X1 = rnorm(N-2*nn)         #generate N numbers with N(0,1)
  X2 = rbinom(N-2*nn,1,0.3)  #generate N numbers with B(1,0.3)
  X3 = runif(N-2*nn,X2-1,1)    #generate N numbers with U(-1,1)
  X4 = rnorm(N-2*nn,0,2)     #generate N numbers with N(0,2)
  X = cbind(1,X1,X2,X3,X4)
  meanY = expit(X%*%beta)
  Y = rbinom(N-2*nn,1,meanY)
  
  allMinisHetero = cbind(X1,X2,X3,X4,Y)
  
  # the heterogeneous site one
  X1.hetero = rnorm(nn,0.1,1) 
  X2.hetero = rbinom(nn,1,0.4)  #generate N numbers with B(1,0.5)
  X3.hetero = runif(nn,X2.hetero-1.1,1.1)    #generate N numbers with U(-2,2)
  X4.hetero = rnorm(nn,0.1,2)     #generate N numbers with N(0,1)
  X.hetero = cbind(1,X1.hetero,X2.hetero,X3.hetero,X4.hetero)
  meanY.hetero = expit(X.hetero%*%beta)
  Y.hetero = rbinom(nn,1,meanY.hetero)
  
  Hetero = cbind(X1.hetero, X2.hetero, X3.hetero, X4.hetero, Y.hetero)
  
  # the heterogeneous cite
  X1.hetero.2 = rnorm(nn,-0.1,1) 
  X2.hetero.2 = rbinom(nn,1,0.5)  #generate N numbers with B(1,0.4)
  X3.hetero.2 = runif(nn,X2.hetero.2-0.9,0.9)    #generate N numbers with U(-2,2)
  X4.hetero.2 = rnorm(nn,-0.1,2)     #generate N numbers with N(0,3)
  X.hetero.2 = cbind(1,X1.hetero.2,X2.hetero.2,X3.hetero.2,X4.hetero.2)
  meanY.hetero.2 = expit(X.hetero.2%*%beta)
  Y.hetero.2 = rbinom(nn,1,meanY.hetero.2)
  
  Hetero.2 = cbind(X1.hetero.2, X2.hetero.2, X3.hetero.2, X4.hetero.2, Y.hetero.2)
    
  return(rbind(allMinisHetero, Hetero, Hetero.2))
}

#simulation 1
simu = function(beta,N,p,k){
  SL = function(beta){
    -Lik(beta,Xlocal,Ylocal) - L%*%beta
  }
  SL2 = function(beta){
    beta_temp = beta - beta0
    #mat_temp = t(t(beta_temp))%*%t(beta_temp)
    -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2%*%t(t(beta_temp))
  }
  SL2_median = function(beta){
    beta_temp = beta - beta0
    #mat_temp = t(t(beta_temp))%*%t(beta_temp)
    -Lik(beta,Xlocal,Ylocal) - L%*%beta - 0.5*t(beta_temp)%*%L2_median%*%t(t(beta_temp))
  }
  #generate the data
  Data = Generate(N,beta,nn = N/k)
  nlocal = round(p*N,1)   #sishewuru to 1e-1
  X = Data[,1:4]
  Y = Data[,5]
  Xlocal = Data[1:nlocal,1:4]
  Ylocal = Data[1:nlocal,5]
  
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
  
  iter = fit$iter+1
  
  #estimatior beta and var
  L = Lgradient(beta0,X,Y)
  L2 = Lgradient2(beta0,X) - Lgradient2(beta0,Xlocal)
  L2_median = Lgradient2_median(beta0,X,k) - Lgradient2(beta0,Xlocal)
  
  #mean of the estimator
  beta_tilde1 = optim(beta0,SL,control = list(maxit = 10000,reltol = 1e-10))$par
  
  beta_tilde2 = optim(beta_tilde1,SL2,control = list(maxit = 10000,reltol = 1e-10))$par
  
  beta_tilde2_median = optim(beta_tilde1,SL2_median,control = list(maxit = 10000,reltol = 1e-10))$par
  # View(beta_tilde1)
  # View(beta_tilde2)
  # var_tilde1 =c(1,1,1,1,1)
  # var_tilde2 =c(1,1,1,1,1)
  
  #var of the estimator
  var_tilde1_mat = Sandwich(beta_tilde1,Xlocal,Ylocal,N)
  var_tilde2_mat = Sandwich(beta_tilde2,Xlocal,Ylocal,N)
  var_tilde2_mat_median = Sandwich(beta_tilde2_median,Xlocal,Ylocal,N)
  
  var_tilde1 = diag(var_tilde1_mat)
  var_tilde2 = diag(var_tilde2_mat)
  var_tilde2_median = diag(var_tilde2_mat_median)
  
  #output the outcome
  total = c(iter,beta0,beta_all,var0,var_all,beta_tilde1,beta_tilde2,beta_tilde2_median,var_tilde1,var_tilde2,var_tilde2_median)
  return(total)
}



####Simulation Setting######
###Overall Setting##########
#1. number of covariates: 4, beta is of dimension 5#

# True beta#
set.seed(1111)
beta_true = c(1,1,-1,1,-1)

#Setting 1 Fix K, increase n
K1 = 10
nn1 = seq(100,1000,by=100)   #number in a site
NN1 = nn1*K1
p1 = 0.1
Nsim = 400
Result_setting1 = matrix(0,nrow =length(NN1) ,ncol =51 )
for(i in 1:length(NN1)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN1[i],p1,K1))  #run func with Nsim times
  set_out1 = apply(out,1,mean)
  Result_setting1[i,] = set_out1
}

colnames(Result_setting1) = col_names_list


#Setting 2 Fix n, increase K
KK2 = seq(2,92,by = 10) # different number of hospitals
n = 1000 # each site has 1000 patients
NN2 = n*KK2
p2 = n/NN2
Nsim = 400
Result_setting2 = matrix(0,nrow =length(NN2) ,ncol =51 )
for(i in 1:length(NN2)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN2[i],p2[i],KK2[i]))
  set_out2 = apply(out,1,mean)
  Result_setting2[i,] = set_out2
}

colnames(Result_setting2) = col_names_list

#Setting 3 Fix N, increase K
N3 = 10000
KK3 = seq(2,92,by = 10)
p3 = 1/KK3
NN3 = rep(N3,10)
Nsim = 400
Result_setting3 = matrix(0,nrow =length(NN3) ,ncol =51 )
for(i in 1:length(NN3)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN3[i],p3[i],KK3[i]))
  set_out3 = apply(out,1,mean)
  Result_setting3[i,] = set_out3
}

colnames(Result_setting3) = col_names_list

#Setting 4 Fix N, fix K, increase p
N4 = 10000
#KK = rep(10,100,by = 10)
p4 = seq(0.05,0.99,length=10)
nn4 = p4*N4
NN4 = rep(N4,10)
Nsim = 400
Result_setting4 = matrix(0,nrow =length(NN4) ,ncol =51 )
for(i in 1:length(NN4)){
  print(i)
  out = replicate(Nsim,simu(beta_true,NN4[i],p4[i],k=10))
  set_out4 = apply(out,1,mean)
  Result_setting4[i,] = set_out4
}

colnames(Result_setting4) = col_names_list


save.image("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/3_hete.Rdata")






