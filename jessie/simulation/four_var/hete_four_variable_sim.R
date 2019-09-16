library(mvtnorm)
library(ggplot2)
library(matlib)
library(meta)

##################################################
###Distributed algorithm for logistic regression##
### homogenrous
### setting one (K = 10, each n = seq(100, 1000, by=100), total N = 1000), 
##################################################
setwd("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation")
#########################
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
Generate = function(N,beta,nn,x1_mean, x1_var,x2_var,x3_par,x4_par,b0,b1,b3,b4){
  
  X1 = rnorm(N-nn)         #generate N numbers with N(0,1)
  X2 = rbinom(N-nn,1,0.3)  #generate N numbers with B(1,0.3)
  X3 = runif(N-nn,-1,1)    #generate N numbers with U(-1,1)
  # X4 = rnorm(N-nn,0,2)     #generate N numbers with N(0,2)
  X4 = rbinom(N-nn,1,0.5)
  
  X = cbind(1,X1,X2,X3,X4)
  # X = cbind(1,X1,X2,X3)
  meanY = expit(X%*%beta)
  Y = rbinom(N-nn,1,meanY)
  
  allMinisHetero = cbind(X1,X2,X3,X4,Y)
  # allMinisHetero = cbind(X1,X2,X3,Y)
  
  # the heterogeneous site
  X1.hetero = rnorm(nn,x1_mean,x1_var)
  X2.hetero = rbinom(nn,1,x2_var)
  X3.hetero = runif(nn,-x3_par,x3_par)
  # X4.hetero = rnorm(nn,1,2)
  X4.hetero = rbinom(nn,1,x4_par)
   
  # beta.hetero = c(1,1,-1,1,-1)
  beta.hetero = c(b0,b1,-1,b3,b4)
  # beta.hetero = c(1.5, 0.8, -1, 1.2, -1.3)
  # beta.hetero = c(1.5, 0.8, -1, 1.2)
  X.hetero = cbind(1,X1.hetero,X2.hetero,X3.hetero,X4.hetero)
  # X.hetero = cbind(1,X1.hetero,X2.hetero,X3.hetero)
  meanY.hetero = expit(X.hetero%*%beta.hetero)
  Y.hetero = rbinom(nn,1,meanY.hetero)
  
  Hetero = cbind(X1.hetero, X2.hetero, X3.hetero, X4.hetero, Y.hetero)
  # Hetero = cbind(X1.hetero, X2.hetero, X3.hetero, Y.hetero)

  
  # return(rbind(allMinisHetero))
  return(rbind(allMinisHetero, Hetero))
}

#simulation
simu = function(beta,N,p,nn,x1_mean, x1_var,x2_var,x3_par,x4_par,b0,b1,b3,b4){
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
  Data = Generate(N,beta,nn,x1_mean, x1_var,x2_var,x3_par,x4_par,b0,b1,b3,b4)
  # nlocal = round(p*N,1)   #round to 1e-1
  X = Data[,1:4]
  Y = Data[,5]
  Xlocal = Data[1:nn,1:4]
  Ylocal = Data[1:nn,5]
  X.homo = Data[1:(N-nn),1:4]
  Y.homo = Data[1:(N-nn),5]
  
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
  beta_meta_fix <-metagen(beta_meta_list, se_meta_list,comb.fixed = TRUE,comb.random = TRUE,prediction=TRUE,sm="OR")$TE.fixed
  beta_meta_random <-metagen(beta_meta_list, se_meta_list,comb.fixed = TRUE,comb.random = TRUE,prediction=TRUE,sm="OR")$TE.random
  
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
  
  # #output the outcome
  # total = c(iter, beta_tilde1_mean, beta_tilde1_median, beta_tilde2_mean_mean, beta_tilde2_mean_median, beta_tilde2_median_mean, beta_tilde2_median_median,
  #           var_tilde1_mean, var_tilde1_median, var_tilde2_mean_mean, var_tilde2_mean_median, var_tilde2_median_mean, var_tilde2_median_median, 
  #           beta0,beta_all,var0,var_all)
  total = c(iter,
            beta_tilde1_mean_homo[3], beta_tilde1_median_homo[3],
            beta_tilde1_mean[3], beta_tilde1_median[3],
            beta_tilde2_mean_mean[3], beta_tilde2_mean_median[3], 
            beta_tilde2_median_mean[3], beta_tilde2_median_median[3],
            beta_meta_fix,beta_meta_random,
            beta0[3],beta_all[3])
  # total = c(iter,beta0,beta_all,var0,var_all,
  #           beta_tilde1,beta_tilde2,beta_tilde2_median,
  #           var_tilde1,var_tilde2,var_tilde2_median)
  return(total)
}


# function to draw box plot
# num_var: number of variable used in the model
# df: data frame
# nn: number of patients in each site
# num_hete: number of hetero 
boxplot_func <- function(num_var,df.tmp,nn,num_hete,K1){
  
  print("I am plotting! Yeah!")
  
  col_names = c("iter",
                "ODAL1\nmean\n(w/o hetero site)", "ODAL1\nmedian\n(w/o hetero site)", 
                "ODAL1\nmean\n(All)","ODAL1\nmedian\n(All)",
                "ODAL2\nmean\nmean", "ODAL2\nmean\nmedian",
                "ODAL2\nmedian\nmean", "ODAL2\nmedian\nmedian",
                "Meta\nfixed\neffect","Meta\nrandom\neffect",
                "Local","Pooled")
  
  
  df.tmp = data.frame(t(df.tmp))
  names(df.tmp) = col_names
  
  df = df.tmp[,-c(1,12,13)]
  
  p <- ggplot(stack(df), aes(x=ind, y=values, fill=ind)) +
    geom_boxplot(notch=TRUE,fill=c("#db3340","#e9868f","#c66e06","#eba833",
                                   "#539e52","#8cc56d","#3f80ba","#95cbdb",
                                   "#9a64af","#cba0db")) +
    labs(title="",x="", y = "Estimated Log Odds Ratio") +
    geom_hline(yintercept =-1,linetype="dotted") + # the beta[x] here would be changed
    theme_classic(base_size=20) + 
    ggtitle(paste(num_var,"variable:", 
                  "Number of sites",K1, 
                  "& Number of patients in each site =",nn, 
                  "& Number of hetero site =",num_hete))
  
  print(p)
}

####Simulation Setting######
###Overall Setting##########
#1. number of covariates: 4, beta is of dimension 5#

####Simulation Setting######
###Overall Setting##########
#1. number of covariates: 4, beta is of dimension 5#
# True beta#
beta_true = c(1,1,-1,1,-1)
Nsim = 100

for (K1 in c(10,50,100)){
  
  nn1 = c(100, 500, 1000)   #number in a site
  NN1 = nn1*K1
  p1 = 0.1

  print(K1)
  
  # setting 1: one hetero
  set.seed(13)
  print("setting 1")
  # 1.1 setting: 100 each1 site
  out1_100 = replicate(Nsim,simu(beta_true,NN1[1],p1,nn = nn1[1],x1_mean=2, x1_var=1.5,x2_var=0.9,x3_par=1.5,x4_par=0.9,b0=1.6,b1=1.8,b3=2,b4=-0.5))  #run func with Nsim times
  # 1.2 setting: 500 each site
  out1_500 = replicate(Nsim,simu(beta_true,NN1[2],p1,nn = nn1[2],x1_mean=2, x1_var=1.5,x2_var=0.9,x3_par=1.5,x4_par=0.9,b0=1.6,b1=1.8,b3=2,b4=-0.5))  #run func with Nsim times
  # 1.3 setting: 1000 each site
  out1_1000 = replicate(Nsim,simu(beta_true,NN1[3],p1,nn = nn1[3],x1_mean=2, x1_var=1.5,x2_var=0.9,x3_par=1.5,x4_par=0.9,b0=1.6,b1=1.8,b3=2,b4=-0.5))  #run func with Nsim times
  
  
  # setting 2: two hetero
  set.seed(13)
  print("setting 2")
  # 2.1 setting: 100 each site
  out2_100 = replicate(Nsim,simu(beta_true,NN1[1],p1,nn = 2*nn1[1],x1_mean=2, x1_var=1.5,x2_var=0.9,x3_par=1.5,x4_par=0.9,b0=1.6,b1=1.8,b3=2,b4=-0.5))  #run func with Nsim times
  # 2.2 setting: 500 each site
  out2_500 = replicate(Nsim,simu(beta_true,NN1[2],p1,nn = 2*nn1[2],x1_mean=2, x1_var=1.5,x2_var=0.9,x3_par=1.5,x4_par=0.9,b0=1.6,b1=1.8,b3=2,b4=-0.5))  #run func with Nsim times
  # 2.3 setting: 1000 each site
  out2_1000 = replicate(Nsim,simu(beta_true,NN1[3],p1,nn = 2*nn1[3],x1_mean=2, x1_var=1.5,x2_var=0.9,x3_par=1.5,x4_par=0.9,b0=1.6,b1=1.8,b3=2,b4=-0.5))  #run func with Nsim times
  
  
  # boxplot
  pdf(paste("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/four_var/very_different_whole_four_var_results_","#sites_",K1,'.pdf',sep = ''), height = 10, width = 16)
  boxplot_func(4,out1_100,nn1[1],1,K1)
  boxplot_func(4,out1_500,nn1[2],1,K1)
  boxplot_func(4,out1_1000,nn1[3],1,K1)
  
  boxplot_func(4,out2_100,nn1[1],2,K1)
  boxplot_func(4,out2_500,nn1[2],2,K1)
  boxplot_func(4,out2_1000,nn1[3],2,K1)
  dev.off()
}



