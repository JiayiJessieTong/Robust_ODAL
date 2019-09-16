Generate = function(N,beta,nn,x1_mu,x1_var,x2_var,b0,b1){
  
  X1 = rnorm(N-nn,0,1)         #generate N numbers with N(0,1)
  X2 = rbinom(N-nn,1,0.3)  #generate N numbers with B(1,0.3)
  
  X = cbind(1,X1,X2)
  meanY = expit(X%*%beta)
  Y = rbinom(N-nn,1,meanY)
  
  allMinisHetero = cbind(X1,X2,Y)
  
  # the heterogeneous site
  X1.hetero = rnorm(nn,x1_mu,x1_var)
  X2.hetero = rbinom(nn,1,x2_var)
  
  beta.hetero = c(b0, b1, -1)
  X.hetero = cbind(1,X1.hetero,X2.hetero)
  meanY.hetero = expit(X.hetero%*%beta.hetero)
  Y.hetero = rbinom(nn,1,meanY.hetero)
  
  Hetero = cbind(X1.hetero, X2.hetero, Y.hetero)
  
  
  return(rbind(allMinisHetero, Hetero))
}

N = 10000
beta = c(-5, 1, -1)
x1_mu = 1
x1_var = 1
x2_var = 0.4
b0 = -3
b1 = 1.2
nn = 1000
Data = Generate(N,beta,nn,x1_mu,x1_var,x2_var,b0,b1)

Xall = Data[,c(1,2)]
Yall = Data[,c(3)]
site = rep(c(1:10),each = 1000)  - 1
