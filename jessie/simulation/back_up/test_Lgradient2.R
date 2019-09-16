expit <- function(x)
{
  exp(x)/(1+exp(x))
}

x1.1 <- c(1,2,3); x1.2 <- c(1,0,1);
x1 <- cbind(1, x1.1, x1.2)

x2.1 <- c(4,5,6); x2.2 <- c(0,1,0);
x2 <- cbind(1, x2.1, x2.2)

x3.1 <- c(7,8,9); x3.2 <- c(1,0,1);
x3 <- cbind(1, x3.1, x3.2)

x.comb <- rbind(x1, x2, x3)

beta = c(-0.5, 1, 1)

result1 <- (t(c(expit(x1 %*% beta)*(1+expit(x1 %*% beta))) * x1) %*% x1 / dim(x1)[1]  +
  t(c(expit(x2 %*% beta)*(1+expit(x2 %*% beta))) * x2) %*% x2 / dim(x2)[1] + 
  t(c(expit(x3 %*% beta)*(1+expit(x3 %*% beta))) * x3) %*% x3 / dim(x3)[1])/3

result2 <- t(c(expit(x.comb %*% beta)*(1+expit(x.comb %*% beta))) * x.comb) %*% x.comb / dim(x.comb)[1]


split(x.comb, (0:nrow(x.comb) %/% 3))

as.data.frame(x.comb)
split(x.comb, sample(1:3, nrow(x.comb), replace=F))

x = runif(24)
ind = sample(rep(1:4,each = length(x)/4))
split(x,ind)

x = seq(1:length(1:dim(x.comb)[1]))
ind = sample(rep(1:3, each = length(x)/3))
sp_ind = split(x,ind)
x1 = x.comb[sp_ind$`1`,]

A <- matrix(c(1:9), 3, 3)  
B <- matrix(c(2:10), 3, 3) 
C <- matrix(c(11:19), 3, 3) 
my.list <- list(A, B, C)
apply(simplify2array(my.list), 1:2, median)
