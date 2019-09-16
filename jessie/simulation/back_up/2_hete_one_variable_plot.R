library(ggplot2)
library(reshape2)
library(grid)
library(Rmisc)
library(dplyr)
library(RColorBrewer)

# Real_beta = c(rep(beta_true[1],10),rep(beta_true[2],10),rep(beta_true[3],10),rep(beta_true[4],10),rep(beta_true[5],10))
# Real_beta <- matrix(Real_beta,nrow = 10 )
cc = c(brewer.pal(9,"Set2")[c(6:1)],"pink","purple")
#cc = c("brown2",cc1[1:3])
##########  process the data#################
# R1_plot = cbind(Result_setting1[,3],Result_setting1[,5],Result_setting1[,7], Result_setting1[,9], Result_setting1[,11], Result_setting1[,13],
#                 Result_setting1[,15],Result_setting1[,17],Result_setting1[,19], Result_setting1[,21], Result_setting1[,23], Result_setting1[,25],
#                 Result_setting1[,29],Result_setting1[,33], nn1)
R1_plot = cbind(Result_setting1[,3],Result_setting1[,5],Result_setting1[,7], 
                Result_setting1[,9], Result_setting1[,11], Result_setting1[,13],
                Result_setting1[,15],Result_setting1[,17],Result_setting1[,22],nn1)

# relative bias: abs(beta - beta_pooled) / beta_pooled
# R1_plot[,1:6] = abs((R1_plot[,1:6] - cbind(R1_plot[,13],R1_plot[,13],R1_plot[,13],R1_plot[,13],R1_plot[,13],R1_plot[,13]))/R1_plot[,13])
R1_plot[,1:9] = R1_plot[,1:9] - matrix(beta_true[2], nrow = 10, ncol = 9)

# ratio of standrad error = sd(beta)/sd(beta_pooled)
# R1_plot[,7:12] = sqrt(R1_plot[,7:12])/sqrt(R1_plot[,14])
# R1_plot[,7:12] = sqrt(R1_plot[,7:12])
# temp = rep(0,15)
# names(temp) =c("ODAL1 mean","ODAL1 median","ODAL2 mean mean", "ODAL2 mean median", "ODAL2 median mean", "ODAL2 median median",
#                "ODAL1 mean","ODAL1 median","ODAL2 mean mean", "ODAL2 mean median", "ODAL2 median mean", "ODAL2 median median",
#                "rl","rl_var","n")
temp = rep(0,10)
names(temp) =c("ODAL1 mean #1-9","ODAL1 median #1-9","ODAL1 mean","ODAL1 median", 
               "ODAL2 mean mean", "ODAL2 mean median", "ODAL2 median mean", "ODAL2 median median",
               "Meta","n")

temp = cbind(t(t(temp)),t(R1_plot))
R1_plot = t(temp)
R1_plot = R1_plot[2:11,]

DR1_PE <- as.data.frame(R1_plot[,c(1:10)])
# DR1_PS <- as.data.frame(R1_plot[,c(7:12,15)])
QU1 = melt(DR1_PE,id.vars ="n")
# QD1 = melt(DR1_PS,id.vars ="n")

colnames(QU1) = c("n","Method","Bias")
# colnames(QD1) = c("n","Method","Std_Err")

PE1 <- ggplot(QU1,aes(x=n,y=Bias,group = Method,colour = Method,shape=Method))+geom_line(size =1) +
  geom_point(size=3) + labs(title = "Bias\n(beta1_hat - beta1_true)",tag = "Setting A",y ="") +
  scale_color_manual(values = cc) + scale_shape_manual(values = c(15:22)) + 
  theme( axis.title=element_text(size=14),
         #legend.title=element_blank(),
         legend.position = "bottom",
         plot.title = element_text(face="bold",size = 20,hjust = 0.5),
         plot.tag = element_text(face="bold",size = 16,angle=90),
         plot.tag.position="left") +
  geom_hline(aes(yintercept=0), linetype="dashed",size=1)
PE1

# PS1 <- ggplot(QD1,aes(x=n,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
#   geom_line(size=1)+geom_point(size=3)+labs(title = "Standard Error\n",y ="") +
#   scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18,19,20)) +
#   theme(legend.position="bottom",axis.title=element_text(size=14),legend.title=element_blank(),
#         plot.title = element_text(face="bold",size = 20,hjust = 0.5))
# multiplot(PE1,PS1, cols=2)
###############END Setting1####################
