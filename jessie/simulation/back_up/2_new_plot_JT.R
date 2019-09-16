library(ggplot2)
library(reshape2)
library(grid)
library(Rmisc)
library(dplyr)
library(RColorBrewer)

# Real_beta = c(rep(beta_true[1],10),rep(beta_true[2],10),rep(beta_true[3],10),rep(beta_true[4],10),rep(beta_true[5],10))
# Real_beta <- matrix(Real_beta,nrow = 10)
cc = c(brewer.pal(9,"Set2")[c(5,3,2,1)],"pink","purple")
#cc = c("brown2",cc1[1:3])
##########  process the data#################
# four variables
R1_plot = cbind(Result_setting1[,3+1],Result_setting1[,8+1],Result_setting1[,13+1], Result_setting1[,18+1], Result_setting1[,23+1], Result_setting1[,28+1],
                Result_setting1[,33+1],Result_setting1[,38+1],Result_setting1[,43+1], Result_setting1[,48+1], Result_setting1[,53+1], Result_setting1[,58+1],
                Result_setting1[,68+1],Result_setting1[,78+1], nn1)
# three variables
# R1_plot = cbind(Result_setting1[,3+1],Result_setting1[,7+1],Result_setting1[,11+1], Result_setting1[,15+1], Result_setting1[,19+1], Result_setting1[,23+1],
#                 Result_setting1[,27+1],Result_setting1[,31+1],Result_setting1[,35+1], Result_setting1[,39+1], Result_setting1[,43+1], Result_setting1[,47+1],
#                 Result_setting1[,55+1],Result_setting1[,63+1], nn1)

# two variables
# R1_plot = cbind(Result_setting1[,3+1],Result_setting1[,6+1],Result_setting1[,9+1], Result_setting1[,12+1], Result_setting1[,15+1], Result_setting1[,18+1],
#                 Result_setting1[,21+1],Result_setting1[,24+1],Result_setting1[,27+1], Result_setting1[,30+1], Result_setting1[,33+1], Result_setting1[,36+1],
#                 Result_setting1[,42+1],Result_setting1[,48+1], nn1)

# relative bias: abs(beta - beta_pooled) / beta_pooled
# R1_plot[,1:6] = abs((R1_plot[,1:6] - cbind(R1_plot[,13],R1_plot[,13],R1_plot[,13],R1_plot[,13],R1_plot[,13],R1_plot[,13]))/R1_plot[,13])
R1_plot[,1:6] = R1_plot[,1:6] - matrix(beta_true[3], ncol = 6, nrow = 10)

# ratio of standrad error = sd(beta)/sd(beta_pooled)
# R1_plot[,7:12] = sqrt(R1_plot[,7:12])/sqrt(R1_plot[,14])
R1_plot[,7:12] = sqrt(R1_plot[,7:12])
temp = rep(0,15)
names(temp) =c("ODAL1 mean","ODAL1 median","ODAL2 mean mean", "ODAL2 mean median", "ODAL2 median mean", "ODAL2 median median",
               "ODAL1 mean","ODAL1 median","ODAL2 mean mean", "ODAL2 mean median", "ODAL2 median mean", "ODAL2 median median",
               "rl","rl_var","n")
temp = cbind(t(t(temp)),t(R1_plot))
R1_plot = t(temp)
R1_plot = R1_plot[2:11,]

DR1_PE <- as.data.frame(R1_plot[,c(1:6,15)])
DR1_PS <- as.data.frame(R1_plot[,c(7:12,15)])
QU1 = melt(DR1_PE,id.vars ="n")
QD1 = melt(DR1_PS,id.vars ="n")

colnames(QU1) = c("n","Method","Bias")
colnames(QD1) = c("n","Method","Std_Err")

PE1 <- ggplot(QU1,aes(x=n,y=Bias,group = Method,colour = Method,shape=Method))+geom_line(size =1) +
  geom_point(size=3) + labs(title = "Bias\n(beta_hat - beta_true)",tag = "Setting A",y ="") +
   scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18,19,20)) + 
  theme( axis.title=element_text(size=14),
         #legend.title=element_blank(),
         legend.position = "bottom",
         plot.title = element_text(face="bold",size = 20,hjust = 0.5),
         plot.tag = element_text(face="bold",size = 16,angle=90),
         plot.tag.position="left") +
  geom_hline(aes(yintercept=0), linetype="dashed",size=1)
PE1

PS1 <- ggplot(QD1,aes(x=n,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1)+geom_point(size=3)+labs(title = "Standard Error\n",y ="") +
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18,19,20)) +
  theme(legend.position="bottom",axis.title=element_text(size=14),legend.title=element_blank(),
        plot.title = element_text(face="bold",size = 20,hjust = 0.5))
multiplot(PE1,PS1, cols=2)
###############END Setting1####################

##########  process the data#################
R2_plot = cbind(Result_setting2[,8],Result_setting2[,23], Result_setting2[,28], Result_setting2[,33],
                Result_setting2[,18],Result_setting2[,38],Result_setting2[,43], Result_setting2[,48],
                Real_beta[, 2],KK2)

R2_plot[,1:4] = R2_plot[,1:4] - cbind(R2_plot[,9],R2_plot[,9],R2_plot[,9],R2_plot[,9])
R2_plot[,5:8] = sqrt(R2_plot[,5:8])/sqrt(R2_plot[,5])
temp = rep(0,10)
names(temp) =c("POOLED","ODAL1","ODAL2 mean", "ODAL2 median",
               "POOLED","ODAL1","ODAL2 mean", "ODAL2 median",
               "rl","K")
temp = cbind(t(t(temp)),t(R2_plot))
R2_plot = t(temp)
R2_plot = R2_plot[2:11,]

DR2_PE <- as.data.frame(R2_plot[,c(1:4,10)])
DR2_PS <- as.data.frame(R2_plot[,c(5:8,10)])
QU2 = melt(DR2_PE,id.vars ="K")
QD2 = melt(DR2_PS,id.vars ="K")

colnames(QU2) = c("K","Method","Bias")
colnames(QD2) = c("K","Method","Std_Err")

PE2 <- ggplot(QU2,aes(x=K,y=Bias,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(tag = "Setting B",y ="") + 
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14),
        plot.tag = element_text(face="bold",size = 16,angle=90),
        plot.tag.position="left") +
  geom_hline(aes(yintercept=0), linetype="dashed")
# PE2

PS2 <- ggplot(QD2,aes(x=K,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + labs(y ="") + geom_point(size=3) +
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14))
# PS2
###############END Setting2####################


##########  process the data#################
R3_plot = cbind(Result_setting3[,8],Result_setting3[,23], Result_setting3[,28], Result_setting3[,33],
                Result_setting3[,18],Result_setting3[,38],Result_setting3[,43], Result_setting3[,48],
                Real_beta[, 2],KK3)

R3_plot[,1:4] = R3_plot[,1:4] - cbind(R3_plot[,9],R3_plot[,9],R3_plot[,9],R3_plot[,9])
R3_plot[,5:8] = sqrt(R3_plot[,5:8])/sqrt(R3_plot[,5])
temp = rep(0,10)
names(temp) =c("POOLED","ODAL1","ODAL2 mean", "ODAL2 median",
               "POOLED","ODAL1","ODAL2 mean", "ODAL2 median",
               "rl","K")
temp = cbind(t(t(temp)),t(R3_plot))
R3_plot = t(temp)
R3_plot = R3_plot[2:11,]

DR3_PE <- as.data.frame(R3_plot[,c(1:4,10)])
DR3_PS <- as.data.frame(R3_plot[,c(5:8,10)])
QU3 = melt(DR3_PE,id.vars ="K")
QD3 = melt(DR3_PS,id.vars ="K")

colnames(QU3) = c("K","Method","Bias")
colnames(QD3) = c("K","Method","Std_Err")

PE3 <- ggplot(QU3,aes(x=K,y=Bias,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(tag = "Setting C",y ="") + 
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14),
        plot.tag = element_text(face="bold",size = 16,angle=90),
        plot.tag.position="left") +
  geom_hline(aes(yintercept=0), linetype="dashed")
# PE3

PS3 <- ggplot(QD3,aes(x=K,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(y ="") + geom_point(size=2) +
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14)) 
# PS3
###############END Setting3####################

##########  process the data#################
R4_plot = cbind(Result_setting4[,8],Result_setting4[,23], Result_setting4[,28], Result_setting4[,33],
                Result_setting4[,18],Result_setting4[,38],Result_setting4[,43], Result_setting4[,48],
                Real_beta[, 2],nn4)

R4_plot[,1:4] = R4_plot[,1:4] - cbind(R4_plot[,9],R4_plot[,9],R4_plot[,9],R4_plot[,9])
R4_plot[,5:8] = sqrt(R4_plot[,5:8])/sqrt(R4_plot[,5])
temp = rep(0,10)
names(temp) =c("POOLED","ODAL1","ODAL2 mean", "ODAL2 median",
               "POOLED","ODAL1","ODAL2 mean", "ODAL2 median",
               "rl","n")
temp = cbind(t(t(temp)),t(R4_plot))
R4_plot = t(temp)
R4_plot = R4_plot[2:11,]

DR4_PE <- as.data.frame(R4_plot[,c(1:4,10)])
DR4_PS <- as.data.frame(R4_plot[,c(5:8,10)])
QU4 = melt(DR4_PE,id.vars ="n")
QD4 = melt(DR4_PS,id.vars ="n")

colnames(QU4) = c("n","Method","Bias")
colnames(QD4) = c("n","Method","Std_Err")

PE4 <- ggplot(QU4,aes(x=n,y=Bias,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(tag = "Setting D",y ="") + 
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="bottom",axis.title=element_text(size=14),
        plot.tag = element_text(face="bold",size = 16,angle=90),
        plot.tag.position="left",legend.key = element_rect(size = 2),
        legend.key.size = unit(1, 'lines'),
        legend.title=element_blank()) +
  geom_hline(aes(yintercept=0), linetype="dashed")
# PE4

PS4 <- ggplot(QD4,aes(x=n,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(y ="") + geom_point(size=2) + 
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="bottom",axis.title=element_text(size=14),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(1, 'lines'),
        legend.title=element_blank()) 
# PS4
###############END Setting4####################

pdf("/Users/Jessie/Desktop/000_UPenn_Research/000_project/000_ODAL/R_code/jessie/simulation/2hetero_diff3_500.pdf", height = 18, width = 14)
multiplot(PE1,PE2,PE3,PE4,PS1,PS2,PS3,PS4, cols=2)
grid.newpage()
pushViewport(viewport(layout = grid.layout(140,51)))
print(PE1, vp = viewport(layout.pos.row = 3:39, layout.pos.col = 1:25))
print(PS1, vp = viewport(layout.pos.row = 3:39, layout.pos.col = 26:50))
print(PE2, vp = viewport(layout.pos.row = 41:70, layout.pos.col = 1:25))
print(PS2, vp = viewport(layout.pos.row = 41:70, layout.pos.col = 26:50))
print(PE3, vp = viewport(layout.pos.row = 72:101, layout.pos.col = 1:25))
print(PS3, vp = viewport(layout.pos.row = 72:101, layout.pos.col = 26:50))
print(PE4, vp = viewport(layout.pos.row = 103:137, layout.pos.col = 1:25))
print(PS4, vp = viewport(layout.pos.row = 103:137, layout.pos.col = 26:50))
dev.off()
