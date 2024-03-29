library(ggplot2)
library(reshape2)
library(grid)
library(Rmisc)
library(dplyr)
library(RColorBrewer)

cc = c(brewer.pal(9,"Set2")[c(3,1)],"pink","purple")
#cc = c("brown2",cc1[1:3])
##########  process the data#################
R1_plot = cbind(Result_setting1[,3],Result_setting1[,8],Result_setting1[,23],
                Result_setting1[,28],Result_setting1[,13],Result_setting1[,18]
                ,Result_setting1[,33],Result_setting1[,38],Real_beta[,2],nn1)
R1_plot[,1:4] = R1_plot[,1:4] - cbind(R1_plot[,9],R1_plot[,9],R1_plot[,9],R1_plot[,9])
R1_plot[,5:8] = sqrt(R1_plot[,5:8])
temp = rep(0,10)
names(temp) =c("LOCAL","POOLED","ODAL1","ODAL2","LOCAL","POOLED","ODAL1",
               "ODAL2","rl","n")
temp = cbind(t(t(temp)),t(R1_plot))
R1_plot = t(temp)
R1_plot = R1_plot[2:11,]

DR1_PE <- as.data.frame(R1_plot[,c(1:4,10)])
DR1_PS <- as.data.frame(R1_plot[,c(5:8,10)])
QU1 = melt(DR1_PE,id.vars ="n")
QD1 = melt(DR1_PS,id.vars ="n")



colnames(QU1) = c("n","Method","Bias")
colnames(QD1) = c("n","Method","Std_Err")

PE1 <- ggplot(QU1,aes(x=n,y=Bias,group = Method,colour = Method,shape=Method))+geom_line(size =1) +
  geom_point(size=3) + labs(title = "Bias\n",tag = "Setting A",y ="") +
   scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) + 
  theme( axis.title=element_text(size=14),
         #legend.title=element_blank(),
         legend.position = "none",
         plot.title = element_text(face="bold",size = 20,hjust = 0.5),
         plot.tag = element_text(face="bold",size = 16,angle=90),
         plot.tag.position="left") +
  geom_hline(aes(yintercept=0), linetype="dashed",size=1)
PS1 <- ggplot(QD1,aes(x=n,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1)+geom_point(size=3)+labs(title = "Standard Error\n",y ="") + 
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14),legend.title=element_blank(),
        plot.title = element_text(face="bold",size = 20,hjust = 0.5))
###############END Setting1####################

##########  process the data#################
R2_plot = cbind(Result_setting2[,3],Result_setting2[,8],Result_setting2[,23],
                Result_setting2[,28],Result_setting2[,13],Result_setting2[,18]
                ,Result_setting2[,33],Result_setting2[,38],Real_beta[,2],KK2)
R2_plot[,1:4] = R2_plot[,1:4] - cbind(R2_plot[,9],R2_plot[,9],R2_plot[,9],R2_plot[,9])
R2_plot[,5:8] = sqrt(R2_plot[,5:8])
temp = rep(0,10)
names(temp) =c("LOCAL","POOLED","ODAL1","ODAL2","LOCAL","POOLED","ODAL1","ODAL2","rl","K")
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

PS2 <- ggplot(QD2,aes(x=K,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + labs(y ="") + geom_point(size=3) +
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14))
###############END Setting2####################


##########  process the data#################
R3_plot = cbind(Result_setting3[,3],Result_setting3[,8],Result_setting3[,23],
                Result_setting3[,28],Result_setting3[,13],Result_setting3[,18]
                ,Result_setting3[,33],Result_setting3[,38],Real_beta[,2],KK3)
R3_plot[,1:4] = R3_plot[,1:4] - cbind(R3_plot[,9],R3_plot[,9],R3_plot[,9],R3_plot[,9])
R3_plot[,5:8] = sqrt(R3_plot[,5:8])
temp = rep(0,10)
names(temp) =c("LOCAL","POOLED","ODAL1","ODAL2","LOCAL","POOLED","ODAL1","ODAL2","rl","K")
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

PS3 <- ggplot(QD3,aes(x=K,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(y ="") + geom_point(size=2) +
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="none",axis.title=element_text(size=14)) 
###############END Setting3####################

##########  process the data#################
R4_plot = cbind(Result_setting4[,3],Result_setting4[,8],Result_setting4[,23],
                Result_setting4[,28],Result_setting4[,13],Result_setting4[,18]
                ,Result_setting4[,33],Result_setting4[,38],Real_beta[,2],nn4)
R4_plot[,1:4] = R4_plot[,1:4] - cbind(R4_plot[,9],R4_plot[,9],R4_plot[,9],R4_plot[,9])
R4_plot[,5:8] = sqrt(R4_plot[,5:8])
temp = rep(0,10)
names(temp) =c("LOCAL","POOLED","ODAL1","ODAL2","LOCAL","POOLED","ODAL1","ODAL2","rl","n")
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

PS4 <- ggplot(QD4,aes(x=n,y=Std_Err,group = Method,colour = Method,shape=Method)) + 
  geom_line(size=1) + geom_point(size=3) + labs(y ="") + geom_point(size=2) + 
  scale_color_manual(values = cc) + scale_shape_manual(values = c(16,17,15,18)) +
  theme(legend.position="bottom",axis.title=element_text(size=14),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(1, 'lines'),
        legend.title=element_blank()) 
###############END setting4####################

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
