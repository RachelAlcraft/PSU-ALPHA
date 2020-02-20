library(ggplot2)
library(gridExtra)
#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################
data <- read.csv(rama_report)

sub_dataA <- subset(data,Chain == "A") 
sub_dataB <- subset(data,Chain == "B") 
sub_dataC <- subset(data,Chain == "C") 
sub_dataD <- subset(data,Chain == "D") 
sub_dataE <- subset(data,Chain == "E") 
sub_dataF <- subset(data,Chain == "F") 

#Phi/Psi
p1 <- ggplot(sub_dataA, aes(x=Phi, y=Psi,color=SS)) + 
  ggtitle(paste("1BVR","A")) + 
  geom_point(size=4) + 
  theme(legend.position="none") +
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))

p2 <- ggplot(sub_dataB, aes(x=Phi, y=Psi,color=SS)) + 
  ggtitle(paste("1BVR","B")) + 
  geom_point(size=4) + 
  theme(legend.position="none") +
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))

p3 <- ggplot(sub_dataC, aes(x=Phi, y=Psi,color=SS)) + 
  ggtitle(paste("1BVR","C")) + 
  geom_point(size=4) + 
  theme(legend.position="none") +
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))

p4 <- ggplot(sub_dataD, aes(x=Phi, y=Psi,color=SS)) + 
  ggtitle(paste("1BVR","D")) + 
  geom_point(size=4) + 
  theme(legend.position="none") + 
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))

p5 <- ggplot(sub_dataE, aes(x=Phi, y=Psi,color=SS)) + 
  ggtitle(paste("1BVR","E")) + 
  geom_point(size=4) + 
  theme(legend.position="none") + 
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))

p6 <- ggplot(sub_dataF, aes(x=Phi, y=Psi,color=SS)) + 
  ggtitle(paste("1BVR","F")) + 
  geom_point(size=4) + 
  scale_colour_discrete(labels = c("right a-helix","ideal beta","extended-glycine","glycine-helix","left a-helix","beta-p","unknown"))+
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))
  

grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2)