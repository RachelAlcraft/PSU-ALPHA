library(ggplot2)
library(gridExtra)

#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################

#!!!!! Warning the heat map can be slow  !!!! ##

#################### USER INPUT ################################
# Enter angstrom cut off for contact map
angstrom = 9
# The data has been created with an arbitrary cutoff of 25A
# This is nice for a heat map, but not for a more recognisable contact map
################################################################



data <- read.csv(calpha_report)
sub_dataA <- subset(data,Chain == "A") 
sub_dataB <- subset(data,Chain == "B") 
sub_dataC <- subset(data,Chain == "C") 
sub_dataD <- subset(data,Chain == "D") 
sub_dataE <- subset(data,Chain == "E") 
sub_dataF <- subset(data,Chain == "F") 

#CAlpha Distance Map
p1 <- ggplot(sub_dataA, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
  geom_point()
p2 <- ggplot(sub_dataB, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
  geom_point()  
p3 <- ggplot(sub_dataC, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
  geom_point()  
p4 <- ggplot(sub_dataD, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
  geom_point()  
p5 <- ggplot(sub_dataE, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
  geom_point()  
p6 <- ggplot(sub_dataF, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
  geom_point()  



grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2)
  

