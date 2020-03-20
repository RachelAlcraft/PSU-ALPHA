library(ggplot2)
library(gridExtra)

#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################

#!!!!! Warning the heat map can be slow  !!!! ##

#################### USER INPUT ################################
# Enter angstrom cut off for contact map
angstrom = 12
# The data has been created with an arbitrary cutoff of 25A
# This is nice for a heat map, but not for a more recognisable contact map
################################################################



data <- read.csv(calpha_report)
data <- subset(data,Distance < angstrom) 

sub_dataA <- subset(data,Chain1 == "A") 
sub_dataA <- subset(sub_dataA,Chain2 == "A")

sub_dataB <- subset(sub_dataA,Id1 > 140) 
sub_dataB <- subset(sub_dataB,Id1 < 170) 
sub_dataB <- subset(sub_dataB,Id2 > 140) 
sub_dataB <- subset(sub_dataB,Id2 < 170)

sub_dataC <- subset(data,Id1 > 180) 
sub_dataC <- subset(sub_dataC,Id1 < 230) 
sub_dataC <- subset(sub_dataC,Id2 > 180) 
sub_dataC <- subset(sub_dataC,Id2 < 230) 


#sub_dataB <- subset(data.Chain1 == "B")# && data.Chain2 == "B") 
#sub_dataC <- subset(data.Chain1 == "C" && data.Chain2 == "C") 
#sub_dataD <- subset(data.Chain1 == "D" && data.Chain2 == "D") 
#sub_dataE <- subset(data.Chain1 == "E" && data.Chain2 == "E") 
#sub_dataF <- subset(data.Chain1 == "F" && data.Chain2 == "F") 

#CAlpha Distance Map
pall <- ggplot(data, aes(x=Id1, y=Id2,color=SS1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Maps")) + 
  scale_colour_discrete(labels = c("right a-helix","ideal beta","extended-glycine","glycine-helix","left a-helix","beta-p","unknown"))+
  labs(x = '1BVR', y='1BVR', color="SecStruc")+
  geom_point()  
p1 <- ggplot(sub_dataA, aes(x=Id1, y=Id2,color=SS1)) + 
  ggtitle(paste("1BVR","Chain A")) + 
  scale_colour_discrete(labels = c("right a-helix","ideal beta","extended-glycine","glycine-helix","left a-helix","beta-p","unknown"))+
  labs(x = '1BVR A', y='1BVR A', color="SecStruc")+
  geom_point()
p2 <- ggplot(sub_dataB, aes(x=Id1, y=Id2,color=SS1)) + 
  ggtitle(paste("1BVR","Chain A Residues 140-170")) + 
  scale_colour_discrete(labels = c("right a-helix","ideal beta","extended-glycine","glycine-helix","left a-helix","beta-p","unknown"))+
  labs(x = '1BVR 140-170', y='1BVR 140-170', color="SecStruc")+
  geom_point()
p3 <- ggplot(sub_dataC, aes(x=Id1, y=Id2,color=Hydro1)) + 
  ggtitle(paste("Model2","Chain A Residues 180-230")) + 
  labs(x = '180-230', y='180-230', color="Hydrophobicity")+
  geom_point()
#p2 <- ggplot(sub_dataB, aes(x=Id1, y=Id2,color=Chemical1)) + 
#  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
#  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
#  geom_point()  
#p3 <- ggplot(sub_dataC, aes(x=Id1, y=Id2,color=Chemical1)) + 
#  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
#  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
#  geom_point()  
#p4 <- ggplot(sub_dataD, aes(x=Id1, y=Id2,color=Chemical1)) + 
#  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
#  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
#  geom_point()  
#p5 <- ggplot(sub_dataE, aes(x=Id1, y=Id2,color=Chemical1)) + 
#  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
#  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
#  geom_point()  
#p6 <- ggplot(sub_dataF, aes(x=Id1, y=Id2,color=Chemical1)) + 
#  ggtitle(paste("1BVR","C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
#  labs(x = PDBFILE, y=PDBFILE, color="Chemical Type")+
#  geom_point()  



grid.arrange(p3,nrow=2)
  

