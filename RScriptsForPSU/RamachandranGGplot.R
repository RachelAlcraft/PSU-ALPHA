library(ggplot2)

setwd("C:/Users/Rache/OneDrive/dev/CPP_PDBStructuralViewer/data/")

############### USER INPUTS ##########################

pdbname <- "2vb1"

############### USER INPUTS ##########################

#Load data fram created from C++ project
ramaname <- paste(pdbname,"_torsion.txt",sep="")
data <- read.csv(ramaname)

#Phi/Psi
par(mfrow=c(1,1))

ggplot(data, aes(x=Phi, y=Psi,color=Chemical)) + 
  ggtitle("Ramachandran by PSU:RAlcraft") + 
  geom_point(size=4) + 
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))
  
#Chi1/Chi2
par(mfrow=c(1,1))

ggplot(data, aes(x=Chi1, y=Chi2,color=Chemical)) + 
  ggtitle("Chi report by PSU:RAlcraft") + 
  geom_point() + 
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))
