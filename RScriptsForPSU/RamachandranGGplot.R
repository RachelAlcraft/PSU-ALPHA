library(ggplot2)

#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################
data <- read.csv(rama_report)

#Phi/Psi
ggplot(data, aes(x=Phi, y=Psi,color=Chemical)) + 
  ggtitle(paste(PDBFILE,"Ramachandran by PSU:Alpha")) + 
  geom_point(size=4) + 
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))
  
#Chi1/Chi2
par(mfrow=c(1,1))

ggplot(data, aes(x=Chi1, y=Chi2,color=Chemical)) + 
  ggtitle(paste(PDBFILE,"Chi report by PSU:Alpha")) + 
  geom_point() + 
  scale_y_continuous(breaks = c(-180,-120,-60,0,60,120,180))+
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180))
