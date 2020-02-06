library(ggplot2)

setwd("C:/Users/Rache/OneDrive/dev/CPP_PDBStructuralViewer/data/")

############### USER INPUTS ##########################

pdbname <- "6j4a"
#pdbname <- "1lyz"
#pdbname <- "2qip"
#pdbname <- "2vb1"

############### USER INPUTS ##########################



#Load data fram created from C++ project
calphaname <- paste(pdbname,"_calpha.txt",sep="")
data <- read.csv(calphaname)
sub_data <- subset(data,Distance <  10) 

#CAlpha Distance Map
par(mfrow=c(2,2))

ggplot(sub_data, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste(pdbname,"C-Alpha Contact Map by PSU:RAlcraft")) + 
  geom_point()  


ggplot(data, aes(x=Id1, y=Id2,color=Distance)) + 
  ggtitle(paste(pdbname,"C-Alpha Contact Map by PSU:RAlcraft")) + 
  scale_colour_gradient(low = "darkgreen", high = "lightgrey", 
                        space = "Lab", na.value = "grey50", 
                        guide = "colourbar",  aesthetics = "colour") + 
  geom_point()  


  

