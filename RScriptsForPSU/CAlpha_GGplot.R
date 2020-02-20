library(ggplot2)

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
sub_data <- subset(data,Distance <  angstrom)
sub_data <- subset(sub_data,Distance >  1) 

#CAlpha Distance Map
ggplot(sub_data, aes(x=Id1, y=Id2,color=SS1)) + 
  ggtitle(paste("1BVR","C-Alpha Contact Map A<",angstrom)) + 
  labs(x = "1BVR", y="1BVR", color="Secondary Structure")+
  scale_colour_discrete(labels = c("right a-helix","ideal beta","extended-glycine","glycine-helix","left a-helix","beta-p","unknown"))+
  geom_point()  

#CAlpha Heat Contact Map
ggplot(data, aes(x=Id1, y=Id2,color=Distance)) + 
  ggtitle(paste(PDBFILE,"C-Alpha Heat Contact Map by PSU:Alpha")) + 
  scale_colour_gradient(low = "darkgreen", high = "lightgrey", 
                        space = "Lab", na.value = "grey50", 
                        guide = "colourbar",  aesthetics = "colour") + 
  labs(x = PDBFILE, y=PDBFILE, color="Distance")+
  geom_point()  


  

