library(ggplot2)

#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################

#!!!!! Warning the heat map can be slow  !!!! ##

#################### USER INPUT ################################
# Enter angstrom cut off for contact map
angstrom = 10
# The data has been created with an arbirary cutoff of 25A
# This is nice for a heat map, but not for a more recognisable contact map
################################################################



data <- read.csv(calpha_report)
sub_data <- subset(data,Distance <  angstrom) 

#CAlpha Distance Map
ggplot(sub_data, aes(x=Id1, y=Id2,color=Chemical1)) + 
  ggtitle(paste(PDBFILE,"C-Alpha Contact Map by PSU:RAlcraft")) + 
  geom_point()  

#CAlpha Heat Contact Map
ggplot(data, aes(x=Id1, y=Id2,color=Distance)) + 
  ggtitle(paste(PDBFILE,"C-Alpha Heat Contact Map by PSU:RAlcraft")) + 
  scale_colour_gradient(low = "darkgreen", high = "lightgrey", 
                        space = "Lab", na.value = "grey50", 
                        guide = "colourbar",  aesthetics = "colour") + 
  geom_point()  


  

