library(ggplot2)

#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################

#!!!!! Warning the heat map can be slow  !!!! ##

#################### USER INPUT ################################
# Enter angstrom cut off for contact map
angstrom = 9
# The data has been created with an arbirary cutoff of 25A
# This is nice for a heat map, but not for a more recognisable contact map
################################################################



data <- read.csv(calpha_report_PPI)
sub_data <- subset(data,Distance <  angstrom) 

#CAlpha Distance Map
ggplot(sub_data, aes(x=Chain1, y=Chain2,color=Amino1)) + 
  ggtitle(paste(PDBFILE,"C-Alpha Contact Map by PSU:Alpha A<",angstrom)) + 
  labs(x = chain1, y=chain2, color="Amino Acid")+
  geom_point()  

#CAlpha Heat Contact Map
ggplot(data, aes(x=Chain1, y=Chain2,color=Distance)) + 
  ggtitle(paste(PDBFILE,"C-Alpha PPI Contact Map by PSU-Alpha")) + 
  scale_colour_gradient(low = "darkgreen", high = "lightgrey", 
                        space = "Lab", na.value = "grey50", 
                        guide = "colourbar",  aesthetics = "colour") + 
  labs(x = chain1, y=chain2, color="Distance")+
  geom_point()  


  

