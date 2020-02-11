library(ggplot2)

#################################################################
# RELIES ON VALUES SET IN SCRIPT PSU_NamesPaths
#################################################################

#!!!!! Warning the heat map can be slow  !!!! ##

#################### USER INPUT ################################
# Enter angstrom cut off for contact map
angstrom = 12
# The data has been created with an arbirary cutoff of 25A
# This is nice for a heat map, but not for a more recognisable contact map
################################################################



data <- read.csv(rmsd_contact_report)
sub_data <- subset(data,Distance <  angstrom) 

#CAlpha Distance Map
ggplot(sub_data, aes(x=PDB1, y=PDB2,color=Hydro1)) + 
  ggtitle(paste(PDBFILE,"C-Alpha PPI Contact Map by PSU-Alpha")) + 
  geom_point()  

#CAlpha Heat Contact Map
ggplot(data, aes(x=PDB1, y=PDB2,color=Distance)) + 
  ggtitle(paste(PDBFILE,"C-Alpha PPI Contact Map by PSU-Alpha")) + 
  scale_colour_gradient(low = "darkgreen", high = "lightgrey", 
                        space = "Lab", na.value = "grey50", 
                        guide = "colourbar",  aesthetics = "colour") + 
  geom_point()  


  

