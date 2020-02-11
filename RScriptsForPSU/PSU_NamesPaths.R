
##############################################################
############### PSU_ALPHA #####################################
##############################################################

#####   USER SETTINGS   #####

## This script can be used to set the names and paths for the subsequent reports ##

RUNID = "test1" 
PDBFILE = "chains"
ROOTPATH = "F:\\PSUA\\Output\\"

##################################################################################
# The working directory and file names will be created from the above (edit as preferred)
workingdirectory = paste(ROOTPATH, RUNID,"\\Reports\\",sep="")
setwd(workingdirectory)

rama_report = paste(PDBFILE,"_torsion.txt",sep="")
calpha_report = paste(PDBFILE,"_calpha.txt",sep="")
rmsd_contact_report = paste(PDBFILE,"_rmsdcontact.txt",sep="")




