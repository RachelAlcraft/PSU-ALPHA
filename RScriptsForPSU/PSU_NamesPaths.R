
##############################################################
############### PSU_ALPHA #####################################
##############################################################

#####   USER SETTINGS   #####

## This script can be used to set the names and paths for the subsequent reports ##

RUNID = "6J4A" 
PDBFILE = "6J4A"
ROOTPATH = "F:\\PSUA\\Output\\"

##################################################################################
# The working directory and file names will be created from the above (edit as preferred)
workingdirectory = paste(ROOTPATH, RUNID,"\\Reports\\",sep="")
setwd(workingdirectory)

chain1 = "A"
chain2 = "A"
rama_report = paste(PDBFILE,"_torsion.txt",sep="")
calpha_report = paste(PDBFILE,"_calpha.txt",sep="")
calpha_report_PPI = paste(PDBFILE,"_calphaChainA.txt",sep="")
rmsd_contact_report = paste(PDBFILE,"_rmsdcontact.txt",sep="")




