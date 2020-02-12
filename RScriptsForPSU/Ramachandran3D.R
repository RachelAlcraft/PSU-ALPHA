
library(scatterplot3d)

#setwd("C:/Users/Rache/Google Drive/ProjectPlan/Education/MScBio/Module_StructuralBio/Dev/")
setwd("C:/Users/Rache/OneDrive/dev/CPP_PDBStructuralViewer/data/")

#Load data fram created from C++ project
data <- read.csv(rama_report)

sub_data_nochi1 <- subset(data,Chi1 !=0) 
sub_data_nochi2 <- subset(data,Chi2 !=0) 
sub_data_nochi3 <- subset(data,Chi3 !=0) 
sub_data_nochi4 <- subset(data,Chi4 !=0) 

#create the colours for the points based on something or other
colours<-vector()
for (i in 1:length(data$AminoAcid)){
  aa <- data$AminoAcid[i]
  # some decision
  col <- "black"
  if (grepl(aa, "ALA:VAL:LEU:ILE") ){col <- "yellow"}
  if (grepl(aa, "PHE:TYR:TRP") ){col <- "orange"}
  if (grepl(aa, "HIS:LYS:ARG") ){col <- "lightblue"}
  if (grepl(aa, "ASP:GLU") ){col <- "pink"}
  if (grepl(aa, "GLY:PRO") ){col <- "darkgrey"}
  if (grepl(aa, "SER:THR:CYS:MET:ASN:GLN") ){col <- "darkgreen"}
  colours <- c(colours,col)
}

# now 4 smaller comparisons with omega and a 3d view
par(mfrow=c(2,2)+.1)
par(mar=c(5,6,4,1)+.1)
plot(sub_data_nochi1$Phi,sub_data_nochi1$Chi1,xlab="Phi",ylab="Chi1",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))
plot(sub_data_nochi2$Phi,sub_data_nochi2$Omega,xlab="Phi", ylab="Chi2",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))
plot(sub_data_nochi3$Omega,sub_data_nochi3$Psi,xlab="Phi", ylab="Chi3",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))
plot(sub_data_nochi4$Omega,sub_data_nochi4$Psi,xlab="Phi", ylab="Chi4",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))


#Plot 3D phi against psi agasint iomega
par(mfrow=c(1,1)+.1)
scatterplot3d(data$Phi, data$Psi, data$Omega, color=colours,xlab="Phi",ylab="Psi",zlab="Omega",main="Phi-Psi-Omega")
