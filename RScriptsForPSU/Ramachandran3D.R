
library(scatterplot3d)

#setwd("C:/Users/Rache/Google Drive/ProjectPlan/Education/MScBio/Module_StructuralBio/Dev/")
setwd("C:/Users/Rache/OneDrive/dev/CPP_PDBStructuralViewer/data/")

#Load data fram created from C++ project
data <- read.csv("6j4a_rama.txt")

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
plot(data$Phi,data$Psi,xlab="Phi",ylab="Psi",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))
#text(data$Phi,data$Psi,xlab="Phi",ylab="Psi",labels=data$AminoAcid)
plot(data$Phi,data$Omega,xlab="Phi", ylab="Omega",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))
#text(data$Phi,data$Omega,xlab="Phi", ylab="Omega",labels=data$Id)
plot(data$Omega,data$Psi,xlab="Omega", ylab="Psi",col=colours,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180))
#text(data$Omega,data$Psi,xlab="Omega", ylab="Psi",labels=data$Id)
#Plot 3D phi against psi agasint iomega
scatterplot3d(data$Phi, data$Psi, data$Omega, color=colours,xlab="Phi",ylab="Psi",zlab="Omega")
