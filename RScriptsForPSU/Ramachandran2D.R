
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


par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.2)

#Plot phi against psi in a large graph
plot(data$Phi,data$Psi,xlab="Phi",ylab="Psi",col=colours,tck=1,pch=19, cex=2,xlim=c(-180,180),ylim=c(-180,180),main="Ferritin Ramachandran (PSU RAlcraft)")
#plot(data$Phi,data$Psi,xlab="Phi",ylab="Psi",col=colours,xaxt = "n")
legend("topright", legend=c("Aliphatic", "Aromatic","Basic","Acidic","Confo","Neutral"),col=c("yellow", "orange","lightblue","pink","darkgrey","darkgreen"), lty=0:0, cex=1.0,pch=19,title="Amino Acid Colours")
#axis(side = 1, at = c(-180,-90,0,90,180),labels = c("-180","-90","0","90","180"),tck=-.02)
#axis(side = 2, lwd=0,at = c(-180,-90,0,90,180))

