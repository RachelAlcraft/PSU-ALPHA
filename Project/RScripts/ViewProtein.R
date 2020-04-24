library(ggplot2)
library(gridExtra)
library(scatterplot3d)



#1) Set the working directory to collect the file with the tct version of a pdb
workingdirectory = "F:\\PSUA\\Code\\PSU-ALPHA\\Project\\RScripts\\"
setwd(workingdirectory)

#2) Get the file containing the pdb codes you want to download
filename="HTest.txt"
proteinfile = read.csv(filename, header = TRUE)

#3) Manipulate look
shapes = c(16, 17, 18,19)
shapes <- shapes[as.numeric(proteinfile$Residue)]

colors <- c("green", "blue", "yellow","red")
colors <- colors[as.numeric(proteinfile$Element)]

#4) Plot
s3d <- scatterplot3d(proteinfile$X, proteinfile$Y, proteinfile$Z,main="Hydrogens",pch=shapes,color=colors,box=FALSE)
legend(legend = levels(proteinfile$Residue), pch = c(16,17,18,19))