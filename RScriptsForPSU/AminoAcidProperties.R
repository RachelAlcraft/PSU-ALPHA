
library(ggplot2)
library(gridExtra)

workingdirectory = paste("F:\\PSUA\\Code\\PSU-ALPHA\\Config\\")
setwd(workingdirectory)

amino_config = "data_aminoinfo.csv"

data <- read.csv(amino_config)



xname = colnames(data)[1]
colnames(data)[1] = "AminoCode" #For some reason the config file has a bad character



p1 <- ggplot(data, aes(x=AminoCode, y=Hydro, fill=Chemical)) + 
  ggtitle(paste("Amino Acid Properties")) + 
  labs(x = "Amino Acid", y="Hydrophobicity", color="Chemical")+
  geom_bar(stat="identity")

p2 <- ggplot(data, aes(x=AminoCode, y=Donicity, fill=Polar)) + 
  ggtitle(paste("Amino Acid Properties")) + 
  labs(x = "Amino Acid", y="Donicity", color="Polar")+
  geom_bar(stat="identity")

grid.arrange(p1,p2,nrow=2)

