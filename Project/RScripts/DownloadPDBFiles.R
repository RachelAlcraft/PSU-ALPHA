



##################################################################################

##################################################################################

#1) Set the working directory to collect the file with the pdb codes in
workingdirectory = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-HighResList\\"
setwd(workingdirectory)

#2) Get the file containing the pdb codes you want to download
filename="highres_and_unique.txt"
pdbcodes = read.csv(filename, header = TRUE)


#3) Now set the workling directory again to the destination database directory
workingdirectory = "F:\\PSUA\\ProjectData\\HighResFiles\\"
setwd(workingdirectory)

#4) Loop through the codes to get each one

for (i in 1:length(pdbcodes$PDB)){

  pdb = paste(pdbcodes$PDB[i])
  url = paste("https://files.rcsb.org/download/",pdb,".pdb",sep="")
  destfile = paste(pdb,".pdb",sep="")
#5) And download to the directory
  download.file(url, destfile)
  print(paste(i,"completed..."))
}

download.file("http://www.bioinf.org.uk/teaching/bbk/biocomp2/project/data/chrom_CDS_9.gz","f:\\chme9.gz")