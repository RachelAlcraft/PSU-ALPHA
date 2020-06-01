



##################################################################################

##################################################################################

#1) Set the working directory to collect the file with the pdb codes in
#workingdirectory = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-HighResList\\"

workingdirectory = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-RemoveSimilarity\\"
setwd(workingdirectory)

#2) Get the file containing the pdb codes you want to download
#filename="highres_and_unique.txt"

filename="2019nonsim95.csv"
pdbcodes = read.csv(filename, header = TRUE)

#3) Now set the workling directory again to the destination database directory
#workingdirectory = "F:\\PSUA\\ProjectData\\HighResFiles\\"

workingdirectory = "F:2\PSUA\22019Data\\pdbs\\"
setwd(workingdirectory)

#4) Loop through the codes to get each one

for (i in 1:length(pdbcodes$PDB)){

  pdb = paste(pdbcodes$PDB[i])
  url = paste("https://files.rcsb.org/download/",pdb,".pdb",sep="")
  destfile = paste(pdb,".pdb",sep="")
  #5) And download to the directory
  try(download.file(url, destfile))
  
  #6) Also download structure factors if they exist
  cif = paste(pdbcodes$PDB[i],"-sf",sep="") #structure factors
  urlcif = paste("https://files.rcsb.org/download/",cif,".cif",sep="")
  destfilecif = paste(cif,".cif",sep="")
  try(download.file(urlcif, destfilecif))
  
  print(paste(i,"completed..."))
}

