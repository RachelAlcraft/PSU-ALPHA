



##################################################################################

##################################################################################

#1) Set the working directory to collect the file with the pdb codes in
#workingdirectory = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-HighResList\\"

# workingdirectory = "F:\\PSUA\\Code\\PSU-ALPHA\\MSC-RemoveSimilarity\\"
workingdirectory = "F:\\Code\\BbkProject\\Thesis\\Method\\02_DataSets\\nonsim_lists\\"
setwd(workingdirectory)

#2) Get the file containing the pdb codes you want to download
#filename="highres_and_unique.txt"

filename="extra_nonsim.csv"
#filename="2018b_nonsim.csv"
#filename="high3_nonsim.csv"

pdbcodes = read.csv(filename, header = TRUE)

#3) Now set the workling directory again to the destination database directory
#workingdirectory = "F:\\PSUA\\ProjectData\\HighResFiles\\"

workingdirectory = "F:\\Code\\BbkTransfer\\pdbfiles\\pdbdata\\"
setwd(workingdirectory)

#4) Loop through the codes to get each one

for (i in 1:length(pdbcodes$PDB)){

  pdb = paste(pdbcodes$PDB[i])
  url = paste("https://files.rcsb.org/download/",pdb,".pdb",sep="")
  destfile = paste(pdb,".pdb",sep="")
  #5) And download to the directory
  if (!file.exists(destfile)){
    try(download.file(url, destfile))
    print(paste(i,"/",length(pdbcodes$PDB),"completed pdb..."))
  }
  else{
    print(paste(i,"/",length(pdbcodes$PDB),"exists pdb..."))  
  }
  
  #6) Also download structure factors if they exist
  cif = paste(pdbcodes$PDB[i],"-sf",sep="") #structure factors
  urlcif = paste("https://files.rcsb.org/download/",cif,".cif",sep="")
  destfilecif = paste(cif,".cif",sep="")
  if (!file.exists(destfilecif)){
    try(download.file(urlcif, destfilecif))
    print(paste(i,"/",length(pdbcodes$PDB),"completed cif..."))
  }
  else{
    print(paste(i,"/",length(pdbcodes$PDB),"exists cif..."))  
  }
  
}

