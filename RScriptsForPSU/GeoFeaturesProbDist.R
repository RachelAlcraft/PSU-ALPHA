
setwd("F:\\PSUA\\GeoDatabase\\")
data <- scan("ALA_geoprobdist.csv",character(),quote="")

row <- data[6]
length(row[1])
vecrow <- unlist(strsplit(row,","))
header <- vecrow[1]
datarow <- as.numeric(vecrow[2:2000])
  
#xkey = "Distance (A)"
xkey = "Degrees"
ykey <- "Probability distribution from 1000 pdbs"


xkey <- paste(xkey,"from", length(datarow),"values")

hist(datarow,main=header,xlab=xkey,ylab=ykey, xlim=c(114,130))#,breaks=7)