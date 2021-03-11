
# Read names of foods and hazards 
ocdata <- read.table("DataOccurrence.txt",header=TRUE,dec=",",stringsAsFactors = FALSE)
attach(ocdata)
foodnames<-names(ocdata[3:dim(ocdata)[2]])
hazardnamesM<-hazardnames[hazardtypes=="M"]
hazardnamesK<-hazardnames[hazardtypes=="K"]
if(sum(hazardtypes=="M")>0){
  # Occurrence Table of Microbials consists of marks "all" or "positives" to denote data type:
OTM<-matrix(NA,sum(hazardtypes=="M"),length(foodnames))  
OTM[1:sum(hazardtypes=="M"),1:length(foodnames)]<-as.matrix(ocdata[hazardtypes=="M",3:length(ocdata[1,])])  
} 
if(sum(hazardtypes=="K")>0){
  # Occurrence Table of Chemicals consists of marks "all" or "positives" to denote data type:
  OTK<-matrix(NA,sum(hazardtypes=="K"),length(foodnames))  
  OTK[1:sum(hazardtypes=="K"),1:length(foodnames)]<-as.matrix(ocdata[hazardtypes=="K",3:length(ocdata[1,])])  
}

# prevalence data for microbials: 
prevdata <- read.table("DataPrevalence.txt",header=TRUE,dec=",",stringsAsFactors = FALSE)
nposM <- matrix(NA,sum(hazardtypes=="M"),length(foodnames))
nsamM <- matrix(NA,sum(hazardtypes=="M"),length(foodnames))
for(i in 1:dim(prevdata)[1]){
  nposM[which(prevdata[i,1]==hazardnamesM),which(prevdata[i,3]==foodnames)] <- prevdata[i,4]
  nsamM[which(prevdata[i,1]==hazardnamesM),which(prevdata[i,3]==foodnames)] <- prevdata[i,5]
}
# prevalence data for chemicals:
nposK <- matrix(NA,sum(hazardtypes=="K"),length(foodnames))
nsamK <- matrix(NA,sum(hazardtypes=="K"),length(foodnames))
for(i in 1:dim(prevdata)[1]){
  nposK[which(prevdata[i,1]==hazardnamesK),which(prevdata[i,3]==foodnames)] <- prevdata[i,4]
  nsamK[which(prevdata[i,1]==hazardnamesK),which(prevdata[i,3]==foodnames)] <- prevdata[i,5]
}

# Read concentration data (original Excel-file saved as '.txt' file)
concen <- read.table("DataConcentrations.txt",header=TRUE,dec=",",stringsAsFactors = FALSE)
attach(concen)
# Read consumption data (original Excel-file saved as '.txt' file)
consum <- read.table("DataConsumptions.txt",header=TRUE,dec=",",stringsAsFactors = FALSE)
attach(consum)

