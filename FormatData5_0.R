
# CONCENTRATION DATA FORMATTING CODE

nf <- length(foodnames)   # Calculate the number of foods
nh <- length(hazardnames) # Calculate the number of hazards
nhK <- sum(hazardtypes=="K") # number of chemical hazards
nhM <- sum(hazardtypes=="M") # number of microbiological hazards
hazardnamesK <- hazardnames[hazardtypes=="K"] # names of chemical hazards
hazardnamesM <- hazardnames[hazardtypes=="M"] # names of microbiological hazards 

# Data frames for all exact observations and all censored observations:
obsdata <- subset(concen,!is.na(Concentration))   # Observations above LOQ
cendataLOQ <- subset(concen,(is.na(Concentration))&(!is.na(LOQ))&(!is.na(LOD))) # Censored observations below LOQ but above LOD
cendataLOD <- subset(concen,(is.na(Concentration))&(is.na(LOQ))&(!is.na(LOD)))  # Censored observations below LOD

# Calculate the number of observations above/below LOQ for each hazard:
if(nhK>0){
nexactK <- matrix(NA,nhK,nf) 
nbelowLOQK <- matrix(NA,nhK,nf) 
nbelowLODK <- matrix(NA,nhK,nf)
}
if(nhM>0){
nexactM <- matrix(NA,nhM,nf) 
nbelowLOQM <- matrix(NA,nhM,nf) 
nbelowLODM <- matrix(NA,nhM,nf)
}

if(nhK>0){
for(i in 1:nhK){
for(j in 1:nf){
nexactK[i,j] <- sum((obsdata$Hazard==hazardnamesK[i])&(obsdata$Type==foodnames[j])) # above LOQ
nbelowLOQK[i,j] <- sum((cendataLOQ$Hazard==hazardnamesK[i])&(cendataLOQ$Type==foodnames[j])) # between LOQ & LOD 
nbelowLODK[i,j] <- sum((cendataLOD$Hazard==hazardnamesK[i])&(cendataLOD$Type==foodnames[j]))
}}
}
if(nhM>0){
for(i in 1:nhM){
for(j in 1:nf){
nexactM[i,j] <- sum((obsdata$Hazard==hazardnamesM[i])&(obsdata$Type==foodnames[j])) # above LOQ
nbelowLOQM[i,j] <- sum((cendataLOQ$Hazard==hazardnamesM[i])&(cendataLOQ$Type==foodnames[j])) # between LOQ & LOD 
nbelowLODM[i,j] <- sum((cendataLOD$Hazard==hazardnamesM[i])&(cendataLOD$Type==foodnames[j]))
}}
}
frac <- 1/100 # absolute lower limit for any value is assumed as = frac*LOD 
# exact log-concentrations and log-LOQ/LOD values for each hazard, each food, each measurement  
if(nhK>0){
logcK <- array(NA,dim=c(nhK,nf,max(nexactK)))  # for exact measurements
logLOQK <- array(NA,dim=c(nhK,nf,max(nbelowLOQK)))  # for LOQ values
logLOQLimK <- array(NA,dim=c(nhK,nf,max(nbelowLOQK))) # for lower-than-LOQ limit with LOQ values
logLODK <- array(NA,dim=c(nhK,nf,max(nbelowLODK)))  # for LOD values
logLODLimK <- array(NA,dim=c(nhK,nf,max(nbelowLODK))) # for lower-than-LOD limit with LOD values
sdpriorlimK <- matrix(NA,nhK,nf)
}
if(nhM>0){
logcM <- array(NA,dim=c(nhM,nf,max(nexactM)))
logLOQM <- array(NA,dim=c(nhM,nf,max(nbelowLOQM)))
logLOQLimM <- array(NA,dim=c(nhM,nf,max(nbelowLOQM)))
logLODM <- array(NA,dim=c(nhM,nf,max(nbelowLODM)))
logLODLimM <- array(NA,dim=c(nhM,nf,max(nbelowLODM)))
sdpriorlimM <- matrix(NA,nhM,nf)
}

if(nhK>0){
for(i in 1:nhK){
for(j in 1:nf){   # natural logarithms:
if(nexactK[i,j]>0){
logcK[i,j,1:nexactK[i,j]] <- 
log(as.numeric(obsdata$Concentration[(obsdata$Hazard==hazardnamesK[i])&(obsdata$Type==foodnames[j])]))
# empirical prior upper bound for SD-parameters (if uniform prior used):
sdpriorlimK[i,j] <- sd(c(min(logcK[i,j,1:nexactK[i,j]])-1,logcK[i,j,1:nexactK[i,j]],max(logcK[i,j,1:nexactK[i,j]])+1))*3 
} 
if(nbelowLOQK[i,j]>0){
logLOQK[i,j,1:nbelowLOQK[i,j]] <- 
log(as.numeric(cendataLOQ$LOQ[(cendataLOQ$Hazard==hazardnamesK[i])&(cendataLOQ$Type==foodnames[j])]))
logLOQLimK[i,j,1:nbelowLOQK[i,j]] <- 
log(as.numeric(cendataLOQ$LOD[(cendataLOQ$Hazard==hazardnamesK[i])&(cendataLOQ$Type==foodnames[j])]))
}
if(nbelowLODK[i,j]>0){
logLODK[i,j,1:nbelowLODK[i,j]] <- 
log(as.numeric(cendataLOD$LOD[(cendataLOD$Hazard==hazardnamesK[i])&(cendataLOD$Type==foodnames[j])]))
logLODLimK[i,j,1:nbelowLODK[i,j]] <- log(frac)+logLODK[i,j,1:nbelowLODK[i,j]] 
}
}}}
if(nhM>0){
for(i in 1:nhM){
for(j in 1:nf){   # natural logarithms:
if(nexactM[i,j]>0){  
logcM[i,j,1:nexactM[i,j]] <- 
log(as.numeric(obsdata$Concentration[(obsdata$Hazard==hazardnamesM[i])&(obsdata$Type==foodnames[j])]))
# empirical prior upper bound for SD-parameters (if uniform prior used):
sdpriorlimM[i,j] <- sd(c(min(logcM[i,j,1:nexactM[i,j]])-1,logcM[i,j,1:nexactM[i,j]],max(logcM[i,j,1:nexactM[i,j]])+1))*3 
}
if(nbelowLOQM[i,j]>0){
logLOQM[i,j,1:nbelowLOQM[i,j]] <- 
log(as.numeric(cendataLOQ$LOQ[(cendataLOQ$Hazard==hazardnamesM[i])&(cendataLOQ$Type==foodnames[j])]))
logLOQLimM[i,j,1:nbelowLOQM[i,j]] <- 
log(as.numeric(cendataLOQ$LOD[(cendataLOQ$Hazard==hazardnamesM[i])&(cendataLOQ$Type==foodnames[j])]))
}
if(nbelowLODM[i,j]>0){
logLODM[i,j,1:nbelowLODM[i,j]] <- 
log(as.numeric(cendataLOD$LOD[(cendataLOD$Hazard==hazardnamesM[i])&(cendataLOD$Type==foodnames[j])]))
logLODLimM[i,j,1:nbelowLODM[i,j]] <- log(frac)+logLODM[i,j,1:nbelowLODM[i,j]] 
}
}}}


# CONSUMPTION DATA FORMATTING CODE

nd <- length(grep(foodnames[1],colnames(consum))) # number of days reported
nr <- length(IDnum)  	   # number of respondents
nf <- length(foodnames)  # number of foods

# Make order matrix that describes the position (column) of studied foods
order <- matrix(NA,nrow = nf, ncol = nd)
for(i in 1:nf){   # food
order[i,] <- grep(foodnames[i],colnames(consum))
}

# Make sure consumption information is numeric
orderv <- c(t(order)) # make order vector
consum[orderv] <- sapply(consum[orderv],as.numeric)

# Make data array for the consumption estimation (to be used in BUGS)
s <- array(NA,c(nr,nd,nf))
sw <- array(NA,c(nr,nd,nf))
Weight <- sapply(Weight,as.numeric)
for(r in 1:nr){
for(i in 1:nf){
s[r,1:nd,i] <- as.numeric(consum[r,order[i,]])  # serving size
sw[r,1:nd,i] <- as.numeric(consum[r,order[i,]]/Weight[r]) # serving size per person weight
}}
s[s == 0] <- NA
sw[is.na(s)] <- NA
logs <- log(s)    # natural logarithms of serving size
logsw <- log(sw)  # natural logarithms of serving size per weight

# To be investigated, not used currently: 
# empirically based wishart prior based on observed data correlations between foods.
# Substitute NAs with mean values, to calculate roughly food correlations from data:
#getlogsw <- matrix(NA,nr*nd,nf)
#counter <- 0
#for(r in 1:nr){
#  for(t in 1:nd){
#    counter <- counter+1
#    getlogsw[counter,1:nf] <- logsw[r,t,1:nf] 
#  }
#}
#for(i in 1:nf){
#  getlogsw[is.na(getlogsw[,i]),i] <- rep(mean(getlogsw[,i],na.rm=TRUE),sum(is.na(getlogsw[,i]))) 
#}
#DI <- solve(cov(getlogsw))  # inverse (imputed-)data covariance matrix for Wishart prior


# Consumption frequencies:
usedays <-  array(NA,dim=c(nr,nf,nd))
for(r in 1:nr){
for(i in 1:nf){
usedays[r,i,1:nd] <- as.numeric(as.numeric(consum[r,order[i,]])>0)
}}

