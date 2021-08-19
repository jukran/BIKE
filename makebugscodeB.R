
# Creates BUGS-code for a model with Markov model for consumption days ('dependent days' model)
file.create("bikemodel.txt")
fileConn<-file("bikemodel.txt")
cat("model{",file="bikemodel.txt",sep="\n")
cat("#-------------",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Concentration code",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# concentration measurements data",file="bikemodel.txt",sep="\n",append=TRUE)

# BUGS code for exact chemical measurements: 
if(nhK > 0){
  for(i in 1:nhK){
    for(j in 1:nf){
      if(OTK[i,j]=="all"){   # both true zeros and true positives possible
        cat(paste("for(k in 1:nexactK[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("logcK[",i,",",j,",k] ~ dnorm(mucK[",i,",",j,"],taucK[",i,",",j,"])"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("positive1K[",i,",",j,",k] ~ dbern(pK[",i,",",j,"])"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("positive1K[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat("}",file="bikemodel.txt",sep="\n",append=TRUE)   
      }
      if(OTK[i,j]=="positives"){ # only positives
        cat(paste("for(k in 1:nexactK[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("logcK[",i,",",j,",k] ~ dnorm(mucK[",i,",",j,"],taucK[",i,",",j,"])"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
      }
}}}

# BUGS code for exact microbiological measurements:
if(nhM > 0){
  for(i in 1:nhM){
    for(j in 1:nf){ 
      if(OTM[i,j]=="all"){
        cat(paste("for(k in 1:nexactM[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("logcM[",i,",",j,",k] ~ dnorm(mucM[",i,",",j,"],taucM[",i,",",j,"])"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("positive1M[",i,",",j,",k] ~ dbern(pM[",i,",",j,"])"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("positive1M[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
      }
      if(OTM[i,j]=="positives"){
        cat(paste("for(k in 1:nexactM[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
        cat(paste("logcM[",i,",",j,",k] ~ dnorm(mucM[",i,",",j,"],taucM[",i,",",j,"])"),file="bikemodel.txt",sep="\n",append=TRUE) 
        cat("}",file="bikemodel.txt",sep="\n",append=TRUE) 
      }
}}}

cat("# Censored data",file="bikemodel.txt",sep="\n",append=TRUE)

# BUGS code for censored chemical measurements between LOQ and LOD
if(nhK > 0){
  for(i in 1:nhK){
    for(j in 1:nf){
      if(nbelowLOQK[i,j]>0){
        if(OTK[i,j]=="all"){  
          cat(paste("for(k in 1:nbelowLOQK[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2K[",i,",",j,",k] ~ dbern(PRintervalK[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRintervalK[",i,",",j,",k]<-pK[",i,",",j,"]*(phi((logLOQK[",i,",",j,",k]-mucK[",i,",",j,"])*pow(taucK[",i,",",j,"],0.5))-phi((logLOQLimK[",i,",",j,",k]-mucK[",i,",",j,"])*pow(taucK[",i,",",j,"],0.5)))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2K[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE) 
        }
        if(OTK[i,j]=="positives"){
          cat(paste("for(k in 1:nbelowLOQK[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2K[",i,",",j,",k] ~ dbern(PRintervalK[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRintervalK[",i,",",j,",k]<-(phi((logLOQK[",i,",",j,",k]-mucK[",i,",",j,"])*pow(taucK[",i,",",j,"],0.5))-phi((logLOQLimK[",i,",",j,",k]-mucK[",i,",",j,"])*pow(taucK[",i,",",j,"],0.5)))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2K[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE)    
        }
      }
}}}

# BUGS code for censored chemical measurements below LOD
if(nhK>0){
  for(i in 1:nhK){
    for(j in 1:nf){
      if(nbelowLODK[i,j]>0){
        if(OTK[i,j]=="all"){   
          cat(paste("for(k in 1:nbelowLODK[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3K[",i,",",j,",k] ~ dbern(PRleftK[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRleftK[",i,",",j,",k]<- 1-pK[",i,",",j,"]*(1-phi((logLODK[",i,",",j,",k]-mucK[",i,",",j,"])*pow(taucK[",i,",",j,"],0.5)))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3K[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
        }
        if(OTK[i,j]=="positives"){
          cat(paste("for(k in 1:nbelowLODK[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3K[",i,",",j,",k] ~ dbern(PRleftK[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRleftK[",i,",",j,",k]<- phi((logLODK[",i,",",j,",k]-mucK[",i,",",j,"])*pow(taucK[",i,",",j,"],0.5))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3K[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE)   
        }
      }
}}}   


# BUGS code for censored microbiological measurements between LOD & LOQ and possibly including true zeros
if(nhM > 0){
  for(i in 1:nhM){
    for(j in 1:nf){
      if(nbelowLOQM[i,j]>0){
        if(OTM[i,j]=="all"){  
          cat(paste("for(k in 1:nbelowLOQM[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2M[",i,",",j,",k] ~ dbern(PRintervalM[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRintervalM[",i,",",j,",k]<-pM[",i,",",j,"]*(phi((logLOQM[",i,",",j,",k]-mucM[",i,",",j,"])*pow(taucM[",i,",",j,"],0.5))-phi((logLOQLimM[",i,",",j,",k]-mucM[",i,",",j,"])*pow(taucM[",i,",",j,"],0.5)))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2M[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE) 
        }
        if(OTM[i,j]=="positives"){
          cat(paste("for(k in 1:nbelowLOQM[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2M[",i,",",j,",k] ~ dbern(PRintervalM[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRintervalM[",i,",",j,",k]<-(phi((logLOQM[",i,",",j,",k]-mucM[",i,",",j,"])*pow(taucM[",i,",",j,"],0.5))-phi((logLOQLimM[",i,",",j,",k]-mucM[",i,",",j,"])*pow(taucM[",i,",",j,"],0.5)))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive2M[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE)    
        }
      } # end of if positives  
}}}

# BUGS code for censored microbiological measurements below LOD and possibly including true zeros
if(nhM > 0){
  for(i in 1:nhM){
    for(j in 1:nf){
      if(nbelowLODM[i,j]>0){
        if(OTM[i,j]=="all"){   
          cat(paste("for(k in 1:nbelowLODM[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3M[",i,",",j,",k] ~ dbern(PRleftM[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRleftM[",i,",",j,",k]<- 1-pM[",i,",",j,"]*(1-phi((logLODM[",i,",",j,",k]-mucM[",i,",",j,"])*pow(taucM[",i,",",j,"],0.5)))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3M[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
        }
        if(OTM[i,j]=="positives"){
          cat(paste("for(k in 1:nbelowLODM[",i,",",j,"]){"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3M[",i,",",j,",k] ~ dbern(PRleftM[",i,",",j,",k])"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("PRleftM[",i,",",j,",k]<- phi((logLODM[",i,",",j,",k]-mucM[",i,",",j,"])*pow(taucM[",i,",",j,"],0.5))"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat(paste("positive3M[",i,",",j,",k] <- 1"),file="bikemodel.txt",sep="\n",append=TRUE)
          cat("}",file="bikemodel.txt",sep="\n",append=TRUE)   
        }
      } # end of if positives  
}}}

if(nhK>0){
  cat("# Priors:",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("for(i in 1:nhK){",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("for(j in 1:nf){",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("mucK[i,j] ~ dunif(-10,10)",file="bikemodel.txt",sep="\n",append=TRUE)
  if(input$priorchoice=="sigma_uniform"){
  cat("taucK[i,j] <- pow(sigcK[i,j],-2)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("sigcK[i,j] ~ dunif(0,sdpriorlimK[i,j])",file="bikemodel.txt",sep="\n",append=TRUE)
  }
  if(input$priorchoice=="tau_gamma"){
  cat("taucK[i,j] ~ dgamma(0.01,0.01)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("sigcK[i,j] <- sqrt(1/taucK[i,j])",file="bikemodel.txt",sep="\n",append=TRUE)  
  }   
  cat("}}",file="bikemodel.txt",sep="\n",append=TRUE)
  for(i in 1:nhK){
    for(j in 1:nf){
      if(OTK[i,j]=="all"){
        cat(paste("pK[",i,",",j,"] ~ dbeta(1,1)"),file="bikemodel.txt",sep="\n",append=TRUE)
      }
      if(OTK[i,j]=="positives"){
        cat(paste("pK[",i,",",j,"] ~ dbeta(",nposK[i,j]+1,",",nsamK[i,j]-nposK[i,j]+1,")"),file="bikemodel.txt",sep="\n",append=TRUE)
      }  
    }}
}
if(nhM>0){
  cat("# Priors:",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("for(i in 1:nhM){",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("for(j in 1:nf){",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("mucM[i,j] ~ dunif(-10,10)",file="bikemodel.txt",sep="\n",append=TRUE)
if(input$priorchoice=="sigma_uniform"){  
  cat("taucM[i,j] <- pow(sigcM[i,j],-2)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("sigcM[i,j] ~ dunif(0,sdpriorlimM[i,j])",file="bikemodel.txt",sep="\n",append=TRUE)
}
if(input$priorchoice=="tau_gamma"){  
  cat("taucM[i,j] ~ dgamma(0.01,0.01)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("sigcM[i,j] <- sqrt(1/taucM[i,j])",file="bikemodel.txt",sep="\n",append=TRUE)
}  
  cat("}}",file="bikemodel.txt",sep="\n",append=TRUE)
  for(i in 1:nhM){
    for(j in 1:nf){
      if(OTM[i,j]=="all"){
        cat(paste("pM[",i,",",j,"] ~ dbeta(1,1)"),file="bikemodel.txt",sep="\n",append=TRUE)
      }
      if(OTM[i,j]=="positives"){
        cat(paste("pM[",i,",",j,"] ~ dbeta(",nposM[i,j]+1,",",nsamM[i,j]-nposM[i,j]+1,")"),file="bikemodel.txt",sep="\n",append=TRUE)
      }  
    }}
}

cat("#------------------------",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Consumption code",file="bikemodel.txt",sep="\n",append=TRUE)

cat("# Body weight model:",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(r in 1:nr){ logWeight[r] ~ dnorm(muw,tauw); logWeight[r] <- log(Weight[r]) }",file="bikemodel.txt",sep="\n",append=TRUE)
cat("muw ~ dunif(-10,10)",file="bikemodel.txt",sep="\n",append=TRUE)
if(input$priorchoice=="sigma_uniform"){
  cat("tauw <- pow(sigw,-2)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("sigw ~ dunif(0,10)",file="bikemodel.txt",sep="\n",append=TRUE)
} 
if(input$priorchoice=="tau_gamma"){
  cat("tauw ~ dgamma(0.01,0.01)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("sigw <- sqrt(1/tauw)",file="bikemodel.txt",sep="\n",append=TRUE)
} 
cat("  ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("  ",file="bikemodel.txt",sep="\n",append=TRUE)


cat("# Consumption measurements data (daily amounts)",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(r in 1:nr){ # individual respondent",file="bikemodel.txt",sep="\n",append=TRUE)
## cat("for(j in 1:nf){ # food",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(t in 1:nd){ # day",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# non-consumption days are here as 'NA' so that this",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# distribution for what the consumption would be if it was positive:",file="bikemodel.txt",sep="\n",append=TRUE)

if(nf>1){
  cat("logsw[r,t,1:nf] ~ dmnorm(mus[r,1:nf],Ts[1:nf,1:nf])",file="bikemodel.txt",sep="\n",append=TRUE)
}
if(nf==1){
  cat("logsw[r,t,1] ~ dnorm(mus[r,1],Ts[1,1])",file="bikemodel.txt",sep="\n",append=TRUE)  
}
cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
if(nf>1){
  cat("# Individual consumption means are correlated among food types",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("# with a population (overall) mean and correlation matrix over food types:",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("mus[r,1:nf] ~ dmnorm(mus0[1:nf],Ts0[1:nf,1:nf])",file="bikemodel.txt",sep="\n",append=TRUE)
}
if(nf==1){
  cat("mus[r,1] ~ dnorm(mus0[1],Ts0[1,1])",file="bikemodel.txt",sep="\n",append=TRUE)  
}
cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Priors for food consumption amounts in general:",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(j in 1:nf){ ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("mus0[j] ~ dunif(-10,10)",file="bikemodel.txt",sep="\n",append=TRUE)
cat("sigs[j] <- sqrt(Ss[j,j])",file="bikemodel.txt",sep="\n",append=TRUE)
cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
if(nf>1){
cat("# Correlations between food type means:",file="bikemodel.txt",sep="\n",append=TRUE)
cat("Ts0[1:nf,1:nf] ~ dwish(DI[1:nf,1:nf],nf2); nf2<-nf+2",file="bikemodel.txt",sep="\n",append=TRUE)
cat("Ts[1:nf,1:nf] ~ dwish(DI[1:nf,1:nf],nf2)",file="bikemodel.txt",sep="\n",append=TRUE)
cat("Ss0[1:nf,1:nf] <- inverse(Ts0[1:nf,1:nf])",file="bikemodel.txt",sep="\n",append=TRUE)
cat("Ss[1:nf,1:nf] <- inverse(Ts[1:nf,1:nf])",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Optional outputs: evaluate correlation matrix between food type amounts: ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("#for(ii in 1:nf){  ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("#for(jj in 1:nf){ ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("#CORlogs0[ii,jj] <- (Ss0[ii,jj])/((sqrt(Ss0[ii,ii]))*(sqrt(Ss0[jj,jj])))",file="bikemodel.txt",sep="\n",append=TRUE)
cat("#CORlogs[ii,jj] <- (Ss[ii,jj])/((sqrt(Ss[ii,ii]))*(sqrt(Ss[jj,jj])))",file="bikemodel.txt",sep="\n",append=TRUE)
cat("#}}",file="bikemodel.txt",sep="\n",append=TRUE)
}
if(nf==1){
  cat("Ts0[1,1] ~ dgamma(1,0.1)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("Ts[1,1] ~ dgamma(1,0.1)",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("Ss0[1,1] <- 1/Ts0[1,1]",file="bikemodel.txt",sep="\n",append=TRUE)
  cat("Ss[1,1] <- 1/Ts[1,1]",file="bikemodel.txt",sep="\n",append=TRUE)
}
cat("# Consumption frequencies data (daily yes/no)",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(r in 1:nr){",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Individual consumption probabilities dependent over days:",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(j in 1:nf){",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(t in 2:nd){",file="bikemodel.txt",sep="\n",append=TRUE)
cat("usedays[r,j,t] ~ dbern(pmarkov[r,j,t]) ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Probability to consume (individual r, food j):",file="bikemodel.txt",sep="\n",append=TRUE)
cat("pmarkov[r,j,t] <- p11[j]*usedays[r,j,t-1]+p01[j]*(1-usedays[r,j,t-1]) ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("}}}",file="bikemodel.txt",sep="\n",append=TRUE)

cat("# Prior for consumption probabilities:",file="bikemodel.txt",sep="\n",append=TRUE)
cat("for(j in 1:nf){",file="bikemodel.txt",sep="\n",append=TRUE)
cat("p01[j] ~ dunif(0,1)",file="bikemodel.txt",sep="\n",append=TRUE)
cat("p11[j] ~ dunif(0,1)",file="bikemodel.txt",sep="\n",append=TRUE)
cat("p10[j] <- 1-p11[j] ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("# Long run consumption probability:",file="bikemodel.txt",sep="\n",append=TRUE)
cat("ppred[j] <- p01[j]/(p01[j]+p10[j]) ",file="bikemodel.txt",sep="\n",append=TRUE)
cat("}",file="bikemodel.txt",sep="\n",append=TRUE)

if(nf>1){
cat("for(i in 1:nf){for(j in 1:nf){ DI[i,j]<-equals(i,j) }}",file="bikemodel.txt",sep="\n",append=TRUE) 
}


cat("}",file="bikemodel.txt",sep="\n",append=TRUE)
close(fileConn)


