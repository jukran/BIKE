
results <- function(iterations){
  # Run the BUGS model code (MCMC sampling of posterior distribution) for given number of iterations.
  # For this: data and initial values for parameters need to be defined,
  # and the names of parameters to be monitored in outputs.
  # --> definitions for 'data', 'inits', 'parameters'.
  # --> call bugs.
  # --> return MCMC outputs.
  
  
  # if both chemical and microbiological hazards:    
  if((nhK>0)&(nhM>0)){
    data <- list("Weight","logcM","logcK","sdpriorlimM","sdpriorlimK","nexactM","nexactK",
                 "nhM","nhK","nf","nr","nd","logsw","usedays")  
    if(sum(nbelowLOQM)>0){ 
      data <- append(data,c("logLOQM","logLOQLimM","nbelowLOQM"))
    }
    if(sum(nbelowLODM)>0){ 
      data <- append(data,c("logLODM","nbelowLODM"))
    }
    if(sum(nbelowLOQK)>0){ 
      data <- append(data,c("logLOQK","logLOQLimK","nbelowLOQK"))
    }
    if(sum(nbelowLODK)>0){ 
      data <- append(data,c("logLODK","nbelowLODK"))
    }
    # initial values:
    initmucK <- matrix(NA,nhK,nf)
    initsigcK <- matrix(NA,nhK,nf)
    for(i in 1:nhK){
      for(j in 1:nf){
        initmucK[i,j] <- mean(logcK[i,j,],na.rm=TRUE)
        initsigcK[i,j] <- sdpriorlimK[i,j]/2
      }
    }
    initmucM <- matrix(NA,nhM,nf)
    initsigcM <- matrix(NA,nhM,nf)
    for(i in 1:nhM){
      for(j in 1:nf){
        initmucM[i,j] <- mean(logcM[i,j,],na.rm=TRUE)
        initsigcM[i,j] <- sdpriorlimM[i,j]/2
      }
    }
  } # end of both
  
  # if only chemical hazards:  
  if((nhK>0)&(nhM==0)){
    data <- list("Weight","logcK","sdpriorlimK","nexactK",
                 "nhK","nf","nr","nd","logsw","usedays")  
    if(sum(nbelowLOQK)>0){ 
      data <- append(data,c("logLOQK","logLOQLimK","nbelowLOQK"))
    }
    if(sum(nbelowLODK)>0){ 
      data <- append(data,c("logLODK","nbelowLODK"))
    }
    # initial values:
    initmucK <- matrix(NA,nhK,nf)
    initsigcK <- matrix(NA,nhK,nf)
    for(i in 1:nhK){
      for(j in 1:nf){
        initmucK[i,j] <- mean(logcK[i,j,],na.rm=TRUE)
        initsigcK[i,j] <- sdpriorlimK[i,j]/2
      }
    }
  } # end of only chemical 
  
  # if only microbiological hazards:  
  if((nhK==0)&(nhM>0)){
    data <- list("Weight","logcM","sdpriorlimM","nexactM",
                 "nhM","nf","nr","nd","logsw","usedays") 
    if(sum(nbelowLOQM)>0){ 
      data <- append(data,c("logLOQM","logLOQLimM","nbelowLOQM"))
    }
    if(sum(nbelowLODM)>0){ 
      data <- append(data,c("logLODM","nbelowLODM"))
    }
    # initial values:
    initmucM <- matrix(NA,nhM,nf)
    initsigcM <- matrix(NA,nhM,nf)
    for(i in 1:nhM){
      for(j in 1:nf){
        initmucM[i,j] <- mean(logcM[i,j,],na.rm=TRUE)
        initsigcM[i,j] <- sdpriorlimM[i,j]/2
      }
    }
  } # end of only microbiological
  
  
  if(input$priorchoice == "sigma_uniform"){ # then inits for sigma are needed
  if((nhK>0)&(nhM>0)){
    inits <- function(){list(muw=0,sigw=1,mucK=initmucK,sigcK=initsigcK,mucM=initmucM,sigcM=initsigcM,pM=matrix(0.5,nhM,nf),mus=matrix(0,nr,nf),mus0=rep(0,nf)) }
    parameters=c("mucK","mucM","sigcK","sigcM","mus0","Ts0","Ts","sigs","Ss0","Ss","ppred","muw","sigw","pM","pK")
  }
  if((nhK>0)&(nhM==0)){
    inits <- function(){list(muw=0,sigw=1,mucK=initmucK,sigcK=initsigcK,mus=matrix(0,nr,nf),mus0=rep(0,nf)) }  
    parameters=c("mucK","sigcK","mus0","Ts0","Ts","sigs","Ss0","Ss","ppred","muw","sigw","pK")
  }
  if((nhK==0)&(nhM>0)){
    inits <- function(){list(muw=0,sigw=1,mucM=initmucM,sigcM=initsigcM,pM=matrix(0.5,nhM,nf),mus=matrix(0,nr,nf),mus0=rep(0,nf)) }  
    parameters=c("mucM","sigcM","mus0","Ts0","Ts","sigs","Ss0","Ss","ppred","muw","sigw","pM")
  }  
  }  # end of prior choice
  if(input$priorchoice == "tau_gamma"){ # then inits for tau are needed
    if((nhK>0)&(nhM>0)){
      inits <- function(){list(muw=0,tauw=1,mucK=initmucK,taucK=initsigcK^(-2),mucM=initmucM,sigcM=initsigcM,pM=matrix(0.5,nhM,nf),mus=matrix(0,nr,nf),mus0=rep(0,nf)) }
      parameters=c("mucK","mucM","sigcK","sigcM","mus0","Ts0","Ts","sigs","Ss0","Ss","ppred","muw","sigw","pM","pK")
    }
    if((nhK>0)&(nhM==0)){
      inits <- function(){list(muw=0,tauw=1,mucK=initmucK,taucK=initsigcK^(-2),mus=matrix(0,nr,nf),mus0=rep(0,nf)) }  
      parameters=c("mucK","sigcK","mus0","Ts0","Ts","sigs","Ss0","Ss","ppred","muw","sigw","pK")
    }
    if((nhK==0)&(nhM>0)){
      inits <- function(){list(muw=0,tauw=1,mucM=initmucM,taucM=initsigcM^(-2),pM=matrix(0.5,nhM,nf),mus=matrix(0,nr,nf),mus0=rep(0,nf)) }  
      parameters=c("mucM","sigcM","mus0","Ts0","Ts","sigs","Ss0","Ss","ppred","muw","sigw","pM")
    }  
  }  # end of prior choice
  
# if prevalence sample data is given, include it:   
if(is.element("positives",OTK)){
  data <- append(data,c("nposK","nsamK"))    
}  
if(is.element("positives",OTM)){
  data <- append(data,c("nposM","nsamM"))    
}
    
dasim = bugs(data,inits,model.file="bikemodel.txt",debug=FALSE,parameters,n.chains=1,n.burnin=burnin,n.iter=iterations,DIC=FALSE,codaPkg=FALSE)
attach.bugs(dasim)

if((nhK>0)&(nhM>0)){
  outputs <- data.frame(mucK,mucM,pM,pK,sigcK,sigcM,sigs,Ts0,Ss0,Ts,Ss,ppred,muw,sigw)
}
if((nhK==0)&(nhM>0)){
  outputs <- data.frame(mucM,sigcM,pM,sigs,Ts0,Ss0,Ts,Ss,ppred,muw,sigw)
}
if((nhK>0)&(nhM==0)){
  outputs <- data.frame(mucK,sigcK,pK,sigs,Ts0,Ss0,Ts,Ss,ppred,muw,sigw)
}

return(outputs)
}


