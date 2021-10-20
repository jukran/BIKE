#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# packages needed:
library(shiny)
library(R2OpenBUGS)
library(mvtnorm)
library(shinyMatrix)

# read data files:
source("ReadData.R")
raspberry <- rgb(208/255,0,111/255) # color definition
burnin <- 1000  # number of burnin iterations for MCMC runs

# Define UI for application that runs BUGS model and draws results
ui <- fluidPage(
  
   # Application title 
   titlePanel("BIKE foodborne exposure model"),
   
   sidebarLayout(
      sidebarPanel(
        
        checkboxGroupInput("thefoodnames","Food types to select:",
                           choices=foodnames,selected=foodnames[1]  ),
        checkboxGroupInput("thehazardnames","Hazards to select:",
                           choices=hazardnames,selected=hazardnames[1]  ),
        checkboxGroupInput("selectresults","Results to view:",
                           choices=c("Concentrations","Consumptions","Exposures","Posterior predictive","Serving correlations","Mean serving correlations"),selected="Concentrations"),
        radioButtons("selectQ",label="Estimated total exposure quantile:",choices=c("None","Q5% Exposure","Q10% Exposure","Q25% Exposure","Q50% Exposure","Q75% Exposure","Q90% Exposure","Q95% Exposure")),
        radioButtons("selectscale",label="Scale:",choices=c("Absolute","Logarithmic")),
        radioButtons("selectdist",label="Distributions:",choices=c("Cumulative","Density")),
        radioButtons("modelchoice",label="Consumption model:",choices=c("Independent days","Dependent days")),
        radioButtons("modelchoice2",label="Between-user variability in consumption frequencies (for 'independent days' model only):",choices=c("No","Yes")),
        radioButtons("priorchoice",label="Priors for variances:",choices=c("tau_gamma","sigma_uniform")),
        radioButtons("diagchoice",label="Visual MCMC sample:",choices=c("None","Concentration parameters","Consumption parameters")),
        sliderInput("Iterations",
                     "Number of MCMC iterations:",
                     min = 2000,
                     max = 150000,
                     value = 4000),
        sliderInput("nV","Variability sample size for Q%",
                    min=50,
                    max=1000,
                    value=100),
        sliderInput("nU","Uncertainty sample size for Q%",
                    min=50,
                    max=1000,
                    value=100),
      titlePanel("Processing factors:"),
                   matrixInput(
                     "factor",class="numeric",
                     value = matrix(1,length(foodnames),length(hazardnames),dimnames=list(substr(foodnames,1,3),substr(hazardnames,1,3))),
                     rows = list(names = TRUE),
                     cols = list(names = TRUE),
                     copy = TRUE,
                     paste = TRUE
                   ),
      titlePanel("Prevalence factors:"),
                   matrixInput(
                   "pfactor",class="numeric",
                   value = matrix(1,length(foodnames),length(hazardnames),dimnames=list(substr(foodnames,1,3),substr(hazardnames,1,3))),
                   rows = list(names = TRUE),
                   cols = list(names = TRUE),
                   copy = TRUE,
                   paste = TRUE
                   )
      ),
      
      # Show a plot of the generated results
      mainPanel(
         plotOutput("distPlot1"),    # concentrations plot
         plotOutput("distPlot2"),    # consumptions plot
         plotOutput("distPlot3"),    # exposures plot
         plotOutput("distPlot4"),    # quantile plot  
         plotOutput("distPlot5"),    # diagnostics plot
         plotOutput("distPlot6"),    # serving correlations plot
         plotOutput("distPlot7"),    # mean serving correlations plot
         textOutput("header"),
         tableOutput("values")       # posterior predictive quantiles table
      )
   )
)

# Define server logic required to run BUGS model and draw results
server <- function(input, output) {
  
  source("FormatData.R",local=TRUE) # read and format data
  
  # compute results from the full model once, then just use as 'current' results for all plots
  # when only post-processing of the MCMC output is needed.
  currentresults  <- reactive({ 
    
    # model assumes independent daily consumptions
    if(input$modelchoice == "Independent days"){  
      # model with user variability in consumption frequency
      if(input$modelchoice2 == "No"){ between.user.pvar <- 0 } 
      # model without user variability in consumption frequency
      if(input$modelchoice2 == "Yes"){ between.user.pvar <- 1 }
      source("makebugscodeA.R",local=TRUE)     # write code for OpenBUGS
      source("RunBUGSfunctionA.R",local=TRUE)  # run BUGS
    } 
    # model assumes that consuming next day depends on consuming previous day 
    if(input$modelchoice == "Dependent days"){  
      source("makebugscodeB.R",local=TRUE)     # write code for OpenBUGS
      source("RunBUGSfunctionB.R",local=TRUE)  # run BUGS
    }
    
    results(round(as.numeric(input$Iterations))) 
    }) 
  
####################################################  
   output$distPlot1 <- renderPlot({
     # generate results based on inputs from ui.R:  Concentrations
     
     currentresults()  
     
     theresults = input$selectresults  
     
     if(is.element("Concentrations",theresults)){
      
     scale = input$selectscale  # absolute or logarithmic  
     foodnamesused = input$thefoodnames # selected foods
     nfused = length(foodnamesused)     # number of selected foods
     foodindex = match(foodnamesused,foodnames) # indexing of selected foods in all foods
     
     hazardnamesused = input$thehazardnames # selected hazards
     hazardtypesused = hazardtypes[is.element(hazardnames,hazardnamesused)] # types of selected hazards (chemical/microbiological)
     nhused = length(hazardnamesused) # number of selected hazards
     hazardnamesK = hazardnames[hazardtypes=="K"] # chemical hazard names
     hazardnamesM = hazardnames[hazardtypes=="M"] # microbiological hazard names
     hazardnamesusedK = hazardnamesused[hazardtypesused=="K"] # selected che hazard names
     hazardnamesusedM = hazardnamesused[hazardtypesused=="M"] # selected mic hazard names
     nhusedK = length(hazardnamesusedK) # number of che hazards selected
     nhusedM = length(hazardnamesusedM) # number of mic hazards selected
     hazardindex = match(hazardnamesused,hazardnames) # 
     hazardindexK = match(hazardnamesusedK,hazardnamesK) # indexing of selected hazards in all che hazards
     hazardindexM = match(hazardnamesusedM,hazardnamesM) # indexing of selected hazards in all mic hazards 
  
     if( (nhused>0)&(nfused>0) ){par(mfrow=c(nhused,nfused),cex.lab=1.3,cex.main=1.3,yaxt="n")}
     
     # Chemical concentrations:
    
     if((nhusedK>0)&(nfused>0)){
       # redefine dimensions if scalars:
       if((nhK==1)&(nf==1)){
         mucK <- array(mucK,dim=c(input$Iterations-burnin,1,1))
         sigcK <- array(sigcK,dim=c(input$Iterations-burnin,1,1))
         pK <- array(pK,dim=c(input$Iterations-burnin,1,1))
       }    
     for(h in 1:nhusedK){
     for(i in 1:nfused){
        if(input$selectdist=="Density"){
        if(scale=="Absolute"){
        cmeanK <- exp(mucK[,hazardindexK[h],foodindex[i]]+0.5*sigcK[,hazardindexK[h],foodindex[i]]^2)
        cmedianK <- exp(mucK[,hazardindexK[h],foodindex[i]])
        maxx = exp(mean(mucK[,hazardindexK[h],foodindex[i]]+2*sigcK[,hazardindexK[h],foodindex[i]]))
        plot(density(cmeanK[cmeanK<maxx]),col="darkgray",lwd=3,main=paste(hazardnamesusedK[h],"in",foodnamesused[i]),xlab="Concentration+",ylab="",xlim=c(0,maxx)) 
        polygon(density(cmeanK[cmeanK<maxx]),col="lightgray")
        lines(density(cmedianK[cmedianK<maxx]),lwd=3)
        # mark data points and possible LOD and LOQ values for censored data:
        rug(exp(logcK[hazardindexK[h],foodindex[i],]),lwd=2.5,col="red",quiet=TRUE)
        rug(exp(logLOQK[hazardindexK[h],foodindex[i],]),lwd=4.5,col="green",quiet=TRUE)
        rug(exp(logLODK[hazardindexK[h],foodindex[i],]),lwd=4.5,col="blue",quiet=TRUE)
        }
        if(scale=="Logarithmic"){
         maxx <- mean(mucK[,hazardindexK[h],foodindex[i]]+2*sigcK[,hazardindexK[h],foodindex[i]])
         minn <- mean(mucK[,hazardindexK[h],foodindex[i]]-2*sigcK[,hazardindexK[h],foodindex[i]])
        plot(density(mucK[,hazardindexK[h],foodindex[i]]/log(10)),col="darkgray",lwd=3,main=paste(hazardnamesusedK[h],"in",foodnamesused[i]),xlab="log Concentration+",ylab="",xlim=c(minn/log(10),maxx/log(10))) 
        polygon(density(mucK[,hazardindexK[h],foodindex[i]]/log(10)),col="lightgray")
        # mark data points and possible LOD and LOQ values for censored data:
        rug(logcK[hazardindexK[h],foodindex[i],]/log(10),lwd=2.5,col="red",quiet=TRUE)
        rug(logLOQK[hazardindexK[h],foodindex[i],]/log(10),lwd=4.5,col="green",quiet=TRUE)
        rug(logLODK[hazardindexK[h],foodindex[i],]/log(10),lwd=4.5,col="blue",quiet=TRUE)
        }
        legend("topright",cex=0.6,paste(c("prev (2.5 %):","prev (50 %):","prev (97.5 %):"),round(quantile(100*pK[,hazardindexK[h],foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  ))  
        } # end of if density
        
        if(input$selectdist=="Cumulative"){
          par(yaxt="s")
          cump <- seq(1,input$Iterations-burnin)
          cump <- cump/length(cump)
          if(scale=="Absolute"){
            cmeanK <- sort(exp(mucK[,hazardindexK[h],foodindex[i]]+0.5*sigcK[,hazardindexK[h],foodindex[i]]^2))
            cmedianK <- sort(exp(mucK[,hazardindexK[h],foodindex[i]]))
            maxx=exp(mean(mucK[,hazardindexK[h],foodindex[i]]+2*sigcK[,hazardindexK[h],foodindex[i]]))
            plot(cmeanK,cump,col="darkgray",lwd=3,main=paste(hazardnamesusedK[h],"in",foodnamesused[i]),xlab="Concentration+",ylab="",xlim=c(0,maxx),type="l") 
            lines(cmedianK,cump,lwd=3)
            # mark data points and possible LOD and LOQ values for censored data:
            rug(exp(logcK[hazardindexK[h],foodindex[i],]),lwd=2.5,col="red",quiet=TRUE)
            rug(exp(logLOQK[hazardindexK[h],foodindex[i],]),lwd=4.5,col="green",quiet=TRUE)
            rug(exp(logLODK[hazardindexK[h],foodindex[i],]),lwd=4.5,col="blue",quiet=TRUE)
            lines(ecdf(
              c(exp(logcK[hazardindexK[h],foodindex[i],]),
                exp(logLOQK[hazardindexK[h],foodindex[i],]),
                exp(logLODK[hazardindexK[h],foodindex[i],])
                )),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
            lines(ecdf(
              c(exp(logcK[hazardindexK[h],foodindex[i],]),
                exp(logLOQLimK[hazardindexK[h],foodindex[i],]),
                exp(logLODLimK[hazardindexK[h],foodindex[i],])
                )),verticals=TRUE,do.points=FALSE,lwd=2,col="blue")
          }
          if(scale=="Logarithmic"){
            maxx <- mean(mucK[,hazardindexK[h],foodindex[i]]+2*sigcK[,hazardindexK[h],foodindex[i]])
            minn <- mean(mucK[,hazardindexK[h],foodindex[i]]-2*sigcK[,hazardindexK[h],foodindex[i]])                                                               
            plot(sort(mucK[,hazardindexK[h],foodindex[i]]/log(10)),cump,lwd=3,main=paste(hazardnamesusedK[h],"in",foodnamesused[i]),xlab="log Concentration+",ylab="",xlim=c(minn/log(10),maxx/log(10)),type="l") 
            # mark data points and possible LOD and LOQ values for censored data:
            rug(logcK[hazardindexK[h],foodindex[i],]/log(10),lwd=2.5,col="red",quiet=TRUE)
            rug(logLOQK[hazardindexK[h],foodindex[i],]/log(10),lwd=4.5,col="green",quiet=TRUE)
            rug(logLODK[hazardindexK[h],foodindex[i],]/log(10),lwd=4.5,col="blue",quiet=TRUE)
            lines(ecdf(
              c(logcK[hazardindexK[h],foodindex[i],]/log(10),
                logLOQK[hazardindexK[h],foodindex[i],]/log(10),
                logLODK[hazardindexK[h],foodindex[i],]/log(10) 
                )),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
            lines(ecdf(
              c(logcK[hazardindexK[h],foodindex[i],]/log(10),
                logLOQLimK[hazardindexK[h],foodindex[i],]/log(10),
                logLODLimK[hazardindexK[h],foodindex[i],]/log(10)
                )),verticals=TRUE,do.points=FALSE,lwd=2,col="blue")
            
          }  
          legend("bottomright",inset=0.04,cex=0.6,paste(c("prev (2.5 %):","prev (50 %):","prev (97.5 %):"),round(quantile(100*pK[,hazardindexK[h],foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%") ))  
        } # end of if cumulative
       
       
       mc <- round(seq(1,input$Iterations-burnin,length=25),0)
         if(input$selectdist=="Density"){
         if(scale=="Absolute"){
         for(n in 1:25){   
         lines(seq(0,maxx,length=100),
               dlnorm(seq(0,maxx,length=100),
                      mucK[mc[n],hazardindexK[h],foodindex[i]],
                      sigcK[mc[n],hazardindexK[h],foodindex[i]]),
               col=raspberry) 
         }}
         if(scale=="Logarithmic"){
         for(n in 1:25){   
         lines(seq(minn/log(10),maxx/log(10),length=100),
               dnorm(seq(minn/log(10),maxx/log(10),length=100),
                      mucK[mc[n],hazardindexK[h],foodindex[i]]/log(10),
                      sigcK[mc[n],hazardindexK[h],foodindex[i]]/log(10)),
               col=raspberry) 
         }} 
         }  # end of if density   
       if(input$selectdist=="Cumulative"){
         par(yaxt="s")
         if(scale=="Absolute"){
           for(n in 1:25){   
             lines(seq(0,maxx*1.1,length=100),
                   plnorm(seq(0,maxx*1.1,length=100),
                          mucK[mc[n],hazardindexK[h],foodindex[i]],
                          sigcK[mc[n],hazardindexK[h],foodindex[i]]),
                   col=raspberry) 
           }}
         if(scale=="Logarithmic"){
           for(n in 1:25){   
             lines(seq(minn/log(10),maxx/log(10),length=100),
                   pnorm(seq(minn/log(10),maxx/log(10),length=100),
                         mucK[mc[n],hazardindexK[h],foodindex[i]]/log(10),
                         sigcK[mc[n],hazardindexK[h],foodindex[i]]/log(10)),
                   col=raspberry)  
           }} 
       }  # end of if cumulative
     }}
     }  # end of if nhusedK nfused
     
     # Microbiological concentrations:
     
     if((nhusedM>0)&(nfused>0)){
       # redefine dimensions if scalars:
       if((nhM==1)&(nf==1)){
         mucM <- array(mucM,dim=c(input$Iterations-burnin,1,1))
         sigcM <- array(sigcM,dim=c(input$Iterations-burnin,1,1))
         pM <- array(pM,dim=c(input$Iterations-burnin,1,1))
       }      
     for(h in 1:nhusedM){
     for(i in 1:nfused){
         if(input$selectdist=="Density"){
         if(scale=="Absolute"){
         cmeanM <- exp(mucM[,hazardindexM[h],foodindex[i]]+0.5*sigcM[,hazardindexM[h],foodindex[i]]^2)
         cmedianM <- exp(mucM[,hazardindexM[h],foodindex[i]])
         maxx = exp(mean(mucM[,hazardindexM[h],foodindex[i]]+2*sigcM[,hazardindexM[h],foodindex[i]]^2))
         plot(density(cmeanM[cmeanM<maxx]),col="darkgray",lwd=3,main=paste(hazardnamesusedM[h],"in",foodnamesused[i]),xlab="Concentration+",ylab="",xlim=c(0,maxx)) 
         polygon(density(cmeanM[cmeanM<maxx]),col="lightgray")
         lines(density(cmedianM[cmedianM<maxx]),lwd=3)
         # mark data points and possible LOD and LOQ values for censored data:
         rug(exp(logcM[hazardindexM[h],foodindex[i],]),lwd=2.5,col="red",quiet=TRUE)
         rug(exp(logLOQM[hazardindexM[h],foodindex[i],]),lwd=4.5,col="green",quiet=TRUE)
         rug(exp(logLODM[hazardindexM[h],foodindex[i],]),lwd=4.5,col="blue",quiet=TRUE)
         }
           
         if(scale=="Logarithmic"){
           maxx <- mean(mucM[,hazardindexM[h],foodindex[i]]+2*sigcM[,hazardindexM[h],foodindex[i]])
           minn <- mean(mucM[,hazardindexM[h],foodindex[i]]-2*sigcM[,hazardindexM[h],foodindex[i]])
         plot(density(mucM[,hazardindexM[h],foodindex[i]]/log(10)),col="darkgray",lwd=3,main=paste(hazardnamesusedM[h],"in",foodnamesused[i]),xlab="log Concentration+",ylab="",xlim=c(minn/log(10),maxx/log(10))) 
         polygon(density(mucM[,hazardindexM[h],foodindex[i]]/log(10)),col="lightgray")
         # mark data points and possible LOD and LOQ values for censored data:
         rug(logcM[hazardindexM[h],foodindex[i],]/log(10),lwd=2.5,col="red",quiet=TRUE)
         rug(logLOQM[hazardindexM[h],foodindex[i],]/log(10),lwd=4.5,col="green",quiet=TRUE)
         rug(logLODM[hazardindexM[h],foodindex[i],]/log(10),lwd=4.5,col="blue",quiet=TRUE)
         } # end of if logarithmic
         legend("topright",cex=0.6,paste(c("prev (2.5 %):","prev (50 %):","prev (97.5 %):"),round(quantile(100*pM[,hazardindexM[h],foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%") ))   
         } # end of if density
       
         if(input$selectdist=="Cumulative"){
           par(yaxt="s")
           cump <- seq(1,input$Iterations-burnin)
           cump <- cump/length(cump)
           if(scale=="Absolute"){
             cmeanM <- sort(exp(mucM[,hazardindexM[h],foodindex[i]]+0.5*sigcM[,hazardindexM[h],foodindex[i]]^2))
             cmedianM <- sort(exp(mucM[,hazardindexM[h],foodindex[i]]))
             maxx=exp(mean(mucM[,hazardindexM[h],foodindex[i]]+2*sigcM[,hazardindexM[h],foodindex[i]]))
             plot(cmeanM,cump,col="darkgray",lwd=3,main=paste(hazardnamesusedM[h],"in",foodnamesused[i]),xlab="Concentration+",ylab="",xlim=c(0,maxx),type="l") 
             lines(cmedianM,cump,lwd=3)
             # mark data points and possible LOD and LOQ values for censored data:
             rug(exp(logcM[hazardindexM[h],foodindex[i],]),lwd=2.5,col="red",quiet=TRUE)
             rug(exp(logLOQM[hazardindexM[h],foodindex[i],]),lwd=4.5,col="green",quiet=TRUE)
             rug(exp(logLODM[hazardindexM[h],foodindex[i],]),lwd=4.5,col="blue",quiet=TRUE)
             lines(ecdf(
               c(exp(logcM[hazardindexM[h],foodindex[i],]),
                 exp(logLOQM[hazardindexM[h],foodindex[i],]),
                 exp(logLODM[hazardindexM[h],foodindex[i],]))
             ),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
             lines(ecdf(
               c(exp(logcM[hazardindexM[h],foodindex[i],]),
                 exp(logLOQLimM[hazardindexM[h],foodindex[i],]),
                 exp(logLODLimM[hazardindexM[h],foodindex[i],]))
             ),verticals=TRUE,do.points=FALSE,lwd=2,col="blue")
             
             
           } # end of if absolute
           
           if(scale=="Logarithmic"){
             maxx <- mean(mucM[,hazardindexM[h],foodindex[i]]+2*sigcM[,hazardindexM[h],foodindex[i]])
             minn <- mean(mucM[,hazardindexM[h],foodindex[i]]-2*sigcM[,hazardindexM[h],foodindex[i]])
             plot(sort(mucM[,hazardindexM[h],foodindex[i]]/log(10)),cump,lwd=3,main=paste(hazardnamesusedM[h],"in",foodnamesused[i]),xlab="log Concentration+",ylab="",xlim=c(minn/log(10),maxx/log(10)),type="l") 
             # mark data points and possible LOD and LOQ values for censored data:
             rug(logcM[hazardindexM[h],foodindex[i],]/log(10),lwd=2.5,col="red",quiet=TRUE)
             rug(logLOQM[hazardindexM[h],foodindex[i],]/log(10),lwd=4.5,col="green",quiet=TRUE)
             rug(logLODM[hazardindexM[h],foodindex[i],]/log(10),lwd=4.5,col="blue",quiet=TRUE)
             lines(ecdf(
               c(logcM[hazardindexM[h],foodindex[i],]/log(10),
                 logLOQM[hazardindexM[h],foodindex[i],]/log(10),
                 logLODM[hazardindexM[h],foodindex[i],]/log(10))
               ),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
             lines(ecdf(
               c(logcM[hazardindexM[h],foodindex[i],]/log(10),
                 logLOQLimM[hazardindexM[h],foodindex[i],]/log(10),
                 logLODLimM[hazardindexM[h],foodindex[i],]/log(10))
             ),verticals=TRUE,do.points=FALSE,lwd=2,col="blue")
           } # end of if logarithmic
           legend("bottomright",inset=0.04,cex=0.6,paste(c("prev (2.5 %):","prev (50 %):","prev (97.5 %):"),round(quantile(100*pM[,hazardindexM[h],foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  ))
         } # end of if cumulative
       
         # plot 25 realisations of variability distributions, from the MCMC sample
         mc <- round(seq(1,input$Iterations-burnin,length=25),0)
         if(input$selectdist=="Density"){
         if(scale=="Absolute"){
         for(n in 1:25){
           lines(seq(0,maxx,length=100),
                 dlnorm(seq(0,maxx*1.1,length=100),
                        mucM[mc[n],hazardindexM[h],foodindex[i]],
                        sigcM[mc[n],hazardindexM[h],foodindex[i]]),
                 col=raspberry)
         }} # end of if absolute
           
         if(scale=="Logarithmic"){
         for(n in 1:25){
           lines(seq(minn/log(10),maxx/log(10),length=100),
                 dnorm(seq(minn/log(10),maxx/log(10),length=100),
                        mucM[mc[n],hazardindexM[h],foodindex[i]]/log(10),
                        sigcM[mc[n],hazardindexM[h],foodindex[i]]/log(10)),
                 col=raspberry) 
         }} # end of if logarithmic
         } # end of if density   
         
         if(input$selectdist=="Cumulative"){
           par(yaxt="s")
           if(scale=="Absolute"){
             for(n in 1:25){
               lines(seq(0,maxx*1.1,length=100),
                     plnorm(seq(0,maxx*1.1,length=100),
                            mucM[mc[n],hazardindexM[h],foodindex[i]],
                            sigcM[mc[n],hazardindexM[h],foodindex[i]]),
                     col=raspberry)  
             }} # end of if absolute
           
           if(scale=="Logarithmic"){
             for(n in 1:25){
               lines(seq(minn/log(10),maxx/log(10),length=100),
                     pnorm(seq(minn/log(10),maxx/log(10),length=100),
                           mucM[mc[n],hazardindexM[h],foodindex[i]]/log(10),
                           sigcM[mc[n],hazardindexM[h],foodindex[i]]/log(10)),
                     col=raspberry)
             }} # end of if logarithmic
         } # end of if cumulative   
       }}
     } # end of if nhusedM nfused  
     } # end of if selectresults == "Concentrations"
   })
   
##################################################################################

      output$distPlot2 <- renderPlot({
     # generate results based on inputs from ui.R:  Consumption amounts
     
     currentresults()
     
     theresults = input$selectresults
     
     if(is.element("Consumptions",theresults)){
       
       foodnamesused = input$thefoodnames
       nfused = length(foodnamesused)
       foodindex = match(foodnamesused,foodnames)
       
       # redefine dimensions as matrix for compatibility if they were scalars:
       if(nf==1){
       mus0 <- matrix(mus0,input$Iterations-burnin,1)  
       sigs <- matrix(sigs,input$Iterations-burnin,1) 
       Ss0 <- array(Ss0,dim=c(input$Iterations-burnin,1,1))
       ppred <- matrix(ppred,input$Iterations-burnin,1) 
       }
       
       OIM <- numeric() # observed individual mean consumptions
       if(nfused>0){
           par(mfrow=c(nfused,2),cex.lab=1.3,cex.main=1.3,yaxt="n")
           for(i in 1:nfused){
             if(input$selectdist=="Density"){
             if(input$selectscale=="Absolute"){
             
             # distributions of chronic consumptions, on consumption days (absolute per bodyweight)
             musmean <- exp(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+0.5*Ss0[,foodindex[i],foodindex[i]])
             musmedian <- exp(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2)
             maxx=exp(mean(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]])))
             plot(density(musmean[musmean<maxx]),col="darkgray",lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="C.consumption/bw+",ylab="",xlim=c(0,maxx)) 
             polygon(density(musmean[musmean<maxx]),col="lightgray")
             lines(density(musmedian[musmedian<maxx]),lwd=3)
             for(r in 1:nr){
                OIM[r]<- mean(exp(logsw[r,1:nd,foodindex[i]]),na.rm=TRUE) 
             } 
               OIM<-OIM[!is.na(OIM)]
               # mark data points: (observed individual means)
               rug(OIM,lwd=2.5,col="red",quiet=TRUE)
             
             mc <- round(seq(1,input$Iterations-burnin,length=25),0)
             for(n in 1:25){
               lines(seq(0,maxx,length=100),
                     dlnorm(seq(0,maxx,length=100),
                            mus0[mc[n],foodindex[i]]+0.5*sigs[mc[n],foodindex[i]]^2,
                            sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])),
                     col=raspberry)
             }
             legend("topright",cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%") )) 
             
             # distribution of acute consumptions, on consumption days (absolute):
             smean <- exp(mus0[,foodindex[i]]+0.5*Ss0[,foodindex[i],foodindex[i]]+0.5*sigs[,foodindex[i]]^2+muw+0.5*sigw^2)
             smedian <- exp(mus0[,foodindex[i]]+muw)
             maxx=exp(mean(mus0[,foodindex[i]]+muw+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]+sigs[,foodindex[i]]^2+sigw^2))) 
             plot(density(smean[smean<maxx]),col="darkgray",lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="A.consumption+",ylab="",xlim=c(0,maxx)) 
             polygon(density(smean[smean<maxx]),col="lightgray")
             lines(density(smedian[smedian<maxx]),lwd=3)
             # mark data points: (individual acute consumptions)
             rug(exp(logs[1:nr,1:nd,foodindex[i]]),lwd=2.5,col="red",quiet=TRUE)
             
             mc <- round(seq(1,input$Iterations-burnin,length=25),0)
             for(n in 1:25){
               lines(seq(0,maxx,length=100),
                     dlnorm(seq(0,maxx,length=100),
                            mus0[mc[n],foodindex[i]]+muw[mc[n]],
                            sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+sigs[mc[n],foodindex[i]]^2+sigw[mc[n]]^2)),
                     col=raspberry)  
             }
             legend("topright",cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  )) 
             } # end of if absolute
               
             if(input$selectscale=="Logarithmic"){
               
               # distributions of chronic consumptions, on consumption days (log per bodyweight)
               logmusmean <- mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2
               maxx <- mean(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]))
               minn <- mean(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2-3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]))
               plot(density(logmusmean/log(10)),col="darkgray",lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="log C.consumption/bw+",ylab="",xlim=c(minn/log(10),maxx/log(10))) 
               polygon(density(logmusmean/log(10)),col="lightgray")
               for(r in 1:nr){
                 OIM[r]<- log(mean(exp(logsw[r,1:nd,foodindex[i]]),na.rm=TRUE)) 
               } 
               OIM<-OIM[!is.na(OIM)]
               # mark data points: (observed individual means, in log-scale)
               rug(OIM/log(10),lwd=2.5,col="red",quiet=TRUE)
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){
                 lines(seq(minn/log(10),maxx/log(10),length=100),
                       dnorm(seq(minn/log(10),maxx/log(10),length=100),
                              (mus0[mc[n],foodindex[i]]+0.5*sigs[mc[n],foodindex[i]]^2)/log(10),
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])/log(10)),
                       col=raspberry)  
               }
               legend("topright",cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  ))    
               
               # distribution of acute consumptions, on consumption days (log):
               logsmean <- mus0[,foodindex[i]]+muw
               maxx <- mean(mus0[,foodindex[i]]+muw+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]+sigs[,foodindex[i]]^2+sigw^2))
               minn <- mean(mus0[,foodindex[i]]+muw-3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]+sigs[,foodindex[i]]^2+sigw^2))
               plot(density(logsmean/log(10)),col="darkgray",lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="log A.consumption+",ylab="",xlim=c(minn/log(10),maxx/log(10))) 
               polygon(density(logsmean/log(10)),col="lightgray")
               # mark data points: (individual acute consumptions, in log-scale)
               rug(logs[1:nr,1:nd,foodindex[i]]/log(10),lwd=2.5,col="red",quiet=TRUE)
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){
                 lines(seq(minn/log(10),maxx/log(10),length=100),
                       dnorm(seq(minn/log(10),maxx/log(10),length=100),
                              (mus0[mc[n],foodindex[i]]+muw[mc[n]])/log(10),
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+sigs[mc[n],foodindex[i]]^2+sigw[mc[n]]^2)/log(10)),
                       col=raspberry)  
               }
      
               legend("topright",cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%") )) 
              
             } # end of if logarithmic
             } # end of if density
             
             if(input$selectdist=="Cumulative"){
               par(yaxt="s")
               cump <- seq(1,input$Iterations-burnin)
               cump <- cump/length(cump)
               if(input$selectscale=="Absolute"){
               
               # distributions of chronic consumptions (absolute per bodyweight)
               musmean <- sort(exp(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+0.5*Ss0[,foodindex[i],foodindex[i]]))
               musmedian <- sort(exp(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2))
               maxx <- exp(mean(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]])))
               plot(musmean,cump,col="darkgray",lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="C.consumption/bw+",ylab="",xlim=c(0,maxx),type="l")
               lines(musmedian,cump,lwd=3)
               for(r in 1:nr){
                 OIM[r]<- mean(exp(logsw[r,1:nd,foodindex[i]]),na.rm=TRUE) 
               } 
               OIM<-OIM[!is.na(OIM)]
               # mark data points: (observed individual means)
               rug(OIM,lwd=2.5,col="red",quiet=TRUE)
               lines(ecdf(OIM),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){
                 lines(seq(0,maxx,length=100),
                       plnorm(seq(0,maxx,length=100),
                              mus0[mc[n],foodindex[i]]+0.5*sigs[mc[n],foodindex[i]]^2,
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])),
                       col=raspberry)  
               }
               legend("bottomright",inset=0.04,cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  )) 
               
               # distribution of acute consumptions (absolute):
               smean <- sort(exp(mus0[,foodindex[i]]+0.5*Ss0[,foodindex[i],foodindex[i]]+0.5*sigs[,foodindex[i]]^2+muw+0.5*sigw^2))
               smedian <- sort(exp(mus0[,foodindex[i]]+muw))
               maxx <- exp(mean(mus0[,foodindex[i]]+muw+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]+sigs[,foodindex[i]]^2+sigw^2)))
               plot(smean,cump,col="darkgray",lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="A.consumption+",ylab="",xlim=c(0,maxx),type="l") 
               lines(smedian,cump,lwd=3)
               # mark data points: (individual acute consumptions)
               rug(exp(logs[1:nr,1:nd,foodindex[i]]),lwd=2.5,col="red",quiet=TRUE)
               lines(ecdf(exp(logs[1:nr,1:nd,foodindex[i]])),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){
                 lines(seq(0,maxx,length=100),
                       plnorm(seq(0,maxx,length=100),
                              mus0[mc[n],foodindex[i]]+muw[mc[n]],
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+sigs[mc[n],foodindex[i]]^2+sigw[mc[n]]^2)),
                       col=raspberry)  
               }
               legend("bottomright",inset=0.04,cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  )) 
             } # end of if absolute
               
             if(input$selectscale=="Logarithmic"){
               # distributions of chronic consumptions (log per bodyweight)
               logmusmean <- sort(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2)
               maxx <- mean(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]))
               minn <- mean(mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2-3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]))
               plot(logmusmean/log(10),cump,lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="log C.consumption/bw+",ylab="",xlim=c(minn/log(10),maxx/log(10)),type="l")
               for(r in 1:nr){
                 OIM[r]<- log(mean(exp(logsw[r,1:nd,foodindex[i]]),na.rm=TRUE)) 
               } 
               OIM<-OIM[!is.na(OIM)]
               # mark data points: (observed individual means, in log-scale)
               rug(OIM/log(10),lwd=2.5,col="red",quiet=TRUE)
               lines(ecdf(OIM/log(10)),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){
                 lines(seq(minn/log(10),maxx/log(10),length=100),
                       pnorm(seq(minn/log(10),maxx/log(10),length=100),
                              (mus0[mc[n],foodindex[i]]+0.5*sigs[mc[n],foodindex[i]]^2)/log(10),
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])/log(10)),
                       col=raspberry)  
               }
               legend("bottomright",inset=0.04,cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  )) 
               # distribution of acute consumptions (log):
               logsmean <- sort(mus0[,foodindex[i]]+muw)
               maxx <- mean(mus0[,foodindex[i]]+muw+3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]+sigs[,foodindex[i]]^2+sigw^2))
               minn <- mean(mus0[,foodindex[i]]+muw-3.1*sqrt(Ss0[,foodindex[i],foodindex[i]]+sigs[,foodindex[i]]^2+sigw^2))
               plot(logsmean/log(10),cump,lwd=3,main=paste(foodnamesused[i],"consumption"),xlab="log A.consumption+",ylab="",xlim=c(minn/log(10),maxx/log(10)),type="l")
               # mark data points: (individual acute consumptions, in log-scale)
               rug(logs[1:nr,1:nd,foodindex[i]]/log(10),lwd=2.5,col="red",quiet=TRUE)
               lines(ecdf(logs[1:nr,1:nd,foodindex[i]]/log(10)),verticals=TRUE,do.points=FALSE,lwd=2,col="red")
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){
                 lines(seq(minn/log(10),maxx/log(10),length=100),
                       pnorm(seq(minn/log(10),maxx/log(10),length=100),
                              (mus0[mc[n],foodindex[i]]+muw[mc[n]])/log(10),
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+sigs[mc[n],foodindex[i]]^2+sigw[mc[n]]^2)/log(10)),
                       col=raspberry) 
               }
               legend("bottomright",inset=0.04,cex=0.6,paste(c("freq (2.5 %):","freq (50 %):","freq (97.5 %):"),round(quantile(100*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  )) 
             } # end of if logarithmic
             } # end of if cumulative
               
           } # end of for nfused
       }  # end of if nfused>0 
       
       
     } # end of if selectresults == "Consumptions"
   })
   
###########################################################################   
   
   output$distPlot3 =  renderPlot({         
     # generate results based on inputs from ui.R:  
     # Exposures
              
     currentresults()
     
     theresults = input$selectresults
     
     if(is.element("Exposures",theresults)){
    
     scale = input$selectscale        
     foodnamesused = input$thefoodnames
     nfused = length(foodnamesused)
     foodindex = match(foodnamesused,foodnames)
     
     hazardnamesused = input$thehazardnames
     hazardtypesused = hazardtypes[is.element(hazardnames,hazardnamesused)]
     nhused = length(hazardnamesused)
     hazardnamesK = hazardnames[hazardtypes=="K"]
     hazardnamesM = hazardnames[hazardtypes=="M"]
     hazardnamesusedK = hazardnamesused[hazardtypesused=="K"]
     hazardnamesusedM = hazardnamesused[hazardtypesused=="M"]
     nhusedK = length(hazardnamesusedK)
     nhusedM = length(hazardnamesusedM)
     hazardindex = match(hazardnamesused,hazardnames)
     hazardindexK = match(hazardnamesusedK,hazardnamesK)
     hazardindexM = match(hazardnamesusedM,hazardnamesM)
     Rall = input$factor  # adjustment factor for concentrations
     Pall = input$pfactor # adjustment factor for prevalences
     
     # redefine dimensions as matrix for compatibility if they were scalars:
     if(nf==1){
       mus0 <- matrix(mus0,input$Iterations-burnin,1)  
       sigs <- matrix(sigs,input$Iterations-burnin,1) 
       Ss0 <- array(Ss0,dim=c(input$Iterations-burnin,1,1))
       ppred <- matrix(ppred,input$Iterations-burnin,1) 
     } 
     
     if( (nhused>0)&(nfused>0) ){par(mfrow=c(nhused,nfused),cex.lab=1.3,cex.main=1.3,yaxt="n")}
     
     # Chemical exposures
     if((nhusedK>0)&(nfused>0)){
       # redefine dimensions if scalars:
       if((nhK==1)&(nf==1)){
         mucK <- array(mucK,dim=c(input$Iterations-burnin,1,1))
         sigcK <- array(sigcK,dim=c(input$Iterations-burnin,1,1))
         pK <- array(pK,dim=c(input$Iterations-burnin,1,1))
       } 
       RK = matrix(NA,nf,nhK) # factors for concentrations
       RK[1:nf,1:nhK] = Rall[1:nf,is.element(hazardnames,hazardnamesusedK)]
       logRK = log(RK)
       PK = matrix(NA,nf,nhK) # factors for prevalence
       PK[1:nf,1:nhK] = Pall[1:nf,is.element(hazardnames,hazardnamesusedK)]
       
       for(h in 1:nhusedK){
         for(i in 1:nfused){
           if(input$selectdist=="Density"){
           ############## exposure.chronicKbw
           # plot posterior of the mean & median exposure/bw 
           # (expected chronic exposure for anyone)
           if(scale=="Absolute"){
           meanexposure <- exp(logRK[foodindex[i],hazardindexK[h]]+
                               mus0[,foodindex[i]]+
                               0.5*sigs[,foodindex[i]]^2+
                               mucK[,hazardindexK[h],foodindex[i]]+
                               0.5*sigcK[,hazardindexK[h],foodindex[i]]^2+
                               0.5*Ss0[,foodindex[i],foodindex[i]])
           maxx <- exp(quantile(logRK[foodindex[i],hazardindexK[h]]+
                                  mus0[,foodindex[i]]+
                                  0.5*sigs[,foodindex[i]]^2+
                                  mucK[,hazardindexK[h],foodindex[i]]+
                                  0.5*sigcK[,hazardindexK[h],foodindex[i]]^2+
                                  2*sqrt(Ss0[,foodindex[i],foodindex[i]]),0.5,names=FALSE))
           medianexposure <- exp(logRK[foodindex[i],hazardindexK[h]]+
                                 mus0[,foodindex[i]]+
                                 0.5*sigs[,foodindex[i]]^2+
                                 mucK[,hazardindexK[h],foodindex[i]]+
                                 0.5*sigcK[,hazardindexK[h],foodindex[i]]^2)
           plot(density(meanexposure[meanexposure<maxx]),main=paste(hazardnamesusedK[h],"from",foodnamesused[i],"(chronic)"),xlab="C.exposure/bw+",ylab="",xlim=c(0,maxx),lwd=3) 
           polygon(density(meanexposure[meanexposure<maxx]),col="lightgray")
           lines(density(medianexposure[medianexposure<maxx]),lwd=3)
           mc <- round(seq(1,input$Iterations-burnin,length=25),0)
           for(n in 1:25){ # variability distributions for chronic exposure
             lines(seq(0,maxx,length=100),
                dlnorm(seq(0,maxx,length=100),
                logRK[foodindex[i],hazardindexK[h]]
                +mus0[mc[n],foodindex[i]]
                +0.5*sigs[mc[n],foodindex[i]]^2
                +mucK[mc[n],hazardindexK[h],foodindex[i]]
                +0.5*sigcK[mc[n],hazardindexK[h],foodindex[i]]^2,
                sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])),
                   col=raspberry) 
           } # end of n
           } # end of if absolute
             
           if(scale=="Logarithmic"){
             logmeanexposure <- logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2+mucK[,hazardindexK[h],foodindex[i]]+0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
             minn <- quantile(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                               +0.5*sigs[,foodindex[i]]^2
                               +mucK[,hazardindexK[h],foodindex[i]]
                               +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
                               -2*sqrt(Ss0[,foodindex[i],foodindex[i]]),0.025,names=FALSE)
             maxx <- quantile(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                              +0.5*sigs[,foodindex[i]]^2    
                              +mucK[,hazardindexK[h],foodindex[i]]
                              +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
                              +2*sqrt(Ss0[,foodindex[i],foodindex[i]]),0.975,names=FALSE) 
             plot(density(logmeanexposure/log(10)),main=paste(hazardnamesusedK[h],"from",foodnamesused[i],"(chronic)"),xlab="log (C.exposure/bw+)",ylab="",xlim=c(minn/log(10),maxx/log(10)),lwd=3) 
             polygon(density(logmeanexposure/log(10)),col="lightgray")
             mc <- round(seq(1,input$Iterations-burnin,length=25),0)   
             for(n in 1:25){ # variability distributions for chronic exposures (log) 
               lines(seq(minn/log(10),maxx/log(10),length=100),
                     dnorm(seq(minn/log(10),maxx/log(10),length=100),
                     (logRK[foodindex[i],hazardindexK[h]]+mus0[mc[n],foodindex[i]]
                      +0.5*sigs[mc[n],foodindex[i]]^2
                      +mucK[mc[n],hazardindexK[h],foodindex[i]]
                      +0.5*sigcK[mc[n],hazardindexK[h],foodindex[i]]^2)/log(10),
                     sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]))/log(10),
                     col=raspberry) 
             } # end of for n
           } # end of if logarithmic
           legend("topright",cex=0.6,paste(c("expo (2.5 %):","expo (50 %):","expo (97.5 %)"),round(quantile(100*PK[foodindex[i],hazardindexK[h]]*pK[,hazardindexK[h],foodindex[i]]*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  ))
           } # end of if density
           
           if(input$selectdist=="Cumulative"){
             par(yaxt="s")
             cump <- seq(1,input$Iterations-burnin)
             cump <- cump/length(cump)
             ############## exposure.chronicKbw
             # plot posterior of the mean & median exposure/bw 
             # (expected chronic exposure for anyone)
             if(scale=="Absolute"){
               meanexposure <- sort(
                 exp(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                     +0.5*sigs[,foodindex[i]]^2
                     +mucK[,hazardindexK[h],foodindex[i]]
                     +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
                     +0.5*Ss0[,foodindex[i],foodindex[i]]))
               maxx <- exp(quantile(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                                    +0.5*sigs[,foodindex[i]]^2
                                    +mucK[,hazardindexK[h],foodindex[i]]
                                    +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
                                    +2*sqrt(Ss0[,foodindex[i],foodindex[i]]),0.5,names=FALSE))
               medianexposure <- sort(
                 exp(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                     +0.5*sigs[,foodindex[i]]^2
                     +mucK[,hazardindexK[h],foodindex[i]]
                     +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2))
               plot(meanexposure,cump,main=paste(hazardnamesusedK[h],"from",foodnamesused[i],"(chronic)"),xlab="C.exposure/bw+",ylab="",xlim=c(0,maxx),lwd=3,type="l",col="darkgray") 
               lines(medianexposure,cump,lwd=3)
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){ # variability distributions for chronic exposure
                 lines(seq(0,maxx*1.1,length=100),
                       plnorm(seq(0,maxx*1.1,length=100),
                              logRK[foodindex[i],hazardindexK[h]]+mus0[mc[n],foodindex[i]]
                              +0.5*sigs[mc[n],foodindex[i]]^2
                              +mucK[mc[n],hazardindexK[h],foodindex[i]]
                              +0.5*sigcK[mc[n],hazardindexK[h],foodindex[i]]^2,
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])),
                       col=raspberry)  
               } # end of n
               # plot empirically generated cumulative exposure distributions
               OIM <- numeric()
               for(r in 1:nr){
                 OIM[r]<- mean(exp(logsw[r,1:nd,foodindex[i]]),na.rm=TRUE) 
               } 
               OIM<-OIM[!is.na(OIM)]
               # collect exact measurements & 
               # and as upper bounds those between LOD-LOQ & <LOD 
               concentrations <- exp(c(logcK[hazardindexK[h],foodindex[i],],
                                       logLOQK[hazardindexK[h],foodindex[i],],
                                       logLODK[hazardindexK[h],foodindex[i],]))
               # and using lower bounds:
               concentrations0 <- exp(c(logcK[hazardindexK[h],foodindex[i],],
                                       logLOQLimK[hazardindexK[h],foodindex[i],],
                                       logLODLimK[hazardindexK[h],foodindex[i],]))
               concentrations <- concentrations[!is.na(concentrations)]
               concentrations0 <- concentrations0[!is.na(concentrations0)]
               
               for(resample in 1:20){   
                 # create 20 replicate ('bootstrap') data with original nsample:
                 sampleOIM <- sample(OIM,length(OIM),replace=TRUE)
                 samplecon <- sample(concentrations,length(concentrations),replace=TRUE)
                 samplecon0 <- sample(concentrations0,length(concentrations0),replace=TRUE)
                 # create 80000 simulations from each replicated data:
                 sampleOIM <- sample(sampleOIM,80000,replace=TRUE)
                 samplecon <- sample(samplecon,80000,replace=TRUE)
                 samplecon0 <- sample(samplecon0,80000,replace=TRUE)
                 lines(ecdf(sampleOIM*mean(samplecon)*RK[foodindex[i],hazardindexK[h]]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="red")
                 lines(ecdf(sampleOIM*mean(samplecon0)*RK[foodindex[i],hazardindexK[h]]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="blue")
               }
             } # end of if absolute
             
             if(scale=="Logarithmic"){
               logmeanexposure <- sort(
                 logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                 +0.5*sigs[,foodindex[i]]^2
                 +mucK[,hazardindexK[h],foodindex[i]]
                 +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2)
               minn <- quantile(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                            +0.5*sigs[,foodindex[i]]^2
                            +mucK[,hazardindexK[h],foodindex[i]]
                            +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
                            -2*sqrt(Ss0[,foodindex[i],foodindex[i]]),0.025,names=FALSE)
               maxx <- quantile(logRK[foodindex[i],hazardindexK[h]]+mus0[,foodindex[i]]
                            +0.5*sigs[,foodindex[i]]^2    
                            +mucK[,hazardindexK[h],foodindex[i]]
                            +0.5*sigcK[,hazardindexK[h],foodindex[i]]^2
                            +2*sqrt(Ss0[,foodindex[i],foodindex[i]]),0.975,names=FALSE)
               plot(logmeanexposure/log(10),cump,main=paste(hazardnamesusedK[h],"from",foodnamesused[i],"(chronic)"),xlab="log (C.exposure/bw+)",ylab="",xlim=c(minn/log(10),maxx/log(10)),lwd=3,type="l") 
   
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)   
               for(n in 1:25){ # variability distributions for chronic exposures (log) 
                 lines(seq(minn/log(10),maxx/log(10),length=100),
                       pnorm(seq(minn/log(10),maxx/log(10),length=100),
                             (logRK[foodindex[i],hazardindexK[h]]+mus0[mc[n],foodindex[i]]
                             +0.5*sigs[mc[n],foodindex[i]]^2
                             +mucK[mc[n],hazardindexK[h],foodindex[i]]
                             +0.5*sigcK[mc[n],hazardindexK[h],foodindex[i]]^2)/log(10),
                             sqrt(Ss0[mc[n],foodindex[i],foodindex[i]])/log(10)),
                       col=raspberry) 
               } # end of for n
               # plot empirically generated cumulative exposure distributions
               OIM <- numeric()
               for(r in 1:nr){
                 OIM[r]<- mean(exp(logsw[r,1:nd,foodindex[i]]),na.rm=TRUE) 
               } 
               OIM<-OIM[!is.na(OIM)]
               # collect exact measurements & 
               # and as upper bounds those between LOD-LOQ & <LOD 
               concentrations <- exp(c(logcK[hazardindexK[h],foodindex[i],],
                                       logLOQK[hazardindexK[h],foodindex[i],],
                                       logLODK[hazardindexK[h],foodindex[i],]))
               # and using lower bounds:
               concentrations0 <- exp(c(logcK[hazardindexK[h],foodindex[i],],
                                       logLOQLimK[hazardindexK[h],foodindex[i],],
                                       logLODLimK[hazardindexK[h],foodindex[i],]))
               concentrations <- concentrations[!is.na(concentrations)]
               concentrations0 <- concentrations0[!is.na(concentrations0)]
               
               for(resample in 1:20){
                 # create 20 replicate ('bootstrap') data with original nsample:
                 sampleOIM <- sample(OIM,length(OIM),replace=TRUE)
                 samplecon <- sample(concentrations,length(concentrations),replace=TRUE)
                 samplecon0 <- sample(concentrations0,length(concentrations),replace=TRUE)
                 # create 80000 simulations from each replicated data:
                 sampleOIM <- sample(sampleOIM,80000,replace=TRUE)
                 samplecon <- sample(samplecon,80000,replace=TRUE)
                 samplecon0 <- sample(samplecon0,80000,replace=TRUE)
                 lines(ecdf(log(sampleOIM*mean(samplecon)*RK[foodindex[i],hazardindexK[h]])/log(10)),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="red")
                 lines(ecdf(log(sampleOIM*mean(samplecon0)*RK[foodindex[i],hazardindexK[h]])/log(10)),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="blue")
               }
             } # end of if logarithmic    
             legend("bottomright",cex=0.6,paste(c("expo (2.5 %):","expo (50 %):","expo (97.5 %):"),round(quantile(100*PK[foodindex[i],hazardindexK[h]]*pK[,hazardindexK[h],foodindex[i]]*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%") ))               
           } # end of if cumulative
         }} # end of for nhused nfused
     } 

     # Microbial exposures:  
     
     if((nhusedM>0)&(nfused>0)){
       # redefine dimensions if scalars:
       if((nhM==1)&(nf==1)){
         mucM <- array(mucM,dim=c(input$Iterations-burnin,1,1))
         sigcM <- array(sigcM,dim=c(input$Iterations-burnin,1,1))
         pM <- array(pM,dim=c(input$Iterations-burnin,1,1))
       } 
       RM = matrix(NA,nf,nhM) # factors for concentration
       RM[1:nf,1:nhM] = Rall[1:nf,is.element(hazardnames,hazardnamesusedM)]
       logRM = log(RM)
       PM = matrix(NA,nf,nhM) # factors for prevalence
       PM[1:nf,1:nhM] = Pall[1:nf,is.element(hazardnames,hazardnamesusedM)]
       
       for(h in 1:nhusedM){
         for(i in 1:nfused){
           if(input$selectdist=="Density"){
           ############## exposure.acuteM
           # plot posterior of the mean & median exposure 
           # (expected acute exposure for anyone)
           if(input$selectscale=="Absolute"){   
           meanexposure <- exp(logRM[foodindex[i],hazardindexM[h]]
                             +mus0[,foodindex[i]]
                             +mucM[,hazardindexM[h],foodindex[i]]
                             +muw
                             +0.5*Ss0[,foodindex[i],foodindex[i]]
                             +0.5*sigs[,foodindex[i]]^2
                             +0.5*sigcM[,hazardindexM[h],foodindex[i]]^2 
                             +0.5*sigw^2 )
           maxx <- exp(quantile(logRM[foodindex[i],hazardindexM[h]]+
                                  mus0[,foodindex[i]]+
                                  mucM[,hazardindexM[h],foodindex[i]]+
                                  muw+2*sqrt(Ss0[,foodindex[i],foodindex[i]]+
                                                 sigs[,foodindex[i]]^2+
                                                 sigcM[,hazardindexM[h],foodindex[i]]^2+
                                                 sigw^2),0.5,names=FALSE))
           medianexposure <- exp(logRM[foodindex[i],hazardindexM[h]]+
                                 mus0[,foodindex[i]]+
                                 mucM[,hazardindexM[h],foodindex[i]]+muw)

           plot(density(meanexposure[meanexposure<maxx]),main=paste(hazardnamesusedM[h],"from",foodnamesused[i],"(acute)"),xlab="A.exposure+",ylab="",xlim=c(0,maxx) ) 
           polygon(density(meanexposure[meanexposure<maxx]),col="lightgray")
           lines(density(medianexposure[medianexposure<maxx]),lwd=3)
                 
           mc <- round(seq(1,input$Iterations-burnin,length=25),0)
           for(n in 1:25){ # variability distributions for acute exposure
            lines(seq(0,maxx*1.1,length=100),
             dlnorm(seq(0,maxx*1.1,length=100),
             logRM[foodindex[i],hazardindexM[h]]+
               mus0[mc[n],foodindex[i]]+
               mucM[mc[n],hazardindexM[h],foodindex[i]]+
               muw[mc[n]],
             sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+
                    sigs[mc[n],foodindex[i]]^2+
                    sigcM[mc[n],hazardindexM[h],foodindex[i]]^2+
                    sigw[mc[n]]^2)),
             col=raspberry) 
           }
           } # end of if absolute
             
           if(input$selectscale=="Logarithmic"){
             logmeanexposure <- logRM[foodindex[i],hazardindexM[h]]+
                              mus0[,foodindex[i]]+mucM[,hazardindexM[h],foodindex[i]]+muw
             maxx = quantile(logRM[foodindex[i],hazardindexM[h]]
                         +mus0[,foodindex[i]]
                         +mucM[,hazardindexM[h],foodindex[i]]
                         +muw+2*sqrt(Ss0[,foodindex[i],foodindex[i]]
                                           +sigs[,foodindex[i]]^2
                                           +sigcM[,hazardindexM[h],foodindex[i]]^2
                                           +sigw^2),0.025,names=FALSE)
             minn = quantile(logRM[foodindex[i],hazardindexM[h]]+mus0[,foodindex[i]]
                         +mucM[,hazardindexM[h],foodindex[i]]
                         +muw-2*sqrt(Ss0[,foodindex[i],foodindex[i]]
                                           +sigs[,foodindex[i]]^2
                                           +sigcM[,hazardindexM[h],foodindex[i]]^2
                                           +sigw^2),0.975,names=FALSE)
             plot(density(logmeanexposure/log(10)),main=paste(hazardnamesusedM[h],"from",foodnamesused[i],"(acute)"),xlab="log (A.exposure+)",ylab="",xlim=c(minn/log(10),maxx/log(10)) ) 
             polygon(density(logmeanexposure/log(10)),col="lightgray")
             mc <- round(seq(1,input$Iterations-burnin,length=25),0)
             for(n in 1:25){ # variability distributions for acute exposure
               lines(seq(minn/log(10),maxx/log(10),length=100),
                     dnorm(seq(minn/log(10),maxx/log(10),length=100),
                            (logRM[foodindex[i],hazardindexM[h]]+
                               mus0[mc[n],foodindex[i]]+
                               mucM[mc[n],hazardindexM[h],foodindex[i]]+
                               muw[mc[n]])/log(10),
                            sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+
                                   sigs[mc[n],foodindex[i]]^2+
                                   sigcM[mc[n],hazardindexM[h],foodindex[i]]^2+
                                   sigw[mc[n]]^2)/log(10)),
                     col=raspberry) 
             }   
           } # end of if logarithmic
           legend("topright",cex=0.6,paste(c("expo (2.5 %):","expo (50 %):","expo (97.5 %):"),round(quantile(100*PM[foodindex[i],hazardindexM[h]]*pM[,hazardindexM[h],foodindex[i]]*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  ))
           } # end of if density
           
           if(input$selectdist=="Cumulative"){
             par(yaxt="s")
             cump=seq(1,input$Iterations-burnin)
             cump=cump/length(cump)
             if(input$selectscale=="Absolute"){
             meanexposure <- sort(exp(logRM[foodindex[i],hazardindexM[h]]
                               +mus0[,foodindex[i]]
                               +mucM[,hazardindexM[h],foodindex[i]]
                               +muw
                               +0.5*Ss0[,foodindex[i],foodindex[i]]
                               +0.5*sigs[,foodindex[i]]^2
                               +0.5*sigcM[,hazardindexM[h],foodindex[i]]^2 
                               +0.5*sigw^2 ))
             maxx <- exp(quantile(logRM[foodindex[i],hazardindexM[h]]+
                                    mus0[,foodindex[i]]+
                                    mucM[,hazardindexM[h],foodindex[i]]+
                                    muw+2*sqrt(Ss0[,foodindex[i],foodindex[i]]+
                                                   sigs[,foodindex[i]]^2+
                                                   sigcM[,hazardindexM[h],foodindex[i]]^2+
                                                   sigw^2),0.5,names=FALSE))
             medianexposure <- sort(exp(logRM[foodindex[i],hazardindexM[h]]+
                                        mus0[,foodindex[i]]+
                                        mucM[,hazardindexM[h],foodindex[i]]+
                                        muw))
             
             plot(meanexposure,cump,main=paste(hazardnamesusedM[h],"from",foodnamesused[i],"(acute)"),xlab="A.exposure+",ylab="",xlim=c(0,maxx),type="l",lwd=3,col="darkgray") 
             lines(medianexposure,cump,lwd=3)
             
             mc <- round(seq(1,input$Iterations-burnin,length=25),0)
             for(n in 1:25){ # variability distributions for acute exposure
               lines(seq(0,maxx*1.1,length=100),
                     plnorm(seq(0,maxx*1.1,length=100),
                            logRM[foodindex[i],hazardindexM[h]]+
                              mus0[mc[n],foodindex[i]]+
                              mucM[mc[n],hazardindexM[h],foodindex[i]]+
                              muw[mc[n]],
                            sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+
                                   sigs[mc[n],foodindex[i]]^2+
                                   sigcM[mc[n],hazardindexM[h],foodindex[i]]^2+
                                   sigw[mc[n]]^2)),
                     col=raspberry) 
             }
             # plot empirically generated cumulative exposure distributions
             servings <- exp(logs[1:nr,1:nd,foodindex[i]])
             # collect exact measurements & 
             # and as upper bounds those between LOD-LOQ & <LOD 
             concentrations <- exp(c(logcM[hazardindexM[h],foodindex[i],],
                                     logLOQM[hazardindexM[h],foodindex[i],],
                                     logLODM[hazardindexM[h],foodindex[i],]))
             # and using lower bounds:
             concentrations0 <- exp(c(logcM[hazardindexM[h],foodindex[i],],
                                      logLOQLimM[hazardindexM[h],foodindex[i],], 
                                      logLODLimM[hazardindexM[h],foodindex[i],]))
             
             servings <- servings[!is.na(servings)]
             concentrations <- concentrations[!is.na(concentrations)]
             concentrations0 <- concentrations0[!is.na(concentrations0)]
             for(resample in 1:20){
               # create 20 replicate ('bootstrap') data with original nsample:
               sampleser <- sample(servings,length(servings),replace=TRUE)
               samplecon <- sample(concentrations,length(concentrations),replace=TRUE)
               samplecon0 <- sample(concentrations0,length(concentrations0),replace=TRUE)
               # create 80000 simulations from each replicated data:
               sampleser <- sample(sampleser,80000,replace=TRUE)
               samplecon <- sample(samplecon,80000,replace=TRUE)
               samplecon0 <- sample(samplecon0,80000,replace=TRUE)
               lines(ecdf(sampleser*samplecon*RM[foodindex[i],hazardindexM[h]]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="red")
               lines(ecdf(sampleser*samplecon0*RM[foodindex[i],hazardindexM[h]]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="blue")
             }
             } # end of if absolute
             
             if(input$selectscale=="Logarithmic"){
               logmeanexposure <- sort(logRM[foodindex[i],hazardindexM[h]]
                                      +mus0[,foodindex[i]]
                                      +mucM[,hazardindexM[h],foodindex[i]]
                                      +muw)
               maxx <- quantile(logRM[foodindex[i],hazardindexM[h]]+mus0[,foodindex[i]]
                         +mucM[,hazardindexM[h],foodindex[i]]
                         +muw+2*sqrt(Ss0[,foodindex[i],foodindex[i]]
                                      +sigs[,foodindex[i]]^2
                                      +sigcM[,hazardindexM[h],foodindex[i]]^2
                                      +sigw^2),0.025,names=FALSE)
               minn <- quantile(logRM[foodindex[i],hazardindexM[h]]+mus0[,foodindex[i]]
                         +mucM[,hazardindexM[h],foodindex[i]]
                         +muw-2*sqrt(Ss0[,foodindex[i],foodindex[i]]
                            +sigs[,foodindex[i]]^2
                            +sigcM[,hazardindexM[h],foodindex[i]]^2
                            +sigw^2),0.975,names=FALSE)
               plot(logmeanexposure/log(10),cump,main=paste(hazardnamesusedM[h],"from",foodnamesused[i],"(acute)"),xlab="log (A.exposure+)",ylab="",xlim=c(minn/log(10),maxx/log(10)),lwd=3,type="l") 
               
               mc <- round(seq(1,input$Iterations-burnin,length=25),0)
               for(n in 1:25){ # variability distributions for acute exposure
                 lines(seq(minn/log(10),maxx/log(10),length=100),
                       pnorm(seq(minn/log(10),maxx/log(10),length=100),
                              (logRM[foodindex[i],hazardindexM[h]]+
                                 mus0[mc[n],foodindex[i]]+
                                 mucM[mc[n],hazardindexM[h],foodindex[i]]+
                                 muw[mc[n]])/log(10),
                              sqrt(Ss0[mc[n],foodindex[i],foodindex[i]]+
                                     sigs[mc[n],foodindex[i]]^2+
                                     sigcM[mc[n],hazardindexM[h],foodindex[i]]^2+
                                     sigw[mc[n]]^2)/log(10)),
                       col=raspberry) 
               } 
               # plot empirically generated cumulative exposure distributions
               servings <- exp(logs[1:nr,1:nd,foodindex[i]])
               # collect exact measurements & 
               # and as upper bounds those between LOD-LOQ & <LOD 
               concentrations <- exp(c(logcM[hazardindexM[h],foodindex[i],],
                                       logLOQM[hazardindexM[h],foodindex[i],],
                                       logLODM[hazardindexM[h],foodindex[i],]))
               # and using lower bounds
               concentrations0 <- exp(c(logcM[hazardindexM[h],foodindex[i],],
                                        logLOQLimM[hazardindexM[h],foodindex[i],],
                                        logLODLimM[hazardindexM[h],foodindex[i],]))
               servings <- servings[!is.na(servings)]
               concentrations <- concentrations[!is.na(concentrations)]
               concentrations0 <- concentrations0[!is.na(concentrations0)]
               
               for(resample in 1:20){
               # create 20 replicate ('bootstrap') data with original nsample:   
               sampleser <- sample(servings,length(servings),replace=TRUE)
               samplecon <- sample(concentrations,length(concentrations),replace=TRUE)
               samplecon0 <- sample(concentrations0,length(concentrations0),replace=TRUE)
               # create 80000 simulations from each replicated data:
               sampleser <- sample(sampleser,80000,replace=TRUE)
               samplecon <- sample(samplecon,80000,replace=TRUE)
               samplecon0 <- sample(samplecon0,80000,replace=TRUE)
               lines(ecdf(log(sampleser*samplecon*RM[foodindex[i],hazardindexM[h]])/log(10)),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="red")
               lines(ecdf(log(sampleser*samplecon0*RM[foodindex[i],hazardindexM[h]])/log(10)),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col="blue")
               }
             } # end of if logarithmic
             legend("bottomright",cex=0.6,paste(c("expo (2.5 %):","expo (50 %):","expo (97.5 %):"),round(quantile(100*PM[foodindex[i],hazardindexM[h]]*pM[,hazardindexM[h],foodindex[i]]*ppred[,foodindex[i]],c(0.025,0.5,0.975),names=FALSE),2),c("%","%","%")  ))
            } # end of if cumulative
         }} # end of nhusedM nfused
     } # end of if nhusedM nfused >0
     
     } # end of if theresults = "Exposures"
   })

   #########################################################
   output$distPlot4 = renderPlot({
     # generate results based on inputs from ui.R:  
     # uncertainties of variability quantiles
     
     currentresults()
     
     theresults = input$selectresults
     nV = input$nV  # variability iterations for simulations
     nU = input$nU  # uncertainty iterations for simulations
     Rall = input$factor  # processing factors
     Pall = input$pfactor  # prevalence factors
     
     if(!is.element("None",input$selectQ)){
       
          # which percentile is selected?: 
          if(is.element("Q5% Exposure",input$selectQ)){theQ=5}
          if(is.element("Q10% Exposure",input$selectQ)){theQ=10}
          if(is.element("Q25% Exposure",input$selectQ)){theQ=25}
          if(is.element("Q50% Exposure",input$selectQ)){theQ=50}
          if(is.element("Q75% Exposure",input$selectQ)){theQ=75}
          if(is.element("Q90% Exposure",input$selectQ)){theQ=90}
          if(is.element("Q95% Exposure",input$selectQ)){theQ=95}
       
           
     foodnamesused = input$thefoodnames
     nfused = length(foodnamesused)
     foodindex = match(foodnamesused,foodnames)
     
     hazardnamesused = input$thehazardnames
     hazardtypesused = hazardtypes[is.element(hazardnames,hazardnamesused)]
     nhused = length(hazardnamesused)
     hazardnamesK = hazardnames[hazardtypes=="K"]
     hazardnamesM = hazardnames[hazardtypes=="M"]
     hazardnamesusedK = hazardnamesused[hazardtypesused=="K"]
     hazardnamesusedM = hazardnamesused[hazardtypesused=="M"]
     nhusedK = length(hazardnamesusedK)
     nhusedM = length(hazardnamesusedM)
     hazardindex = match(hazardnamesused,hazardnames)
     hazardindexK = match(hazardnamesusedK,hazardnamesK)
     hazardindexM = match(hazardnamesusedM,hazardnamesM)
     
     if( (nhused>0)&(nfused>0) ){
       par(mfrow=c(nhused,1),cex.lab=1.3,cex.main=1.3,yaxt="n")
       }
     
    # generate nU variability distributions (each with nV variability simulations), 
    # evaluate quantiles for each of those variability distributions:
     mc <- round(seq(1,input$Iterations-burnin,length=nU),0)
     
     if((nhusedK>0)&(nfused>0)){ # Chemical exposure quantiles
       # redefine dimensions if scalar
       if(nf==1){
         mus0 <- matrix(mus0,input$Iterations-burnin,1)
         sigs <- matrix(sigs,input$Iterations-burnin,1)
         ppred <- matrix(ppred,input$Iterations-burnin,1)
         mucK <- array(mucK,dim=c(input$Iterations-burnin,1,1))
         sigcK <- array(sigcK,dim=c(input$Iterations-burnin,1,1))
         pK <- array(pK,dim=c(input$Iterations-burnin,1,1))
         Ts0 <- array(Ts0,dim=c(input$Iterations-burnin,1,1)) 
         Ts <- array(Ts,dim=c(input$Iterations-burnin,1,1)) 
       } 
       RK = matrix(NA,nf,nhK) # concentration factors
       RK[1:nf,1:nhK] = Rall[1:nf,is.element(hazardnames,hazardnamesusedK)]
       logRK = log(RK)
       PK = matrix(NA,nf,nhK) # prevalence factors
       PK[1:nf,1:nhK] = Pall[1:nf,is.element(hazardnames,hazardnamesusedK)]
       
     if(input$modelchoice == "Independent days"){
       if(nf==1){ # redefine dimensions if scalars
       logitp0 <- matrix(logitp0,input$Iterations-burnin,1)
       Tp <- array(Tp,dim=c(input$Iterations-burnin,1,1))
       }
       if(input$modelchoice2 == "No"){ between.user.pvar <- 0 } 
       if(input$modelchoice2 == "Yes"){ between.user.pvar <- 1 }
       logitpmc <- matrix(NA,nV,nf); pmc<-logitpmc; musmc <- pmc
       Eemc <- array(NA,dim=c(nV,nhusedK,nfused))  # for all days
       Eemcconuse <- array(NA,dim=c(nV,nhusedK,nfused)) # for contaminated consumption days
       Eetotmc <- matrix(NA,nV,nhusedK)
       Eetotmcconuse <- matrix(NA,nV,nhusedK)
       Q <- matrix(NA,nU,nhusedK) # for chronic exposure all days
       Qplus <- matrix(NA,nU,nhusedK) # for chronic exposure from consumption days
       thin <- 0 # for indexing a thinned sample of (simulated) variability distributions
       exposurevarsample <- matrix(NA,ceiling(nU/5),nV) # for thinned uncertainty sample
       for(u in 1:nU){ # for nU parameter sets
         for(v in 1:nV){ # for nV variable values per each parameter set
           if(nf>1){ # if many foods
             if(between.user.pvar==1){ # variability between users' frequencies    
               logitpmc[v,1:nf] <- rmvnorm(1,logitp0[mc[u],1:nf],solve(Tp[mc[u],1:nf,1:nf]))
             }
             if(between.user.pvar==0){ # no variability between users' frequencies   
               logitpmc[v,1:nf] <- logitp0[mc[u],1:nf]
             }
           pmc[v,1:nf] <- exp(logitpmc[v,1:nf])/(1+exp(logitpmc[v,1:nf])) # individual use probability
           musmc[v,1:nf] <- rmvnorm(1,mus0[mc[u],1:nf],solve(Ts0[mc[u],1:nf,1:nf])) # individual mean amount
           }
           if(nf==1){ # if only one food
             if(between.user.pvar==1){ # variability between users' frequencies  
               logitpmc[v,1] <- rnorm(1,logitp0[mc[u],1],sqrt(1/Tp[mc[u],1,1]))
             }
             if(between.user.pvar==0){ # no variability between users' frequencies  
               logitpmc[v,1] <- logitp0[mc[u],1]
             }
             pmc[v,1] <- exp(logitpmc[v,1])/(1+exp(logitpmc[v,1])) # individual use probability
             musmc[v,1] <- rnorm(1,mus0[mc[u],1],sqrt(1/Ts0[mc[u],1,1])) # individual mean amount  
           }
           
           for(h in 1:nhusedK){
             for(i in 1:nfused){
               Eemc[v,h,i]<-pK[mc[u],hazardindexK[h],foodindex[i]]*PK[foodindex[i],hazardindexK[h]]*
                 pmc[v,foodindex[i]]*exp(logRK[foodindex[i],hazardindexK[h]]
                                          +musmc[v,foodindex[i]]
                                          +0.5*sigs[mc[u],foodindex[i]]^2
                                          +mucK[mc[u],hazardindexK[h],foodindex[i]]
                                          +0.5*sigcK[mc[u],hazardindexK[h],foodindex[i]]^2)
               
               Eemcconuse[v,h,i] <- exp(logRK[foodindex[i],hazardindexK[h]]
                                     +musmc[v,foodindex[i]]
                                     +0.5*sigs[mc[u],foodindex[i]]^2
                                     +mucK[mc[u],hazardindexK[h],foodindex[i]]
                                     +0.5*sigcK[mc[u],hazardindexK[h],foodindex[i]]^2) 
             }
             # simulated total chronic exposure for individual:
             Eetotmc[v,h] <- sum(Eemc[v,h,1:nfused]) 
             # simulated total chronic exposure for individual for contaminated consumption days 
             Eetotmcconuse[v,h] <- sum(Eemcconuse[v,h,1:nfused])   
            }
         } # end of v (variability) 
         for(h in 1:nhusedK){
           Q[u,h]<-quantile(Eetotmc[,h],theQ/100,names=FALSE)
           Qplus[u,h]<-quantile(Eetotmcconuse[,h],theQ/100,names=FALSE)
           }
         #######################################################
         # pick out thinned sample:
         if(ceiling(u/5)==floor(u/5)){ thin<-thin+1; exposurevarsample[thin,1:nV]<- t(Eetotmcconuse[1:nV,h]) }
         #######################################################
       } # end of u (uncertainty)
       for(h in 1:nhusedK){
       if(input$selectscale=="Absolute"){
       #########################################
       plot(ecdf(exposurevarsample[1,1:nV]),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(min(exposurevarsample[1:thin,1:nV]),max(exposurevarsample[1:thin,1:nV])),lwd=1,lty=3,col=rgb(0.816,0.004,0.435),xlab="C.exposure/bw+",ylab="",main=paste("Exposure:",hazardnamesusedK[h],"total from",nfused,"foods (chronic)"))
       for(a in 2:thin){
       lines(ecdf(exposurevarsample[a,1:nV]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
       }
         quplim <- quantile(Qplus[,h],0.95,names=FALSE)
         qlolim <- quantile(Qplus[,h],0.05,names=FALSE)
       lines(density(Qplus[,h],from=qlolim,to=quplim)$x,density(Qplus[,h],from=qlolim,to=quplim)$y/max(density(Qplus[,h],from=qlolim,to=quplim)$y),lwd=3)   
       lines(quantile(Qplus[,h],c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
       lines(quantile(Qplus[,h],c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
       legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),c(round(quantile(Q[,h],c(0.05,0.5,0.95),names=FALSE),2)) )) 
       legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(Qplus[,h],c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
       }
       if(input$selectscale=="Logarithmic"){
       plot(ecdf(log(exposurevarsample[1,1:nV])/log(10)),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(min(log(exposurevarsample[1:thin,1:nV])/log(10)),max(log(exposurevarsample[1:thin,1:nV])/log(10))),lwd=1,lty=3,col=raspberry,ylab="",xlab="log (C.exposure/bw+)",main=paste("Exposure:",hazardnamesusedK[h],"total from",nfused,"foods (chronic)"))
       for(a in 2:thin){
       lines(ecdf(log(exposurevarsample[a,1:nV])/log(10)),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
       }
         quplim <- quantile(log10(Qplus[,h]),0.95,names=FALSE)
         qlolim <- quantile(log10(Qplus[,h]),0.05,names=FALSE)
       lines(density(log10(Qplus[,h]),from=qlolim,to=quplim)$x,density(log10(Qplus[,h]),from=qlolim,to=quplim)$y/max(density(log10(Qplus[,h]),from=qlolim,to=quplim)$y),lwd=3)   
       lines(quantile(log10(Qplus[,h]),c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
       lines(quantile(log10(Qplus[,h]),c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
       legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),c(round(quantile(log10(Q[,h]),c(0.05,0.5,0.95),names=FALSE),2)) )) 
       legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log10(Qplus[,h]),c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
       }    
       } # uncertainty distribution of Q% exposure for hth chemical
       
     } # end of if independent days
       
     if(input$modelchoice == "Dependent days"){
       musmc <- matrix(NA,nV,nf)
       Eemc <- array(NA,dim=c(nV,nhusedK,nfused)) # for all days
       Eemcconuse <- array(NA,dim=c(nV,nhusedK,nfused)) # for contaminated consumption days
       Eetotmc <- matrix(NA,nV,nhusedK)
       Eetotmcconuse <- matrix(NA,nV,nhusedK)
       Q <- matrix(NA,nU,nhusedK) # for chronic exposure all days
       Qplus <- matrix(NA,nU,nhusedK) # for chronic exposure from consumption days
       thin <- 0 # for indexing a thinned sample of (simulated) variability distributions
       exposurevarsample <- matrix(NA,ceiling(nU/5),nV) # for thinned uncertainty sample
       for(u in 1:nU){ # for nU parameter sets
         for(v in 1:nV){ # for nV variable values per each parameter set
           if(nf>1){ # if several foods
           musmc[v,1:nf] <- rmvnorm(1,mus0[mc[u],1:nf],solve(Ts0[mc[u],1:nf,1:nf])) # individual mean amount
           }
           if(nf==1){ # if only one food
           musmc[v,1] <- rnorm(1,mus0[mc[u],1],sqrt(1/Ts0[mc[u],1,1])) # individual mean amount    
           }
           for(h in 1:nhusedK){
             for(i in 1:nfused){
           Eemc[v,h,i] <- pK[mc[u],hazardindexK[h],foodindex[i]]*PK[foodindex[i],hazardindexK[h]]*
                          ppred[mc[u],foodindex[i]]*exp(logRK[foodindex[i],hazardindexK[h]]
                                              +musmc[v,foodindex[i]]
                                              +0.5*sigs[mc[u],foodindex[i]]^2
                                              +mucK[mc[u],hazardindexK[h],foodindex[i]]
                                              +0.5*sigcK[mc[u],hazardindexK[h],foodindex[i]]^2)
           Eemcconuse[v,h,i] <- exp(logRK[foodindex[i],hazardindexK[h]]
                                 +musmc[v,foodindex[i]]
                                 +0.5*sigs[mc[u],foodindex[i]]^2
                                 +mucK[mc[u],hazardindexK[h],foodindex[i]]
                                 +0.5*sigcK[mc[u],hazardindexK[h],foodindex[i]]^2)
             } # end of i
             # simulated total chronic exposure for individual: 
             Eetotmc[v,h] <- sum(Eemc[v,h,1:nfused])    
             # simulated total chronic exposure for individual for contaminated consumption days:
             Eetotmcconuse[v,h] <- sum(Eemcconuse[v,h,1:nfused]) 
           } # end of h 
         } # end of v (variability)   
         for(h in 1:nhusedK){
           Q[u,h]<-quantile(Eetotmc[,h],theQ/100,names=FALSE)
           Qplus[u,h]<-quantile(Eetotmcconuse[,h],theQ/100,names=FALSE)
           }
         #######################################################
         # pick out thinned sample:
         if(ceiling(u/5)==floor(u/5)){ thin<-thin+1; exposurevarsample[thin,1:nV]<- t(Eetotmcconuse[1:nV,h]) }
         #######################################################
       } # end of u (uncertainty)
       for(h in 1:nhusedK){
         if(input$selectscale=="Absolute"){
         plot(ecdf(exposurevarsample[1,1:nV]),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(min(exposurevarsample[1:thin,1:nV]),max(exposurevarsample[1:thin,1:nV])),lwd=1,lty=3,col=raspberry,xlab="C.exposure/bw+",ylab="",main=paste("Uncertainty of distribution:",hazardnamesusedK[h],"total from",nfused,"foods (chronic)"))
         for(a in 2:thin){
            lines(ecdf(exposurevarsample[a,1:nV]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
         }
           quplim <- quantile(Qplus[,h],0.95,names=FALSE)
           qlolim <- quantile(Qplus[,h],0.05,names=FALSE)
           lines(density(Qplus[,h],from=qlolim,to=quplim)$x,density(Qplus[,h],from=qlolim,to=quplim)$y/max(density(Qplus[,h],from=qlolim,to=quplim)$y),lwd=3)   
           lines(quantile(Qplus[,h],c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
           lines(quantile(Qplus[,h],c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
         legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5%","50%","95%"),round(quantile(Q[,h],c(0.05,0.5,0.95),names=FALSE),2) ))
         legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5%","50%","95%"),round(quantile(Qplus[,h],c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
         }
         if(input$selectscale=="Logarithmic"){
         plot(ecdf(log(exposurevarsample[1,1:nV])/log(10)),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(min(log(exposurevarsample[1:thin,1:nV])/log(10)),max(log(exposurevarsample[1:thin,1:nV])/log(10))),lwd=1,lty=3,col=raspberry,ylab="",xlab="log (C.exposure/bw+)",main=paste("Uncertainty of distribution:",hazardnamesusedK[h],"total from",nfused,"foods (chronic)"))
         for(a in 2:thin){
          lines(ecdf(log(exposurevarsample[a,1:nV])/log(10)),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
         }    
           quplim <- quantile(log10(Qplus[,h]),0.95,names=FALSE)
           qlolim <- quantile(log10(Qplus[,h]),0.05,names=FALSE)
           lines(density(log10(Qplus[,h]),from=qlolim,to=quplim)$x,density(log10(Qplus[,h]),from=qlolim,to=quplim)$y/max(density(log10(Qplus[,h]),from=qlolim,to=quplim)$y),lwd=3)  
           lines(quantile(log10(Qplus[,h]),c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
           lines(quantile(log10(Qplus[,h]),c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
         legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log(Q[,h])/log(10),c(0.05,0.5,0.95),names=FALSE),2) ))
         legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log10(Qplus[,h]),c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
         }
       } # end of uncertainty distribution of Q% exposure for hth chemical
       
     } # end of if dependent days
     } # end of if nhusedK>0 nfused>0
     
     if((nhusedM>0)&(nfused>0)){  # Microbial exposure quantiles
       # redefine dimensions if scalar
       if(nf==1){
         mus0 <- matrix(mus0,input$Iterations-burnin,1)
         sigs <- matrix(sigs,input$Iterations-burnin,1)
         ppred <- matrix(ppred,input$Iterations-burnin,1)
         mucM <- array(mucM,dim=c(input$Iterations-burnin,1,1))
         sigcM <- array(sigcM,dim=c(input$Iterations-burnin,1,1)) 
         pM <- array(pM,dim=c(input$Iterations-burnin,1,1))
         Ts0 <- array(Ts0,dim=c(input$Iterations-burnin,1,1)) 
         Ts <- array(Ts,dim=c(input$Iterations-burnin,1,1)) 
       } 
       RM = matrix(NA,nf,nhM) # factors for concentrations
       RM[1:nf,1:nhM] = Rall[1:nf,is.element(hazardnames,hazardnamesusedM)]
       logRM = log(RM)
       PM = matrix(NA,nf,nhM) # factors for prevalence
       PM[1:nf,1:nhM] = Pall[1:nf,is.element(hazardnames,hazardnamesusedM)]
       
       if(input$modelchoice == "Independent days"){
         if(nf==1){ # redefine dimensions if scalar
         logitp0 <- matrix(logitp0,input$Iterations-burnin,1)     
         Tp <- array(Tp,dim=c(input$Iterations-burnin,1,1))
         }
         if(input$modelchoice2 == "No"){ between.user.pvar <- 0 } 
         if(input$modelchoice2 == "Yes"){ between.user.pvar <- 1 }
         wmc <- numeric()
         logitpmc <- matrix(NA,nV,nf); pmc<-logitpmc; 
         musmc <- matrix(NA,nV,nf)
         Umc<-musmc; smc<-musmc
         Imc <- array(NA,dim=c(nV,nhM,nf)); cmc<-Imc
         Eemc <- array(NA,dim=c(nV,nhusedM,nfused))
         Eemcconuse <- Eemc 
         Eetotmc <- matrix(NA,nV,nhusedM);Eetotmcplus <- matrix(NA,nV,nhusedM); nplus<-matrix(NA,nU,nhusedM)
         Eetotmczeros <- matrix(NA,nV,nhusedM)
         Q <- matrix(NA,nU,nhusedM); Qplus <- matrix(NA,nU,nhusedM) 
         thin <- 0 # for indexing a thinned sample of (simulated) variability distributions
         exposurevarsample <- matrix(NA,ceiling(nU/5),nV) # for thinned uncertainty sample
         for(u in 1:nU){ # for nU parameter sets
           for(v in 1:nV){ # for nV variable values per each parameter set
             wmc[v] <- rlnorm(1,muw,sigw) # bodyweight for v:th individual
             if(nf>1){ # if many foods
             if(between.user.pvar==1){ # variability between users' frequencies    
             logitpmc[v,1:nf] <- rmvnorm(1,logitp0[mc[u],1:nf],solve(Tp[mc[u],1:nf,1:nf]))
             }
             if(between.user.pvar==0){ # no variability between users' frequencies   
             logitpmc[v,1:nf] <- logitp0[mc[u],1:nf]
             }    
             pmc[v,1:nf] <- exp(logitpmc[v,1:nf])/(1+exp(logitpmc[v,1:nf]))
             Umc[v,1:nf] <- rbinom(nf,rep(1,nf),pmc[v,1:nf]) # actual use
             musmc[v,1:nf] <- rmvnorm(1,mus0[mc[u],1:nf],solve(Ts0[mc[u],1:nf,1:nf]))
             smc[v,1:nf] <- exp(rmvnorm(1,musmc[v,1:nf],solve(Ts[mc[u],1:nf,1:nf])))  # actual amount
             }
             if(nf==1){ # if only one food
               if(between.user.pvar==1){ # variability between users' frequencies  
               logitpmc[v,1] <- rnorm(1,logitp0[mc[u],1],sqrt(1/Tp[mc[u],1,1]))
               }
               if(between.user.pvar==0){ # no variability between users' frequencies  
               logitpmc[v,1] <- logitp0[mc[u],1]
               }
               pmc[v,1] <- exp(logitpmc[v,1])/(1+exp(logitpmc[v,1]))
               Umc[v,1] <- rbinom(1,rep(1,1),pmc[v,1]) # actual use
               musmc[v,1] <- rnorm(1,mus0[mc[u],1],sqrt(1/Ts0[mc[u],1,1]))  
               smc[v,1] <- exp(rnorm(1,musmc[v,1],sqrt(1/Ts[mc[u],1,1])))  # actual amount   
             }
             for(h in 1:nhM){
             Imc[v,h,1:nf] <- rbinom(nf,rep(1,nf),pM[v,h,1:nf]*PM[1:nf,h]) # actual contamination yes/no
             cmc[v,h,1:nf] <- rlnorm(nf,mucM[mc[u],h,1:nf],sigcM[mc[u],h,1:nf]) # actual contamination level
             } 
             for(h in 1:nhusedM){
                 for(i in 1:nfused){
                 Eemc[v,h,i] <- Imc[v,hazardindexM[h],foodindex[i]]*Umc[v,foodindex[i]]*
                                smc[v,foodindex[i]]*
                                RM[foodindex[i],hazardindexM[h]]*
                                cmc[v,hazardindexM[h],foodindex[i]]*wmc[v] 
                 Eemcconuse[v,h,i] <- smc[v,foodindex[i]]*
                                RM[foodindex[i],hazardindexM[h]]*
                                cmc[v,hazardindexM[h],foodindex[i]]*wmc[v]
                 }
               if(sum(Eemc[v,h,1:nfused])<=100000){
                 Eetotmczeros[v,h] <- rpois(1,sum(Eemc[v,h,1:nfused]))}  # sum of all serving exposures, incl. zeros  
               if(sum(Eemc[v,h,1:nfused])>100000){
                 Eetotmczeros[v,h] <- round(rnorm(1,sum(Eemc[v,h,1:nfused]),sqrt(sum(Eemc[v,h,1:nfused]))))}  # sum of all serving exposures, incl. zeros  
               if(sum(Eemcconuse[v,h,1:nfused])<=100000){
               Eetotmc[v,h] <- rpois(1,sum(Eemcconuse[v,h,1:nfused]))} # simulated total acute exposure for individual when used and contaminated 
               if(sum(Eemcconuse[v,h,1:nfused])>100000){
               Eetotmc[v,h] <- round(rnorm(1,sum(Eemcconuse[v,h,1:nfused]),sqrt(sum(Eemcconuse[v,h,1:nfused]))))} # simulated total acute exposure for individual when used and contaminated    
             } # end of h
           } # end of v (variability) 
           for(h in 1:nhusedM){
             nplus[u,h] <- sum(Eetotmc[1:nV,h]>0)
             if(nplus[u,h]==0){  
               Eetotmcplus[1,h] <- NA
               Qplus[u,h] <- NA  # quantile from positive servings
             }
             if(nplus[u,h]>0){
               Eetotmcplus[1:nplus[u,h],h] <- Eetotmc[Eetotmc[1:nV,h]>0,h]
               Qplus[u,h]<-quantile(Eetotmcplus[1:nplus[u,h],h],theQ/100,names=FALSE) # quantile from pos servings
             }
             Q[u,h] <- quantile(Eetotmczeros[,h],theQ/100,names=FALSE) # quantile from ALL servings
           }
           #######################################################
           # pick out thinned sample of positive acute exposures:
           if(ceiling(u/5)==floor(u/5)){ 
             thin<-thin+1 
             if(nplus[u,h]==0){
               exposurevarsample[thin,1]<-NA  
             }
             if(nplus[u,h]>0){
               exposurevarsample[thin,1:nplus[u,h]] <- t(Eetotmcplus[1:nplus[u,h],h]) 
             }
             }
           #######################################################
         } # end of u (uncertainty)
         
         for(h in 1:nhusedM){
           if(input$selectscale=="Absolute"){
             if(sum(!is.na(exposurevarsample[1,]))>0){
             plot(ecdf(exposurevarsample[1,!is.na(exposurevarsample[1,])]),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(min(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE),0.1*max(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE)),lwd=1,lty=3,col=raspberry,xlab="A.dose+",ylab="",main=paste("Exposure:",hazardnamesusedM[h],"total from",nfused,"foods (acute)"))
             }
             for(a in 2:thin){
             if(sum(!is.na(exposurevarsample[a,]))>0){   
             lines(ecdf(exposurevarsample[a,!is.na(exposurevarsample[a,])]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
             }
             }
             quplim <- quantile(Qplus[,h],0.95,names=FALSE,na.rm=TRUE)
             qlolim <- quantile(Qplus[,h],0.05,names=FALSE,na.rm=TRUE)
             lines(density(Qplus[,h],na.rm=TRUE,from=qlolim,to=quplim)$x,density(Qplus[,h],na.rm=TRUE,from=qlolim,to=quplim)$y/max(density(Qplus[,h],na.rm=TRUE,from=qlolim,to=quplim)$y),lwd=3)   
             lines(quantile(Qplus[,h],c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             lines(quantile(Qplus[,h],c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
           legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(Q[,h],c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
           legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(Qplus[,h],c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
           }
           
           if(input$selectscale=="Logarithmic"){
             if(sum(!is.na(exposurevarsample[1,]))>0){
               plot(ecdf(log10(exposurevarsample[1,!is.na(exposurevarsample[1,])])),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(log10(min(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE)),log10(max(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE))),lwd=1,lty=3,col=raspberry,xlab="log A.dose+",ylab="",main=paste("Exposure:",hazardnamesusedM[h],"total from",nfused,"foods (acute)"))
             }
             for(a in 2:thin){
               if(sum(!is.na(exposurevarsample[a,]))>0){   
                 lines(ecdf(log10(exposurevarsample[a,!is.na(exposurevarsample[a,])])),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
               }
             }
             quplim <- quantile(log10(Qplus[,h]),0.95,names=FALSE,na.rm=TRUE)
             qlolim <- quantile(log10(Qplus[,h]),0.05,names=FALSE,na.rm=TRUE)
             lines(density(log10(Qplus[,h]),na.rm=TRUE,from=qlolim,to=quplim)$x,density(log10(Qplus[,h]),na.rm=TRUE,from=qlolim,to=quplim)$y/max(density(log10(Qplus[,h]),na.rm=TRUE,from=qlolim,to=quplim)$y),lwd=3)
             lines(quantile(log10(Qplus[,h]),c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             lines(quantile(log10(Qplus[,h]),c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log10(Q[,h]),c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
             legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log10(Qplus[,h]),c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
           }
           } # end of uncertainty distribution of Q% exposure for hth microbe
       } # end of if independent days
       
       if(input$modelchoice == "Dependent days"){
         wmc <- numeric()
         pmc <- matrix(NA,nV,nf) 
         musmc <- matrix(NA,nV,nf)
         Umc<-musmc; smc<-musmc
         Imc <- array(NA,dim=c(nV,nhM,nf)); cmc<-Imc
         Eemc <- array(NA,dim=c(nV,nhusedM,nfused))
         Eemcconuse <- Eemc 
         Eetotmc <- matrix(NA,nV,nhusedM);Eetotmcplus <- matrix(NA,nV,nhusedM); nplus<-matrix(NA,nU,nhusedM)
         Eetotmczeros <- matrix(NA,nV,nhusedM)
         Q <- matrix(NA,nU,nhusedM); Qplus <- matrix(NA,nU,nhusedM) 
         thin <- 0 # for indexing a thinned sample of (simulated) variability distributions
         exposurevarsample <- matrix(NA,ceiling(nU/5),nV) # for thinned uncertainty sample
         for(u in 1:nU){ # for nU parameter sets
           for(v in 1:nV){ # for nV variable values per each parameter set
             wmc[v] <- rlnorm(1,muw,sigw) # bodyweight
             Umc[v,1:nf] <- rbinom(nf,rep(1,nf),ppred[mc[u],1:nf]) # actual use
             if(nf>1){ # if many foods
             musmc[v,1:nf] <- rmvnorm(1,mus0[mc[u],1:nf],solve(Ts0[mc[u],1:nf,1:nf]))
             smc[v,1:nf] <- exp(rmvnorm(1,musmc[v,1:nf],solve(Ts[mc[u],1:nf,1:nf]))) # actual amount
             }
             if(nf==1){ # if only one food
             musmc[v,1] <- rnorm(1,mus0[mc[u],1],sqrt(1/Ts0[mc[u],1,1]))
             smc[v,1] <- exp(rnorm(1,musmc[v,1],sqrt(1/Ts[mc[u],1,1]))) # actual amount  
             }
             for(h in 1:nhM){
               Imc[v,h,1:nf] <- rbinom(nf,rep(1,nf),pM[v,h,1:nf]*PM[1:nf,h]) # actual contamination yes/no
               cmc[v,h,1:nf] <- rlnorm(nf,mucM[mc[u],h,1:nf],sigcM[mc[u],h,1:nf]) # actual contamination level
             } 
             for(h in 1:nhusedM){
                 for(i in 1:nfused){
                 Eemc[v,h,i] <- Imc[v,hazardindexM[h],foodindex[i]]*Umc[v,foodindex[i]]*
                               smc[v,foodindex[i]]*
                               RM[foodindex[i],hazardindexM[h]]*
                               cmc[v,hazardindexM[h],foodindex[i]]*wmc[v] 
                 Eemcconuse[v,h,i] <- smc[v,foodindex[i]]*
                               RM[foodindex[i],hazardindexM[h]]*
                               cmc[v,hazardindexM[h],foodindex[i]]*wmc[v] 
                 }
               if(sum(Eemc[v,h,1:nfused])<=100000){
                 Eetotmczeros[v,h] <- rpois(1,sum(Eemc[v,h,1:nfused]))}  # sum of all serving exposures, incl. zeros 
               if(sum(Eemc[v,h,1:nfused])>100000){
                 Eetotmczeros[v,h] <- rnorm(1,sum(Eemc[v,h,1:nfused]),sqrt(sum(Eemc[v,h,1:nfused])))}
               if(sum(Eemcconuse[v,h,1:nfused])<=100000){
                 Eetotmc[v,h] <- rpois(1,sum(Eemcconuse[v,h,1:nfused]))} # simulated total acute exposure for individual when used and contaminated 
               if(sum(Eemcconuse[v,h,1:nfused])>100000){
                 Eetotmc[v,h] <- round(rnorm(1,sum(Eemcconuse[v,h,1:nfused]),sqrt(sum(Eemcconuse[v,h,1:nfused]))))} # simulated total acute exposure for individual when used and contaminated
             }  # end of h
           } # end of v (variability)  
           for(h in 1:nhusedM){
             nplus[u,h] <- sum(Eetotmc[1:nV,h]>0)
             if(nplus[u,h]==0){  
               Eetotmcplus[1,h] <- NA
               Qplus[u,h] <- NA  # quantile from positive servings
             }
             if(nplus[u,h]>0){
               Eetotmcplus[1:nplus[u,h],h] <- Eetotmc[Eetotmc[1:nV,h]>0,h]
               Qplus[u,h]<-quantile(Eetotmcplus[1:nplus[u,h],h],theQ/100,names=FALSE) # quantile from pos servings
             }
             Q[u,h] <- quantile(Eetotmczeros[,h],theQ/100,names=FALSE) # quantile from ALL servings
             }
           #######################################################
           # pick out thinned sample of positive acute exposures:
           if(ceiling(u/5)==floor(u/5)){ 
             thin<-thin+1 
             if(nplus[u,h]==0){
               exposurevarsample[thin,1]<-NA  
             }
             if(nplus[u,h]>0){
               exposurevarsample[thin,1:nplus[u,h]] <- t(Eetotmcplus[1:nplus[u,h],h]) 
             }
           }
           #######################################################
         } # end of u (uncertainty)
         for(h in 1:nhusedM){
           if(input$selectscale=="Absolute"){
             if(sum(!is.na(exposurevarsample[1,]))>0){
               plot(ecdf(exposurevarsample[1,!is.na(exposurevarsample[1,])]),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(min(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE),0.1*max(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE)),lwd=1,lty=3,col=raspberry,xlab="A.dose+",ylab="",main=paste("Exposure:",hazardnamesusedM[h],"total from",nfused,"foods (acute)"))
             }
             for(a in 2:thin){
               if(sum(!is.na(exposurevarsample[a,]))>0){   
                 lines(ecdf(exposurevarsample[a,!is.na(exposurevarsample[a,])]),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)     
               }
             }
             quplim <- quantile(Qplus[,h],0.95,names=FALSE,na.rm=TRUE)
             qlolim <- quantile(Qplus[,h],0.05,names=FALSE,na.rm=TRUE)
             lines(density(Qplus[,h],na.rm=TRUE,from=qlolim,to=quplim)$x,density(Qplus[,h],na.rm=TRUE,from=qlolim,to=quplim)$y/max(density(Qplus[,h],na.rm=TRUE,from=qlolim,to=quplim)$y),lwd=3)   
             lines(quantile(Qplus[,h],c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             lines(quantile(Qplus[,h],c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(Q[,h],c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
             legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(Qplus[,h],c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
           }
           
           if(input$selectscale=="Logarithmic"){
             if(sum(!is.na(exposurevarsample[1,]))>0){
               plot(ecdf(log10(exposurevarsample[1,!is.na(exposurevarsample[1,])])),verticals=TRUE,do.points=FALSE,yaxt="s",xlim=c(log10(min(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE)),log10(max(exposurevarsample[1:thin,1:max(nplus)],na.rm=TRUE))),lwd=1,lty=3,col=raspberry,xlab="log A.dose+",ylab="",main=paste("Exposure:",hazardnamesusedM[h],"total from",nfused,"foods (acute)"))
             }
             for(a in 2:thin){
               if(sum(!is.na(exposurevarsample[a,]))>0){   
                 lines(ecdf(log10(exposurevarsample[a,!is.na(exposurevarsample[a,])])),verticals=TRUE,do.points=FALSE,lwd=1,lty=3,col=raspberry)    
               }
             }
             quplim <- quantile(log10(Qplus[,h]),0.95,names=FALSE,na.rm=TRUE)
             qlolim <- quantile(log10(Qplus[,h]),0.05,names=FALSE,na.rm=TRUE)
             lines(density(log10(Qplus[,h]),na.rm=TRUE,from=qlolim,to=quplim)$x,density(log10(Qplus[,h]),na.rm=TRUE,from=qlolim,to=quplim)$y/max(density(log10(Qplus[,h]),na.rm=TRUE,from=qlolim,to=quplim)$y),lwd=3)
             lines(quantile(log10(Qplus[,h]),c(0.05,0.05),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             lines(quantile(log10(Qplus[,h]),c(0.95,0.95),names=FALSE,na.rm=TRUE),c(0,1),lwd=3)
             legend("topright",title=paste("Total",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log10(Q[,h]),c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
             legend("bottomright",title=paste("Positives",theQ,"% quantile:"),paste(c("5 %:","50 %:","95 %:"),round(quantile(log10(Qplus[,h]),c(0.05,0.5,0.95),names=FALSE,na.rm=TRUE),2) ))
           }
           } # end of uncertainty distribution of Q% exposure for hth microbe
       } # end of if dependent days   
       
     } # end of if nhusedM>0 nfused>0
     
     } # end of if theresults
   })  # end of renderPlot
   ######################################################################
      
   output$distPlot5 <- renderPlot({      
     # generate results based on inputs from ui.R:  
     # MCMC diagnostic plots
     
     currentresults()
     
     thediagnostics = input$diagchoice
     if(sum(is.element(c("Concentration parameters"),thediagnostics))>0){
       
       foodnamesused = input$thefoodnames
       nfused = length(foodnamesused)
       foodindex = match(foodnamesused,foodnames)
       
       hazardnamesused = input$thehazardnames
       hazardtypesused = hazardtypes[is.element(hazardnames,hazardnamesused)]
       nhused = length(hazardnamesused)
       hazardnamesK = hazardnames[hazardtypes=="K"]
       hazardnamesM = hazardnames[hazardtypes=="M"]
       hazardnamesusedK = hazardnamesused[hazardtypesused=="K"]
       hazardnamesusedM = hazardnamesused[hazardtypesused=="M"]
       nhusedK = length(hazardnamesusedK)
       nhusedM = length(hazardnamesusedM)
       hazardindex = match(hazardnamesused,hazardnames)
       hazardindexK = match(hazardnamesusedK,hazardnamesK)
       hazardindexM = match(hazardnamesusedM,hazardnamesM)
       
       if((nhusedK>0)&(nhusedM==0) ){
         par(mar=rep(2,4),mfrow=c(nhusedK*3,nfused),cex.lab=1.3,cex.main=1.3)
       }
       if((nhusedM>0)&(nhusedK==0) ){
         par(mar=rep(2,4),mfrow=c(nhusedM*3,nfused),cex.lab=1.3,cex.main=1.3)
       }
       if((nhusedM>0)&(nhusedK>0) ){
         par(mar=rep(2,4),mfrow=c((nhusedK*3+nhusedM*3),nfused),cex.lab=1.3,cex.main=1.3)
       }
       
       if((nhusedK>0)&(nfused>0)){
         # redefine dimensions if scalars
         if((nhK==1)&(nf==1)){
           mucK <- array(mucK,dim=c(input$Iterations-burnin,1,1))
           sigcK <- array(sigcK,dim=c(input$Iterations-burnin,1,1))
           pK <- array(pK,dim=c(input$Iterations-burnin,1,1))
         } 
         for(h in 1:nhusedK){
           for(i in 1:nfused){
             plot(mucK[,hazardindexK[h],foodindex[i]]/log(10),pch=16,cex=0.5,main=bquote(.(hazardnamesusedK[h])~"in"~.(foodnamesused[i])~":"~mu),xlab="Iteration",col=raspberry)
             lines(0.3*(input$Iterations-burnin)/max(density(mucK[,hazardindexK[h],foodindex[i]]/log(10))$y)*density(mucK[,hazardindexK[h],foodindex[i]]/log(10))$y,density(mucK[,hazardindexK[h],foodindex[i]]/log(10))$x,lwd=3)
             plot(sigcK[,hazardindexK[h],foodindex[i]]/log(10),pch=16,cex=0.5,main=bquote(.(hazardnamesusedK[h])~"in"~.(foodnamesused[i])~":"~sigma),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(sigcK[,hazardindexK[h],foodindex[i]]/log(10))$y)*density(sigcK[,hazardindexK[h],foodindex[i]]/log(10))$y,density(sigcK[,hazardindexK[h],foodindex[i]]/log(10))$x,lwd=3)
             plot(pK[,hazardindexK[h],foodindex[i]],pch=16,cex=0.5,main=bquote(.(hazardnamesusedK[h])~"in"~.(foodnamesused[i])~":"~q),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(pK[,hazardindexK[h],foodindex[i]])$y)*density(pK[,hazardindexK[h],foodindex[i]])$y,density(pK[,hazardindexK[h],foodindex[i]])$x,lwd=3)
          }} # for, for
       } # if
       if( (nhusedM>0)&(nfused>0) ){
         # redefine dimensions if scalars
         if((nhM==1)&(nf==1)){
           mucM <- array(mucM,dim=c(input$Iterations-burnin,1,1))
           sigcM <- array(sigcM,dim=c(input$Iterations-burnin,1,1))
           pM <- array(pM,dim=c(input$Iterations-burnin,1,1))
         } 
          for(h in 1:nhusedM){
           for(i in 1:nfused){
             plot(mucM[,hazardindexM[h],foodindex[i]]/log(10),pch=16,cex=0.5,main=bquote(.(hazardnamesusedM[h])~"in"~.(foodnamesused[i])~":"~mu),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(mucM[,hazardindexM[h],foodindex[i]]/log(10))$y)*density(mucM[,hazardindexM[h],foodindex[i]]/log(10))$y,density(mucM[,hazardindexM[h],foodindex[i]]/log(10))$x,lwd=3)
             plot(sigcM[,hazardindexM[h],foodindex[i]]/log(10),pch=16,cex=0.5,main=bquote(.(hazardnamesusedM[h])~"in"~.(foodnamesused[i])~":"~sigma),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(sigcM[,hazardindexM[h],foodindex[i]]/log(10))$y)*density(sigcM[,hazardindexM[h],foodindex[i]]/log(10))$y,density(sigcM[,hazardindexM[h],foodindex[i]]/log(10))$x,lwd=3)
             plot(pM[,hazardindexM[h],foodindex[i]],pch=16,cex=0.5,main=bquote(.(hazardnamesusedM[h])~"in"~.(foodnamesused[i])~":"~q),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(pM[,hazardindexM[h],foodindex[i]])$y)*density(pM[,hazardindexM[h],foodindex[i]])$y,density(pM[,hazardindexM[h],foodindex[i]])$x,lwd=3)
           }} # for, for
       } # if
       
       } # if concentration parameters
     
     if(is.element("Consumption parameters",thediagnostics)){
       
       foodnamesused = input$thefoodnames
       nfused = length(foodnamesused)
       foodindex = match(foodnamesused,foodnames)
       
       hazardnamesused = input$thehazardnames
       hazardtypesused = hazardtypes[is.element(hazardnames,hazardnamesused)]
       nhused = length(hazardnamesused)
       
       
       if( (nhused>0)&(nfused>0) ){
         par(mar=rep(2,4),mfrow=c(3,nfused),cex.lab=1.3,cex.main=1.3)
       }
       
       if(nfused>0){
         # redefine dimensions if scalar
         if(nf==1){
           mus0 <- matrix(mus0,input$Iterations-burnin,1)
           sigs <- matrix(sigs,input$Iterations-burnin,1)
           ppred <- matrix(ppred,input$Iterations-burnin,1)
         } 
           for(i in 1:nfused){
             plot(mus0[,foodindex[i]]/log(10),pch=16,cex=0.5,main=bquote(.(foodnamesused[i])~":"~mu),xlab="Iteration",col=raspberry)  
             lines(0.3*(input$Iterations-burnin)/max(density(mus0[,foodindex[i]]/log(10))$y)*density(mus0[,foodindex[i]]/log(10))$y,density(mus0[,foodindex[i]]/log(10))$x,lwd=3)
             plot(sigs[,foodindex[i]]/log(10),pch=16,cex=0.5,main=bquote(.(foodnamesused[i])~":"~sigma),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(sigs[,foodindex[i]]/log(10))$y)*density(sigs[,foodindex[i]]/log(10))$y,density(sigs[,foodindex[i]]/log(10))$x,lwd=3)
             plot(ppred[,foodindex[i]],pch=16,cex=0.5,main=bquote(.(foodnamesused[i])~":"~p),xlab="Iteration",col=raspberry) 
             lines(0.3*(input$Iterations-burnin)/max(density(ppred[,foodindex[i]])$y)*density(ppred[,foodindex[i]])$y,density(ppred[,foodindex[i]])$x,lwd=3)
           } # for
       } # if
       
     } # if consumption parameters 
     
     })
   
   
   ##################################################
   
   output$distPlot6 <- renderPlot({      
     # generate results based on inputs from ui.R: 
     # Correlation plots for consumptions
     
     currentresults()
    
     theresults = input$selectresults
     
     if(is.element("Serving correlations",theresults)){
       
       foodnamesused = input$thefoodnames
       nfused = length(foodnamesused)
       foodindex = match(foodnamesused,foodnames)
       
       if(nfused>1){  # generate a sample of positive consumptions and plot in pairs
         nsample <- 500
         sampledmus <- matrix(NA,nsample,nf)
         sampledsw <- matrix(NA,nsample,nf)
         mcsample <- round(seq(1,(input$Iterations-burnin),length=nsample))
         for(i in 1:nsample){
         sampledmus[i,1:nf] <- rmvnorm(1,mus0[mcsample[i],1:nf],solve(Ts0[mcsample[i],1:nf,1:nf]))
         sampledsw[i,1:nf] <- exp(rmvnorm(1,sampledmus[i,1:nf],solve(Ts[mcsample[i],1:nf,1:nf])))
         }
         
         datasw <- matrix(NA,nr*nd,nf)
         rowindex <- 0
         for(r in 1:nr){
           for(t in 1:nd){
             rowindex <- rowindex+1
             datasw[rowindex,1:nf]<-exp(logsw[r,t,1:nf])
           }
         }
         group <- c(rep(1,nr*nd),rep(2,nsample)) # groups for data values and simulated values
         DF1 <- data.frame(log10(datasw[1:(nr*nd),foodindex]))
         DF2 <- data.frame(log10(sampledsw[1:nsample,foodindex]))
         colnames(DF1) <- foodnamesused
         colnames(DF2) <- foodnamesused
         pairs(rbind(DF1,DF2),
               main="Pairwise scatterplots of log (consumption/bw+)",
               cex=c(1.8,0.5)[group],pch=c(8,16)[group],col=c("blue",raspberry)[group])
       } # nfused >1 
     } # if serving correlations
   })
   
   ################################################
   
   output$distPlot7 <- renderPlot({     
     # generate results based on inputs from ui.R: 
     # Correlation plots for mean consumptions
     
     currentresults()
     
     theresults = input$selectresults
     
     if(is.element("Mean serving correlations",theresults)){
       
       foodnamesused = input$thefoodnames
       nfused = length(foodnamesused)
       foodindex = match(foodnamesused,foodnames)
       
       if(nfused>1){  # generate a sample of positive consumptions and plot in pairs
         nsample <- 500
         sampledmus <- matrix(NA,nsample,nf) # for the means in log-scale
         sampledmeans <- matrix(NA,nsample,nf) # for the means in absolute scale
         mcsample <- round(seq(1,(input$Iterations-burnin),length=nsample))
         for(i in 1:nsample){
           sampledmus[i,1:nf] <- rmvnorm(1,mus0[mcsample[i],1:nf],solve(Ts0[mcsample[i],1:nf,1:nf]))
           sampledmeans[i,1:nf] <- exp(sampledmus[i,1:nf]+0.5*sigs[mcsample[i],1:nf]^2) 
         }
         datameansw <- matrix(NA,nr,nf)
         for(r in 1:nr){
           for(i in 1:nf){
             datameansw[r,i]<- mean(exp(logsw[r,1:nd,i]),na.rm=TRUE)
           }
         }
         group <- c(rep(1,nr),rep(2,nsample)) # groups for data values and simulated values
         DF1 <- data.frame(log10(datameansw[1:nr,foodindex]))
         DF2 <- data.frame(log10(sampledmeans[1:nsample,foodindex]))
         colnames(DF1) <- foodnamesused
         colnames(DF2) <- foodnamesused
         pairs(rbind(DF1,DF2),
               main="Pairwise scatterplots of log (E(consumption/bw+))",
               cex=c(1.8,0.5)[group],pch=c(8,16)[group],col=c("blue",raspberry)[group])
       } 
     } # if mean serving correlations
   })
   
   #######################################################################
   
   resultValues <- reactive({
     # generate results based on inputs from ui.R: 
     # create data frame containing posterior predictive summaries
     
     currentresults()

     theresults = input$selectresults
       
     foodnamesused = input$thefoodnames
     nfused = length(foodnamesused)
     foodindex = match(foodnamesused,foodnames)
     
     hazardnamesused = input$thehazardnames
     hazardtypesused = hazardtypes[is.element(hazardnames,hazardnamesused)]
     nhused = length(hazardnamesused)
     hazardnamesK = hazardnames[hazardtypes=="K"]
     hazardnamesM = hazardnames[hazardtypes=="M"]
     hazardnamesusedK = hazardnamesused[hazardtypesused=="K"]
     hazardnamesusedM = hazardnamesused[hazardtypesused=="M"]
     nhusedK = length(hazardnamesusedK)
     nhusedM = length(hazardnamesusedM)
     hazardindex = match(hazardnamesused,hazardnames)
     hazardindexK = match(hazardnamesusedK,hazardnamesK)
     hazardindexM = match(hazardnamesusedM,hazardnamesM)
     Rall <- input$factor # optional processing factors
     Pall <- input$pfactor  # optional prevalence factors
     
     if(is.element("Posterior predictive",theresults)){
       
       if(is.element("Concentrations",theresults)){
         
         DF <- data.frame(Results="No food-hazard selected")
         
         if((nhusedK>0)&(nfused>0)){  # if some chemical hazard in some food selected
           # redefine dimensions if scalars:
           if((nhK==1)&(nf==1)){
             mucK <- array(mucK,dim=c(input$Iterations-burnin,1,1))
             sigcK <- array(sigcK,dim=c(input$Iterations-burnin,1,1))
             pK <- array(pK,dim=c(input$Iterations-burnin,1,1))
           }    
         cKmc <- array(NA,dim=c(input$Iterations-burnin,nhusedK,nfused))   
         for(mc in 1:(input$Iterations-burnin)){  # simulate posterior predictive concentrations
           for(h in 1:nhusedK){
             cKmc[mc,h,1:nfused] <- rlnorm(nfused,mucK[mc,hazardindexK[h],foodindex[1:nfused]],sigcK[mc,hazardindexK[h],foodindex[1:nfused]]) # actual contamination level
           } # end of h
         } # end of mc
         hazardnamesusedKinfoodnamesused <- array(NA,nhusedK*nfused) # collect the names of used hazard-food combinations
         hlo98cK <- numeric()
         hup98cK <- numeric()
         hlo90cK <- numeric()
         hup90cK <- numeric()
         hlo80cK <- numeric()
         hup80cK <- numeric()
         hmediancK <- numeric()
         counterK <- 0
         for(i in 1:nhusedK){
           for(j in 1:nfused){
             counterK <- counterK +1
             hazardnamesusedKinfoodnamesused[counterK] <- paste0(hazardnamesusedK[i]," in ",foodnamesused[j],":") 
             hlo98cK[counterK] <- quantile(cKmc[,i,j],c(0.01),names=FALSE) # calculate quantile
             hup98cK[counterK] <- quantile(cKmc[,i,j],c(0.99),names=FALSE) # calculate quantile
             hlo90cK[counterK] <- quantile(cKmc[,i,j],c(0.05),names=FALSE) # calculate quantile
             hup90cK[counterK] <- quantile(cKmc[,i,j],c(0.95),names=FALSE) # calculate quantile
             hlo80cK[counterK] <- quantile(cKmc[,i,j],c(0.10),names=FALSE) # calculate quantile
             hup80cK[counterK] <- quantile(cKmc[,i,j],c(0.90),names=FALSE) # calculate quantile
             hmediancK[counterK] <- quantile(cKmc[,i,j],c(0.5),names=FALSE) # calculate quantile
           }
         } 
           DFKconcentrations <- data.frame(
             Quantity = paste(hazardnamesusedKinfoodnamesused,"concentr+"),
             Q01 = as.character(round(hlo98cK[1:counterK],2)),
             Q05 = as.character(round(hlo90cK[1:counterK],2)),
             Q10 = as.character(round(hlo80cK[1:counterK],2)),
             Median = as.character(round(hmediancK[1:counterK],2)),
             Q90 = as.character(round(hup80cK[1:counterK],2)),
             Q95 = as.character(round(hup90cK[1:counterK],2)),
             Q99 = as.character(round(hup98cK[1:counterK],2)),
             stringsAsFactors=FALSE)
         } # end of if nhusedK nfused  
           
         if((nhusedM>0)&(nfused>0)){  # if some microbial hazard in some food selected
           # redefine dimensions if scalars:
           if((nhM==1)&(nf==1)){
             mucM <- array(mucM,dim=c(input$Iterations-burnin,1,1))
             sigcM <- array(sigcM,dim=c(input$Iterations-burnin,1,1))
             pM <- array(pM,dim=c(input$Iterations-burnin,1,1))
           } 
           cMmc <- array(NA,dim=c(input$Iterations-burnin,nhusedM,nfused))
           for(mc in 1:(input$Iterations-burnin)){  # simulate posterior predictive concentrations
             for(h in 1:nhusedM){
               cMmc[mc,h,1:nfused] <- rlnorm(nfused,mucM[mc,hazardindexM[h],foodindex[1:nfused]],sigcM[mc,hazardindexM[h],foodindex[1:nfused]]) # actual contamination level
             } # end of h
           } # end of mc
           hazardnamesusedMinfoodnamesused <- array(NA,nhusedM*nfused) # collect the names of used hazard-food combinations
           hlo98cM <- numeric()
           hup98cM <- numeric()
           hlo90cM <- numeric()
           hup90cM <- numeric()
           hlo80cM <- numeric()
           hup80cM <- numeric()
           hmediancM <- numeric()
           counterM <- 0
           for(i in 1:nhusedM){
             for(j in 1:nfused){
               counterM <- counterM +1
               hazardnamesusedMinfoodnamesused[counterM] <- paste0(hazardnamesusedM[i]," in ",foodnamesused[j],":") 
               hlo98cM[counterM] <- quantile(cMmc[,i,j],c(0.01),names=FALSE) # calculate quantile
               hup98cM[counterM] <- quantile(cMmc[,i,j],c(0.99),names=FALSE) # calculate quantile
               hlo90cM[counterM] <- quantile(cMmc[,i,j],c(0.05),names=FALSE) # calculate quantile
               hup90cM[counterM] <- quantile(cMmc[,i,j],c(0.95),names=FALSE) # calculate quantile
               hlo80cM[counterM] <- quantile(cMmc[,i,j],c(0.10),names=FALSE) # calculate quantile
               hup80cM[counterM] <- quantile(cMmc[,i,j],c(0.90),names=FALSE) # calculate quantile
               hmediancM[counterM] <- quantile(cMmc[,i,j],c(0.5),names=FALSE) # calculate quantile
             }
           } 
           DFMconcentrations <- data.frame(
             Quantity = paste(hazardnamesusedMinfoodnamesused,"concentr+"),
             Q01 = as.character(round(hlo98cM[1:counterM],2)),
             Q05 = as.character(round(hlo90cM[1:counterM],2)),
             Q10 = as.character(round(hlo80cM[1:counterM],2)),
             Median = as.character(round(hmediancM[1:counterM],2)),
             Q90 = as.character(round(hup80cM[1:counterM],2)),
             Q95 = as.character(round(hup90cM[1:counterM],2)),
             Q99 = as.character(round(hup98cM[1:counterM],2)),
             stringsAsFactors=FALSE)
         } # end of if nhusedK nfused  
         
       } # end of concentrations
       
     if(is.element("Consumptions",theresults)|is.element("Exposures",theresults)){
     
       # redefine dimensions if scalars:
       if(nf==1){
         mus0 <- matrix(mus0,input$Iterations-burnin,1)  
         sigs <- matrix(sigs,input$Iterations-burnin,1) 
         ppred <- matrix(ppred,input$Iterations-burnin,1)
         if(input$modelchoice=="Independent days"){
         if(input$modelchoice2=="Yes"){ # model with between user variability in use frequencies
           Tp <- array(Tp,dim=c(input$Iterations-burnin,1,1)) 
         }
         logitp0 <- matrix(logitp0,input$Iterations-burnin,1)
         }
         Ts <- array(Ts,dim=c(input$Iterations-burnin,1,1))
         Ts0 <- array(Ts0,dim=c(input$Iterations-burnin,1,1))
       }   
     DF <- data.frame(Results="No food-hazard selected")
     
     if(nfused>0){
       
     if(nhusedK>0){ # formatting for posterior predictive distributions for chemical hazards
       logitpmc <- matrix(NA,input$Iterations-burnin,nf)
       pmc <- matrix(NA,input$Iterations-burnin,nf)
       musmc <- matrix(NA,input$Iterations-burnin,nf)
       EemcK <- array(NA,dim=c(input$Iterations-burnin,nhusedK,nfused)) 
       chronictotbwK <- matrix(NA,(input$Iterations-burnin),nhusedK)
       RK = matrix(NA,nf,nhK)
       RK[1:nf,1:nhK] = Rall[1:nf,is.element(hazardnames,hazardnamesusedK)]
       logRK = log(RK)
       PK = matrix(NA,nf,nhK)
       PK[1:nf,1:nhK] = Pall[1:nf,is.element(hazardnames,hazardnamesusedK)]
       hlo90totbwK <- numeric(nhusedK)
       hup90totbwK <- numeric(nhusedK)
       hmediantotbwK <- numeric(nhusedK)
       hlo80totbwK <- numeric(nhusedK)
       hup80totbwK <- numeric(nhusedK)
       hlo98totbwK <- numeric(nhusedK)
       hup98totbwK <- numeric(nhusedK)
     } # end of nhusedK
     if(nhusedM>0){ # formatting for posterior predictive distributions for microbial hazards
       logitpmc <- matrix(NA,input$Iterations-burnin,nf)
       pmc <- matrix(NA,input$Iterations-burnin,nf)
       musmc <- matrix(NA,input$Iterations-burnin,nf)
       wmc <- numeric()
       Umc <- matrix(NA,input$Iterations-burnin,nf)
       smc <- matrix(NA,input$Iterations-burnin,nf)
       Imc <- array(NA,dim=c(input$Iterations-burnin,nhM,nf))
       cmc <- array(NA,dim=c(input$Iterations-burnin,nhM,nf))
       EemcM <- array(NA,dim=c(input$Iterations-burnin,nhusedM,nfused))    
       acutetotM <- matrix(NA,(input$Iterations-burnin),nhusedM)
       RM = matrix(NA,nf,nhM)
       RM[1:nf,1:nhM] = Rall[1:nf,is.element(hazardnames,hazardnamesusedM)]
       logRM = log(RM)
       PM = matrix(NA,nf,nhM)
       PM[1:nf,1:nhM] = Pall[1:nf,is.element(hazardnames,hazardnamesusedM)]
       hlo90totacuteM <- numeric(nhusedM)
       hup90totacuteM <- numeric(nhusedM)
       hmediantotacuteM <- numeric(nhusedM)
       hlo80totacuteM <- numeric(nhusedM)
       hup80totacuteM <- numeric(nhusedM)
       hlo98totacuteM <- numeric(nhusedM)
       hup98totacuteM <- numeric(nhusedM)
     } # end of nhusedM
         
     if(nhusedK>0){ # simulate posterior predictive distributions for chemical hazards
     # redefine dimensions if scalars:
     if((nhK==1)&(nf==1)){
         mucK <- array(mucK,dim=c(input$Iterations-burnin,1,1))
         sigcK <- array(sigcK,dim=c(input$Iterations-burnin,1,1))
         pK <- array(pK,dim=c(input$Iterations-burnin,1,1))
     }       
     if(input$modelchoice == "Independent days"){   
     for(mc in 1:(input$Iterations-burnin)){
       if(nf>1){ # many foods
         if(input$modelchoice2=="Yes"){ # variability between users frequencies   
         logitpmc[mc,1:nf] <- rmvnorm(1,logitp0[mc,1:nf],solve(Tp[mc,1:nf,1:nf]))   
         }
         if(input$modelchoice2=="No"){ # no variability between users frequencies
         logitpmc[mc,1:nf] <- logitp0[mc,1:nf]       
         }   
       pmc[mc,1:nf] <- exp(logitpmc[mc,1:nf])/(1+exp(logitpmc[mc,1:nf])) # individual use probability
       musmc[mc,1:nf] <- rmvnorm(1,mus0[mc,1:nf],solve(Ts0[mc,1:nf,1:nf])) # individual mean log amount
       }
       if(nf==1){ # only one food
         if(input$modelchoice2=="Yes"){ # variability between users frequencies   
         logitpmc[mc,1] <- rnorm(1,logitp0[mc,1],sqrt(1/Tp[mc,1,1]))
         }
         if(input$modelchoice2=="No"){ # no variability between users frequencies
         logitpmc[mc,1] <- logitp0[mc,1]    
         }   
       pmc[mc,1] <- exp(logitpmc[mc,1])/(1+exp(logitpmc[mc,1])) # individual use probability
       musmc[mc,1] <- rnorm(1,mus0[mc,1],sqrt(1/Ts0[mc,1,1])) # individual mean log amount     
       }
       
       for(h in 1:nhusedK){
         for(i in 1:nfused){
           # chronic (mean) exposure for a random consumer, hazard h, food i:
         EemcK[mc,h,i] <- pK[mc,hazardindexK[h],foodindex[i]]*PK[foodindex[i],hazardindexK[h]]*
                          pmc[mc,foodindex[i]]*exp(logRK[foodindex[i],hazardindexK[h]]
                                                  +musmc[mc,foodindex[i]]
                                                  +0.5*sigs[mc,foodindex[i]]^2
                                                  +mucK[mc,hazardindexK[h],foodindex[i]]
                                                  +0.5*sigcK[mc,hazardindexK[h],foodindex[i]]^2) 
         }
         chronictotbwK[mc,h] <- sum(EemcK[mc,h,1:nfused]) # sum over foods
       } # end of h
     } # end of mc
     } # end of independent days
       
     if(input$modelchoice == "Dependent days"){   
     for(mc in 1:(input$Iterations-burnin)){
       if(nf>1){ # many foods
       musmc[mc,1:nf] <- rmvnorm(1,mus0[mc,1:nf],solve(Ts0[mc,1:nf,1:nf])) # individual mean log amount
       }
       if(nf==1){ # only one food
       musmc[mc,1] <- rnorm(1,mus0[mc,1],sqrt(1/Ts0[mc,1,1])) # individual mean log amount   
       }
       for(h in 1:nhusedK){
            for(i in 1:nfused){
               # chronic (mean) exposure for a random consumer, hazard h, food i:
               EemcK[mc,h,i] <- pK[mc,hazardindexK[h],foodindex[i]]*PK[foodindex[i],hazardindexK[h]]*
                                ppred[mc,foodindex[i]]*exp(logRK[foodindex[i],hazardindexK[h]]
                                                    +musmc[mc,foodindex[i]]
                                                    +0.5*sigs[mc,foodindex[i]]^2
                                                    +mucK[mc,hazardindexK[h],foodindex[i]]
                                                    +0.5*sigcK[mc,hazardindexK[h],foodindex[i]]^2) 
            }
             chronictotbwK[mc,h] <- sum(EemcK[mc,h,1:nfused]) # sum over foods
       } # end of h
     } # end of mc
     } # end of dependent days
            
     } # end of nhusedK 
       
     if(nhusedM>0){  # simulate posterior predictive distributions for microbial hazards
       # redefine dimensions if scalars:
       if((nhM==1)&(nf==1)){
         mucM <- array(mucM,dim=c(input$Iterations-burnin,1,1))
         sigcM <- array(sigcM,dim=c(input$Iterations-burnin,1,1))
         pM <- array(pM,dim=c(input$Iterations-burnin,1,1))
       }  
       
     if(input$modelchoice == "Independent days"){    
     for(mc in 1:(input$Iterations-burnin)){
       if(nf>1){  # many foods
       if(input$modelchoice2=="Yes"){ # variability between users frequencies    
       logitpmc[mc,1:nf] <- rmvnorm(1,logitp0[mc,1:nf],solve(Tp[mc,1:nf,1:nf]))
       }
       if(input$modelchoice2=="No"){ # no variability between users frequencies
       logitpmc[mc,1:nf] <- logitp0[mc,1:nf]   
       }   
       pmc[mc,1:nf] <- exp(logitpmc[mc,1:nf])/(1+exp(logitpmc[mc,1:nf])) # individual use probability
       musmc[mc,1:nf] <- rmvnorm(1,mus0[mc,1:nf],solve(Ts0[mc,1:nf,1:nf])) # individual mean amount
       } 
       if(nf==1){ # only one food
       if(input$modelchoice2=="Yes"){ # variability between users frequencies
       logitpmc[mc,1] <- rnorm(1,logitp0[mc,1],sqrt(1/Tp[mc,1,1])) 
       }
       if(input$modelchoice2=="No"){ # no variability between users frequencies
       logitpmc[mc,1] <- logitp0[mc,1]     
       }   
       pmc[mc,1] <- exp(logitpmc[mc,1])/(1+exp(logitpmc[mc,1])) # individual use probability
       musmc[mc,1] <- rnorm(1,mus0[mc,1],sqrt(1/Ts0[mc,1,1])) # individual mean amount  
       }
       wmc[mc] <- rlnorm(1,muw,sigw) # bodyweight for random individual
       Umc[mc,1:nf] <- rbinom(nf,rep(1,nf),pmc[mc,1:nf]) # actual random use
       smc[mc,1:nf] <- rlnorm(nf,musmc[mc,1:nf],sigs[mc,1:nf]) # actual random amount
       for(h in 1:nhM){
         Imc[mc,h,1:nf] <- rbinom(nf,rep(1,nf),pM[mc,h,1:nf]*PM[1:nf,h]) # actual contamination yes/no
         cmc[mc,h,1:nf] <- rlnorm(nf,mucM[mc,h,1:nf],sigcM[mc,h,1:nf]) # actual contamination level
       } 
       
       for(h in 1:nhusedM){ #Predict final count with poisson distribution:
         for(i in 1:nfused){
           # acute exposure for a random consumer, hazard h, food i:
           EemcM[mc,h,i] <- Imc[mc,hazardindexM[h],foodindex[i]]*Umc[mc,foodindex[i]]*
            smc[mc,foodindex[i]]*RM[foodindex[i],hazardindexM[h]]*
            cmc[mc,hazardindexM[h],foodindex[i]]*wmc[mc]
         }  
         # sum(exposure.acuteM[mc,hazardindexM[h],foodindex[1:nfused]])<=10000
         if(sum(EemcM[mc,h,1:nfused])<=10000){ # use Poisson when the mean is 'small'
           acutetotM[mc,h] <- rpois(1,sum(EemcM[mc,h,1:nfused])) 
         }
         # sum(exposure.acuteM[mc,hazardindexM[h],foodindex[1:nfused]])>10000
         if(sum(EemcM[mc,h,1:nfused])>10000){ # use rounded Normal when the mean is 'large'
           acutetotM[mc,h] <- round(rnorm(1,sum(EemcM[mc,h,1:nfused]),sqrt(sum(EemcM[mc,h,1:nfused])))) 
         }
       } # end of h
       } # end of mc
       } # end of independent days
       if(input$modelchoice == "Dependent days"){    
         for(mc in 1:(input$Iterations-burnin)){
           if(nf>1){ # many foods
           musmc[mc,1:nf] <- rmvnorm(1,mus0[mc,1:nf],solve(Ts0[mc,1:nf,1:nf])) # individual mean log amount
           }
           if(nf==1){ # only one food
           musmc[mc,1] <- rnorm(1,mus0[mc,1],sqrt(1/Ts0[mc,1,1])) # individual mean log amount
           }
           wmc[mc] <- rlnorm(1,muw,sigw) # bodyweight for random individual
           Umc[mc,1:nf] <- rbinom(nf,rep(1,nf),ppred[mc,1:nf]) # actual random use
           smc[mc,1:nf] <- rlnorm(nf,musmc[mc,1:nf],sigs[mc,1:nf]) # actual random amount
           for(h in 1:nhM){
             Imc[mc,h,1:nf] <- rbinom(nf,rep(1,nf),pM[mc,h,1:nf]*PM[1:nf,h]) # actual contamination yes/no
             cmc[mc,h,1:nf] <- rlnorm(nf,mucM[mc,h,1:nf],sigcM[mc,h,1:nf]) # actual contamination level
           } 
           
           for(h in 1:nhusedM){ #Predict final count with poisson distribution:
             for(i in 1:nfused){
               # acute exposure for a random consumer, hazard h, food i:
               EemcM[mc,h,i] <- Imc[mc,hazardindexM[h],foodindex[i]]*Umc[mc,foodindex[i]]*
                smc[mc,foodindex[i]]*RM[foodindex[i],hazardindexM[h]]*
                cmc[mc,hazardindexM[h],foodindex[i]]*wmc[mc]
             }  
             # sum(exposure.acuteM[mc,hazardindexM[h],foodindex[1:nfused]])<=10000
             if(sum(EemcM[mc,h,1:nfused])<=10000){ # use Poisson when the mean is 'small'
               acutetotM[mc,h] <- rpois(1,sum(EemcM[mc,h,1:nfused])) 
             }
             # sum(exposure.acuteM[mc,hazardindexM[h],foodindex[1:nfused]])>10000
             if(sum(EemcM[mc,h,1:nfused])>10000){ # use rounded Normal when the mean is 'large'
               acutetotM[mc,h] <- round(rnorm(1,sum(EemcM[mc,h,1:nfused]),sqrt(sum(EemcM[mc,h,1:nfused]))))
             }
           } # end of h
         } # end of mc
       } # end of dependent days
       
     } # end of nhusedM
       
     ############ Get posterior predictive results into data frames:  #################
       
     if(nhusedK>0){
     for(h in 1:nhusedK){ # posterior predictive summaries (quantiles) of individual exposures (chronic)
       hlo90totbwK[h] <- quantile(chronictotbwK[,h],0.05,names=FALSE)
       hmediantotbwK[h] <- quantile(chronictotbwK[,h],0.5,names=FALSE)
       hup90totbwK[h] <- quantile(chronictotbwK[,h],0.95,names=FALSE)
       hlo80totbwK[h] <- quantile(chronictotbwK[,h],0.10,names=FALSE)
       hup80totbwK[h] <- quantile(chronictotbwK[,h],0.90,names=FALSE)
       hlo98totbwK[h] <- quantile(chronictotbwK[,h],0.01,names=FALSE)
       hup98totbwK[h] <- quantile(chronictotbwK[,h],0.99,names=FALSE)
     }
     }
     if(nhusedM>0){
     for(h in 1:nhusedM){  # posterior predictive summaries (quantiles) of individual exposures (acute)
       hlo90totacuteM[h] <- round(quantile(acutetotM[,h],0.05,names=FALSE))
       hmediantotacuteM[h] <- round(quantile(acutetotM[,h],0.5,names=FALSE))
       hup90totacuteM[h] <- round(quantile(acutetotM[,h],0.95,names=FALSE))
       hlo80totacuteM[h] <- round(quantile(acutetotM[,h],0.10,names=FALSE))
       hup80totacuteM[h] <- round(quantile(acutetotM[,h],0.90,names=FALSE))
       hlo98totacuteM[h] <- round(quantile(acutetotM[,h],0.01,names=FALSE))
       hup98totacuteM[h] <- round(quantile(acutetotM[,h],0.99,names=FALSE))
     }
     }
     if(is.element("Consumptions",theresults)){
       musmc <- matrix(NA,input$Iterations-burnin,nf)
       fconslo90bw <- numeric(nfused)
       fconsup90bw <- numeric(nfused)
       fconsmedianbw <- numeric(nfused)
       fconslo80bw <- numeric(nfused)
       fconsup80bw <- numeric(nfused)
       fconslo98bw <- numeric(nfused)
       fconsup98bw <- numeric(nfused)
       fconslo90 <- numeric(nfused)
       fconsup90 <- numeric(nfused)
       fconsmedian <- numeric(nfused)
       fconslo80 <- numeric(nfused)
       fconsup80 <- numeric(nfused)
       fconslo98 <- numeric(nfused)
       fconsup98 <- numeric(nfused)
     for(mc in 1:(input$Iterations-burnin)){
       if(nf>1){ # many foods
       musmc[mc,1:nf] <- rmvnorm(1,mus0[mc,1:nf],solve(Ts0[mc,1:nf,1:nf]))
       }
       if(nf==1){ # only one food
       musmc[mc,1] <- rnorm(1,mus0[mc,1],sqrt(1/Ts0[mc,1,1]))   
       }
     }
     for(i in 1:nfused){ # posterior predictive summaries (quantiles) of individual chronic consumptions (/bw and absolute)
       fconslo90bw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.05,names=FALSE)
       fconslo90[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.05,names=FALSE)
       fconsup90bw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.95,names=FALSE)
       fconsup90[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.95,names=FALSE) 
       fconslo80bw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.10,names=FALSE)
       fconslo80[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.10,names=FALSE)
       fconsup80bw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.90,names=FALSE)
       fconsup80[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.90,names=FALSE)
       fconslo98bw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.01,names=FALSE) 
       fconslo98[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.01,names=FALSE)
       fconsup98bw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.99,names=FALSE) 
       fconsup98[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.99,names=FALSE)
       fconsmedianbw[i] <- quantile(exp(musmc[,foodindex[i]]+0.5*sigs[,foodindex[i]]^2),0.5,names=FALSE)
       fconsmedian[i] <- quantile(exp(musmc[,foodindex[i]]+muw+0.5*sigs[,foodindex[i]]^2+0.5*sigw^2),0.5,names=FALSE)
     } # end of nfused
     } # end of if consumptions
     
     # Compose data frame for chemical exposure
     if(nhusedK>0){
     DF1K <- data.frame(
       Quantity = paste(hazardnamesusedK,"total chronic exposure/bw"),
       Q01 = as.character(round(hlo98totbwK[1:nhusedK],2)),
       Q05 = as.character(round(hlo90totbwK[1:nhusedK],2)),
       Q10 = as.character(round(hlo80totbwK[1:nhusedK],2)),
       Median = as.character(round(hmediantotbwK[1:nhusedK],2)),
       Q90 = as.character(round(hup80totbwK[1:nhusedK],2)),
       Q95 = as.character(round(hup90totbwK[1:nhusedK],2)),
       Q99 = as.character(round(hup98totbwK[1:nhusedK],2)),
       stringsAsFactors=FALSE)
     }
     # Compose data frame for microbial exposure   
     if(nhusedM>0){
     DF1M <- data.frame(
       Quantity = paste(hazardnamesusedM,"total acute exposure"),
       Q01 = as.character(round(hlo98totacuteM[1:nhusedM],2)),
       Q05 = as.character(round(hlo90totacuteM[1:nhusedM],2)),
       Q10 = as.character(round(hlo80totacuteM[1:nhusedM],2)),
       Median = as.character(round(hmediantotacuteM[1:nhusedM],2)),
       Q90 = as.character(round(hup80totacuteM[1:nhusedM],2)),
       Q95 = as.character(round(hup90totacuteM[1:nhusedM],2)),
       Q99 = as.character(round(hup98totacuteM[1:nhusedM],2)),
       stringsAsFactors=FALSE)
     }
     
     # Compose data frame for consumptions   
     if(is.element("Consumptions",theresults)){   
     DF4 <- data.frame(
       Quantity = paste(foodnamesused,"mean daily use/bw+"),
       Q01 = as.character(round(fconslo98bw[1:nfused],2)),
       Q05 = as.character(round(fconslo90bw[1:nfused],2)),
       Q10 = as.character(round(fconslo80bw[1:nfused],2)),
       Median = as.character(round(fconsmedianbw[1:nfused],2)),
       Q90 = as.character(round(fconsup80bw[1:nfused],2)),
       Q95 = as.character(round(fconsup90bw[1:nfused],2)),
       Q99 = as.character(round(fconsup98bw[1:nfused],2)),
       stringsAsFactors=FALSE)
     DF5 <- data.frame(
       Quantity = paste(foodnamesused,"mean daily use+"),
       Q01 = as.character(round(fconslo98[1:nfused],2)),
       Q05 = as.character(round(fconslo90[1:nfused],2)),
       Q10 = as.character(round(fconslo80[1:nfused],2)),
       Median = as.character(round(fconsmedian[1:nfused],2)),
       Q90 = as.character(round(fconsup80[1:nfused],2)),
       Q95 = as.character(round(fconsup90[1:nfused],2)),
       Q99 = as.character(round(fconsup98[1:nfused],2)),
       stringsAsFactors=FALSE)
     }
     
     } # end of nfused>0
     } # end of consumption or exposure
          
     if(!is.element("Concentrations",theresults)){   
     if(is.element("Consumptions",theresults)&!is.element("Exposures",theresults)){
      DF <- rbind.data.frame(DF4,DF5)
     }
     if(is.element("Exposures",theresults)&!is.element("Consumptions",theresults)){
       if((nhusedK>0)&(nhusedM>0)){ DF <- rbind.data.frame(DF1K,DF1M) }
       if((nhusedK>0)&(nhusedM==0)){ DF <- rbind.data.frame(DF1K) }
       if((nhusedK==0)&(nhusedM>0)){ DF <- rbind.data.frame(DF1M) }
     }
     if(is.element("Exposures",theresults)&is.element("Consumptions",theresults)){
       if((nhusedK>0)&(nhusedM>0)){ DF <- rbind.data.frame(DF1K,DF1M,DF4,DF5) }
       if((nhusedK>0)&(nhusedM==0)){ DF <- rbind.data.frame(DF1K,DF4,DF5) }
       if((nhusedK==0)&(nhusedM>0)){ DF <- rbind.data.frame(DF1M,DF4,DF5) }
     }
     } # end of if !Concentrations
       if(is.element("Concentrations",theresults)){
         if((nhusedK>0)&(nhusedM==0)){DF <- rbind.data.frame(DFKconcentrations)}
         if((nhusedK==0)&(nhusedM>0)){DF <- rbind.data.frame(DFMconcentrations)}
         if((nhusedK>0)&(nhusedM>0)){DF <- rbind.data.frame(DFKconcentrations,DFMconcentrations)}
         
         if(is.element("Consumptions",theresults)&!is.element("Exposures",theresults)){
           DF <- rbind.data.frame(DF,DF4,DF5)
         }
         if(is.element("Exposures",theresults)&!is.element("Consumptions",theresults)){
           if((nhusedK>0)&(nhusedM>0)){ DF <- rbind.data.frame(DF,DF1K,DF1M) }
           if((nhusedK>0)&(nhusedM==0)){ DF <- rbind.data.frame(DF,DF1K) }
           if((nhusedK==0)&(nhusedM>0)){ DF <- rbind.data.frame(DF,DF1M) }
         }
         if(is.element("Exposures",theresults)&is.element("Consumptions",theresults)){
           if((nhusedK>0)&(nhusedM>0)){ DF <- rbind.data.frame(DF,DF1K,DF1M,DF4,DF5) }
           if((nhusedK>0)&(nhusedM==0)){ DF <- rbind.data.frame(DF,DF1K,DF4,DF5) }
           if((nhusedK==0)&(nhusedM>0)){ DF <- rbind.data.frame(DF,DF1M,DF4,DF5) }
         }
       } # end of if Concentrations
       
       
     DF # results collected as data frame
     } # end of theresults selection (posterior predictive)    
   }) 
 
   
   headertext <- reactive({ 
     currentresults()
     theresults = input$selectresults
     if(is.element("Posterior predictive",theresults)&(is.element("Concentrations",theresults)|is.element("Consumptions",theresults)|is.element("Exposures",theresults)) ){"POSTERIOR PREDICTIVE DISTRIBUTION SUMMARIES"}
   })
  
   # Show the values using an HTML table
   
    output$values <- renderTable({
      resultValues()
   })
   
   output$header <- renderText({
      headertext()
    })
       
}
 
# Run the application 
shinyApp(ui = ui, server = server)


