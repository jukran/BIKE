# BIKE
Foodborne exposure assessment model for biological and chemical hazards

This is a Bayesian model for estimating foodborne exposures based on dietary data on food consumption frequencies and amounts, and corresponding occurrence data on hazard concentrations (chemical and microbiological) and hazard prevalence (microbiological). The model can be used as a shiny app providing a user interface. Data need to be specified in files that are read by R-code, which then creates the necessary BUGS code of the Bayesian model which is simulated using OpenBUGS in the background. Results are processed in R and presented in the shiny app. Users need to have both R and OpenBUGS installed. 

Read the full paper: https://www.mdpi.com/2304-8158/10/11/2520 
