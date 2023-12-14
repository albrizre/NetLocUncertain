library(nimble)
library(spatstat)
library(spdep)
library(sp)
library(sf)

# Set the working directory...
# setwd("")

source("uncertain_locations_model_sim_gamma.R")

load("Data/grafo_intersecciones.rda")
load("Data/X1.rda")
load("Data/X2.rda")
K=4

for (error in c(0,5,10,20)){
  for (sim in 1:100){
    if (!file.exists(paste0("SimModels/uncertain_locations_model_sim_",error,"_sim",sim,"_gamma_par.rda"))){
  
      load(paste0("SimData/distances_",K,"_",error,"_sim",sim,".rda"))
      load(paste0("SimData/segments_",K,"_",error,"_sim",sim,".rda"))
      load(paste0("SimData/weights_",K,"_",error,"_sim",sim,".rda"))
      load(paste0("SimData/weights_squared_",K,"_",error,"_sim",sim,".rda"))
      
      # Model construction and call
      
      lines=st_as_sf(as.psp(grafo_intersecciones))
      segment_lengths=st_length(lines)[-1] # remove window length
      constants <- list(
        
        N = length(X1)+nrow(distances),
        N_events = nrow(distances),
        N_network = length(X1),
        w=c(rep(0,nrow(distances)),segment_lengths), 
        I=c(rep(1,nrow(distances)),rep(0,length(X1))),
        dirch_alpha=c(1,1),
        lambda0 = 1, #nrow(distances)/sum(segment_lengths),
        distances = distances,
        C=10000,
        weights =  weights,
        weights_squared =  weights_squared,
        alpha = -4
    
      )
      
      # Fit model
      
      set.seed(12345)
      data <- list(zeros = rep(0,constants$N), 
                   segments = segments,
                   X1 = X1,
                   X2 = X2)
      inits <- function() list(beta1 = 0,
                               beta2 = 0,
                               gamma_par = 2,
                               segment_rand = rep(1,constants$N_events),
                               s = as.numeric(segments[,1]))
      
      mcmc.output <- nimbleMCMC(uncertain_locations_model_sim_gamma, data = data, inits = inits, constants = constants,
                                monitors = c("beta1",
                                             "beta2",
                                             "gamma_par",
                                             "s",
                                             "lambda"), thin = 10,
                                niter = 20000, nburnin = 10000, nchains = 1, 
                                summary = TRUE, WAIC = TRUE)
      save(mcmc.output,file=paste0("SimModels/uncertain_locations_model_sim_",error,"_sim",sim,"_gamma_par.rda"))
    }
  }
}

