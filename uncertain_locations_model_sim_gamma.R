uncertain_locations_model_sim_gamma <- nimbleCode({
  
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  gamma_par ~ dunif(0.1,10)

  for (i in 1:N_events) {
    
    aux_cat[1:4] <- pow(weights[i,1:4],gamma_par)
    segment_rand[i] ~ dcat(aux_cat[1:4])
    s[i] <- segments[i,segment_rand[i]]
    lambda[i] <- lambda0*exp(alpha + beta1*X1[s[i]] + beta2*X2[s[i]])
    log(L[i]) <- I[i]*log(lambda[i])-w[i]*lambda[i] # I_i = 1 if event, 0 otherwise; w[i] = 'segment length' if segment point, 0 otherwise
    phi[i] <- -log(L[i])+C
    zeros[i] ~ dpois(phi[i])
    
  }
  
  for (i in (N_events+1):N) {
    
    lambda[i] <- lambda0*exp(alpha + beta1*X1[i-N_events] + beta2*X2[i-N_events])
    log(L[i]) <- I[i]*log(lambda[i])-w[i]*lambda[i] # I_i = 1 if event, 0 otherwise; w[i] = 'segment length' if segment point, 0 otherwise
    phi[i] <- -log(L[i])+C
    zeros[i] ~ dpois(phi[i])
    
  }
  
})