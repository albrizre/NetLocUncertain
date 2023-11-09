nearest_segment_model_sim <- nimbleCode({
  
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)

  for (i in 1:N_events) {
    
    s[i] <- segments[i,1]
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