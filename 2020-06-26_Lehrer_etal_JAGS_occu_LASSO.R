model {
  
  # Priors on betas
  # LASSO on parameters
  # double exponential distribution
  for(spec in 1:nspec){
    b0[spec] ~ dnorm(0,2) # priors on occupancy intercept
    a0[spec] ~ dnorm(0,2) # prior on detection intercept
    a_pcp[spec] ~ ddexp(0, lambda[spec]) # prior on precip detection coefficient
    lambda[spec] ~ dgamma(1,1)
    for(beta in 1:ncof){
      b[spec,beta] ~ ddexp(0, lambda[spec]) # prior on occupancy betas
    }
  }
  
  # Priors on detection parameter for detectors
  for(detect in 1:ndetector){
    a_det[detect] ~ ddexp(0, lambda_det)
  }
  lambda_det ~ dgamma(1,1)
  
  # Prior for session
  for(spec in 1:nspec){
    sigma[spec] ~ dgamma(1,1)
    tau[spec] <- 1/sqrt(sigma[spec])
    for(session in 1:nsession){
      u[spec,session] ~ dnorm(0,tau[spec])
    }
  }
      
  
  # Ecological model
  for(spec in 1:nspec){
    for(site in 1:nsite){
      for(session in 1:nsession){
        z[spec,site,session] ~ dbern(psi[spec,site,session])
        logit(psi[spec,site,session]) <- b0[spec] + inprod(b[spec,], sitecovs[site,]) + u[spec,session]
      }
    }
  }
  
  # # Observational Model
  for(spec in 1:nspec){
    for(site in 1:nsite){
      for(session in 1:nsession){
        logit(p[spec,site,session]) <- a0[spec] + a_det[detector[site,session]] +
          a_pcp[spec] * precip[site,session]
        y[spec,site,session] ~ dbin(z[spec,site,session]*p[spec,site,session], jmat[site,session])
      }
    }
  }
  
  # species richness
  for(site in 1:nsite){
    for(session in 1:nsession){
      rich[site,session] <- sum(z[,site,session])
    }
  }
  
}