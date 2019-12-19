
# slip = data frame of slip events for EACH SEQUENCE 
likelihood2 <- function(param){
  slips <- param[1]
  p.enter <- param[2]
  p.stay  <- param[3]
  
  llh <- c()
  for (i in 1:length(slips)){
    s <- slips[[i]]
    x <- sum(s == 0)
    y <- sum(s != 0)
    z <- sum(s[which(s!=0)] - 1)
    
    # Log likelihood of each tip/anc pair
    #(1 - p.slip)^x * p.slip^y * (1-p.stay)^y * p.stay^z
    llh[i] <- 2 * y * x * z * log(1-p.slip) * log(p.slip) * log(1-p.stay) * log(p.stay)
  }
  sum(llh)
  
}


prior <- function(param){
  slips <- param[1]
  p.enter <- param[2]
  p.stay  <- param[3]
  
  sum.slips <- sum(slips)
  
  prior.pe <- dunif(p.enter, log = T)
  prior.ps <- dunif(p.stay, log = T)
  
  if (sum.slips > 0){
    prior.s <- slips
  }else{
    prior.s <- slips[sample(length(slips), sum.slips, replace=T)]
  }
    
  return(prior)
}


posterior <- function(slip){
  prior(slip) + likelihood(slip)
}

proposalFunction <- function(slip){
  return(rnorm(1,mean=slip, sd=0.0009))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,1))
  chain[1,] <- startvalue
  
  
  for (i in 1:iterations){
    
    proposal <- proposalFunction(chain[i,])
    prop <- exp(posterior(proposal) - posterior(chain[i,]))
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
  *  if (runif(1) < prop) {
      chain[i+1,] <- proposal
      
      # if the proportion is less than the random uniform sample, REJCECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
    }
    if (i %% 100 == 0){
      print(paste0("STATE ",i,": ", chain[i,1]))
    }
  }
  return(chain)
  
}

# RUN MCMC
startvalue <- 0.5
chain2 <- runMCMC(startvalue, 3000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- 10000
acceptance <- 1 - mean(duplicated(chain2[-(1:burnin),]))

# PLOTTING 
hist(chain[-(1:burnin),1],nclass=30, main="Posterior of x", xlab="True value = red line" )
abline(v = median(chain[-(1:burnin),1]))

plot(chain[-(1:burnin),1], type = "l", xlab="MCMC Steps" , main = "Chain values of x")
abline(h = median(chain[-(1:burnin),1]), col="red")