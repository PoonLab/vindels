
# #slips <- matrix(rep(rep(0,121),10), nrow=10,ncol=121)
# slips <- c(rep(0,121000))
# slips[sample(121000,100)] <- sample(3,100,replace = T) * 3
slips <- nchar(insertions[insertions$Count!=0, "Seq"])
slips <- slips[-167]
slips <- c(slips, rep(0, sum(nchar(insertions$Vseq))))

# randomly shuffle all the entries
rnd <- sample(length(slips),length(slips))
slips <- slips[rnd]





# EXPERIMENTAL -------------------------------
# slip = data frame of slip events for EACH SEQUENCE 

# GOAL: 
# make a likelihood function that will compute the probability of a tip/anc pair 
# it will do the following: 
# a) count the number of substitution differences 

llh.pair <- function(tip, anc, slips){
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  
}

# attempted to find the log - likelihood of an affine model 
likelihood2 <- function(slip, param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  x <- sum(slip == 0)
  y <- sum(slip != 0)
  z <- sum(slip[which(slip!=0)] - 1)
  
  llh <- c()
  # Log likelihood of each tip/anc pair
  #(1 - p.slip)^x * p.slip^y * (1-p.stay)^y * p.stay^z
  llh <- 2 * y * x * z * log(1-p.enter) * log(p.enter) * log(1-p.stay) * log(p.stay)
  sum(llh)
  
}
# --------------------------------------------------

likelihood <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  (1-p.enter)^x*(p.enter)^y*(1-p.stay)^y*(p.stay)^z
}



prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  prior.pe <- dlnorm(p.enter,meanlog=-40,sdlog=5)
  prior.ps <- dlnorm(p.stay,meanlog=-0.3,sdlog=0.15)
  
  return(prior.pe + prior.ps)
}


posterior <- function(param){
  print(prior(param))
  print(likelihood(param))
  prior(param) + likelihood(param)
}

proposalFunction <- function(param){
  p.enter <- rlnorm(1,meanlog=param[1],sdlog=0.1)
  p.stay <- rlnorm(1,meanlog=param[2],sdlog=0.01)
  return(c(p.enter,p.stay))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,2))
  chain[1,] <- startvalue
  
  
  for (i in 1:iterations){
    
    proposal <- proposalFunction(chain[i,])
    print(proposal)
    print(posterior(proposal))
    print(posterior(chain[i,]))
    prop <- exp(posterior(proposal) - posterior(chain[i,]))
    print(prop)
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    *  if (runif(1) < prop) {
      chain[i+1,] <- proposal
      
      # if the proportion is less than the random uniform sample, REJCECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
    }
    if (i %% 100 == 0){
      print(paste0("STATE ",i,": ", chain[i,]))
    }
  }
  return(chain)
  
}

# RUN MCMC
startvalue <- c(0.0000000000001,0.6)
chain2 <- runMCMC(startvalue, 5000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- 10000
acceptance <- 1 - mean(duplicated(chain2[-(1:burnin),]))

# PLOTTING 
hist(chain[-(1:burnin),1],nclass=30, main="Posterior of x", xlab="True value = red line" )
abline(v = median(chain[-(1:burnin),1]))

plot(chain[-(1:burnin),1], type = "l", xlab="MCMC Steps" , main = "Chain values of x")
abline(h = median(chain[-(1:burnin),1]), col="red")