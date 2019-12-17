require(bbmle)

insertions <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/ins-nosep-all.csv",row.names=1)

#slips <- matrix(rep(rep(0,121),10), nrow=10,ncol=121)
slips <- c(rep(0,121000))
slips[sample(121000,100)] <- sample(3,100,replace = T) * 3


objf <- function(p.slip, p.stay){
  -affinell(p.slip, p.stay, slips)
}

result <- mle2(objf, start=list(p.slip=1,p.stay=1), method = "L-BFGS-B", lower=1e-12, upper = 1)


geomll <- function(p.copy){
  N <- length(counts)
  # log likelihood of the geometric distribution
  N * log(p.copy) + sum(counts) * log(1-p.copy)
}


# CUSTOM MCMC IN R


likelihood <- function(slip){
  if (slip <= 0){
    return(0)
  }else if (slip > 1){
    return(0)
  }
  sum(dgeom(counts,prob=slip, log=T))
}

prior <- function(slip){
  prior <- dunif(slip, log = T)
  
  return(prior)
}


posterior <- function(slip){
  prior(slip) + likelihood(slip)
}

proposalFunction <- function(slip){
  return(rnorm(1,mean=slip, sd=0.003))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,1))
  chain[1,] <- startvalue
  
  
  for (i in 1:iterations){

    proposal <- proposalFunction(chain[i,])
    prop <- exp(posterior(proposal) - posterior(chain[i,]))
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
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