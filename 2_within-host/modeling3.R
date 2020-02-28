source("~/vindels/2_within-host/utils.r")
path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

# PROBLEMATIC CASE: remove instances with gaps in the ancestor but NO INSERTION
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]

# CASE: remove instances with insertion position 0
insertions <- insertions[-c(which(insertions$Pos==0)),]
res <- as.data.frame(t(unname(mapply(restoreDel,insertions$Vseq, insertions$Anc, insertions$Seq, insertions$Pos))))
insertions$Vseq <- as.character(res[,1])
insertions$Pos <- as.numeric(as.character(res[,2]))

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

# attempted to find the log - likelihood of an affine model 
likelihood <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  # Log likelihood of each tip/anc pair
  #(1 - p.slip)^x * p.slip^y * (1-p.stay)^y * p.stay^z
  llh <-  x*log(1-p.enter) + y*log(p.enter) + y*log(1-p.stay) + z*log(p.stay)
  llh
  
}
# --------------------------------------------------

prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  prior.pe <- dlnorm(p.enter,meanlog=-10,sdlog=2, log=T)
  prior.ps <- dlnorm(p.stay,meanlog=-0.15,sdlog=0.05,log=T)
  
  return(prior.pe + prior.ps)
}

posterior <- function(param){
  #print(prior(param))
  #print(likelihood(param))
  prior(param) + likelihood(param)
}

proposalFunction <- function(param){
  p.enter <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
  p.stay <- rlnorm(1,meanlog=log(param[2]),sdlog=0.02)
  return(c(p.enter,p.stay))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,2))
  chain[1,] <- startvalue
  for (i in 1:iterations){
    
    proposal <- proposalFunction(chain[i,])
    #print(proposal)
    #print(posterior(proposal))
    #print(posterior(chain[i,]))
    prop <- exp(posterior(proposal) - posterior(chain[i,]))
    #print(prop)
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
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
startvalue <- c(1e-5,0.7)
chain <- runMCMC(startvalue, 1000000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- 100000
acceptance <- round(1 - mean(duplicated(chain[-(1:burnin),])),2)


#chain <- read.csv("~/slip-model-parsed.csv", row.names = 1, header=F)


med1 <- round(median(chain[-(1:burnin),1]),6)
med2 <- round(median(chain[-(1:burnin),2]),3)

par(mfrow=c(1,2), mar=c(5,5,4,1))
# PLOTTING 
hist(chain[-(1:burnin),1],nclass=30, main="Posterior of Enter", xlab="Prob(Enter)",ylab="Frequency",col="lightskyblue")
abline(v = med1, col='red',lwd=2)
text(0.000168, 70000, paste0("Median = ", med1))
text(0.000168, 60000, paste0("Acceptance = ", as.character(acceptance)))
hist(chain[-(1:burnin),2],nclass=30, main="Posterior of Stay", xlab="Prob(Stay)", ylab="Frequency",col="lightskyblue",xlim=c(0.87,0.91))
text(0.905, 70000, paste0("Median = ", med2))
abline(v = med2, col='red',lwd=2)
#text(0.00017, 400, paste0("Acceptance = ",as.character(acceptance)))
dev.off()

len <- length(chain[,1])
thinned <- seq(burnin, len, 100)

par(mfrow=c(1,2), mar=c(5,5,4,1))
plot(chain[thinned,1], type = "l", xlab="MCMC Steps" , ylab="Prob(Enter)",main = "Chain values of Enter")
abline(h = med1, col="red")
plot(chain[thinned,2], type = "l", xlab="MCMC Steps" , ylab="Prob(Stay)",main = "Chain values of Stay")
abline(h = med2, col="red")