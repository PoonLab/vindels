# require(bbmle)
# 
# insertions <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/ins-nosep-all.csv",row.names=1, stringsAsFactors = F)
# 
# # #slips <- matrix(rep(rep(0,121),10), nrow=10,ncol=121)
# # slips <- c(rep(0,121000))
# # slips[sample(121000,100)] <- sample(3,100,replace = T) * 3
# counts <- nchar(insertions[insertions$Count!=0, "Seq"])
# counts <- counts[-167]
# counts <- c(counts, rep(0, sum(nchar(insertions$Vseq))))
# 
# # randomly shuffle all the entries
# rnd <- sample(length(counts),length(counts))
# counts <- counts[rnd]

source("~/vindels/2_within-host/utils.r")
path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

# PROBLEMATIC CASE: remove instances with gaps in the ancestor but NO INSERTION
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]

# CASE: remove instances with insertion position 0
insertions <- insertions[-c(which(insertions$Pos==0)),]
res <- as.data.frame(t(unname(mapply(restoreTipDel,insertions$Vseq, insertions$Anc, insertions$Seq, insertions$Pos))))
insertions$Vseq <- as.character(res[,1])
insertions$Pos <- as.numeric(as.character(res[,2]))

# #slips <- matrix(rep(rep(0,121),10), nrow=10,ncol=121)
# slips <- c(rep(0,121000))
# slips[sample(121000,100)] <- sample(3,100,replace = T) * 3
slips <- nchar(insertions[insertions$Count!=0, "Seq"])
#slips <- slips[-167]
slips <- c(slips, rep(0, sum(nchar(insertions$Vseq))))

# randomly shuffle all the entries
rnd <- sample(length(slips),length(slips))
slips <- slips[rnd]
counts <- slips
rm(slips)

# 
# objf <- function(p.slip, p.stay){
#   -affinell(p.slip, p.stay, slips)
# }
# 
# result <- mle2(objf, start=list(p.slip=1,p.stay=1), method = "L-BFGS-B", lower=1e-12, upper = 1)
# 
# 
# geomll <- function(p.copy){
#   N <- length(counts)
#   # log likelihood of the geometric distribution
#   N * log(p.copy) + sum(counts) * log(1-p.copy)
# }


# CUSTOM MCMC IN R


likelihood <- function(slip){
  if (slip <= 0){
    return(-Inf)
  }else if (slip > 1){
    return(-Inf)
  }
  #sum(dgeom(counts,prob=(1-slip), log=T))
  sum(dpois(counts,lambda=slip, log=T))
}

prior <- function(slip){
  prior <- dunif(slip, min=0, max=10, log = T)
  
  return(prior)
}


posterior <- function(slip){
  prior(slip) + likelihood(slip)
}

proposalFunction <- function(slip){
  return(rnorm(1,mean=slip, sd=0.01))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,2))
  chain[1,1] <- startvalue
  chain[1,2] <- posterior(startvalue)
  
  for (i in 1:iterations){
    p.current <- posterior(chain[i,1])
    
    proposal <- proposalFunction(chain[i,1])
    p.next <- posterior(proposal) 
    prop <- exp(p.next - p.current)
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- c(proposal, p.next)
      
    # if the proportion is less than the random uniform sample, REJCECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
    }
    if (i %% 100 == 0){
      print(paste(c("STATE ",i,": ", chain[i,]), collapse=" "))
    }
  }
  return(chain)
  
}

# RUN MCMC
startvalue <- 0.5
chain <- runMCMC(startvalue, 25000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- ceiling((nrow(chain)-1)*0.1)
acceptance <- 1 - mean(duplicated(chain[-(1:burnin),]))
print(paste0("Acceptance: ", acceptance))

#png(file="~/vindels/Figures/within-host/mcmc-posterior.png",width=800,height=600, res=120)
#hist(chain[-(1:burnin),1],nclass=30, main="Posterior of x", xlab="Probability of Copy (1 - slip)" )
#abline(v = median(chain[-(1:burnin),1]), col="red")
#dev.off()

# PLOTTING 
#hist(chain[-(1:burnin),1],nclass=30, main="Posterior of x", xlab="True value = red line" )
#abline(v = median(chain[-(1:burnin),1]), col="red")

#png(file="~/vindels/Figures/within-host/mcmc-trace.png",width=800,height=600, res=120)
#plot(chain[-(1:burnin),1], type = "l", xlab="MCMC Steps" , main = "Chain values of x")
#abline(h = median(chain[-(1:burnin),1]), col="red")
#dev.off()
