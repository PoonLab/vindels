source("~/vindels/2_within-host/utils.r")
source("~/vindels/2_within-host/slip-model-utils.r")
path <- "~/PycharmProjects/hiv-withinhost/"

# ---- Real Data ---- 
insertions <- read.csv(paste0(path,"10_nucleotide/tips/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

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



# ----- Simulated Data ----

genSeq <- function(len){
  seq <- c()
  for (n in 1:len){
    num <- runif(1)
    if (num < f[1]){
      nucl <- "A"
    }else if (num > f[1] && num < (f[2]+f[1])){
      nucl <- "C"
    }else if (num > (f[2]+f[1]) && num < (f[3]+f[2]+f[1])){
      nucl <- "G"
    }else{
      nucl <- "T"
    }
    seq[n] <- nucl
  }
  paste(seq, collapse="")
}


nt <- c("A", "C", "G", "T")

lens <- c(75, 126, 105,  90 , 33)
f <- c(0.4158261, 0.1649874, 0.1761463, 0.2423612 )
names(f) <- nt

simPair <- function(rate){
  vlen <- lens[sample(1:5, 1)]

  # ------ INDELS -------
  # determine the number of insertions that occur 
  counts <- rpois(vlen,lambda=rate)
  return(counts)
}

# CUSTOM MCMC IN R

simSeqs <- function(iter, rate){
  all.counts <- sapply(1:iter, function(n){
    if (n %% 1000 == 0 ){
      print(n)
    }
    pair <- simPair(rate)
    # VALUE 1 = Tip, VALUE 2 = Ancestor, VALUE 3 = Branch length
    return(pair)
  })
  
  return(unlist(all.counts))
}

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
  prior <- dunif(slip, min=0, max=1, log = T)
  
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


counts <- simSeqs(25000, 0.001)

# RUN MCMC
startvalue <- 0.5
chain <- runMCMC(startvalue, 100000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- ceiling((nrow(chain)-1)*0.1)
acceptance <- 1 - mean(duplicated(chain[-(1:burnin),]))
print(paste0("Acceptance: ", acceptance))


# ----- Plotting ---- 
par(mar=c(6,6,3,1.5),las=1)
hist(chain[-(1:burnin),1], 
     freq=F,
     col="dodgerblue",
     breaks=15,
     yaxt="n",
     ylab="",
     main="",
     xlab="Posterior of Lambda",cex.axis=1.2,cex.lab=1.5)
abline(v=0.001,col='black',lwd=4,lty=2)
lines(xy.coords(x=c(0,1), y=c(1,1)),col="red",lwd=2.5)
title(ylab="Density", line=4, cex.lab=1.5)
axis(2,labels=c("0","5e3", "1e4", "1.5e4"),
     at=c(0,5000,10000,15000),
     cex.axis=1.2)

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
