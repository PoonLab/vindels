require(bbmle)
require(ape)
# used to fill in deletion gaps found in the tip sequences 

source("~/vindels/2_within-host/utils.r")

createSlips <- function(anc, ins, pos){
  # start out with a base vector containing nchar number of zeros 
  # remove the gap characters from the ancestral sequence 
  anc <- gsub("-","",anc)
  base <- rep(0,nchar(anc)+1)
  # if there is no insertion, simply add a zero to complete the vector 
  if (ins == ""){
    return (base)
    
    # if there is an insertion, add the slip count  at the appropriate position
  }else{
    len <- nchar(ins)
    pos <- as.numeric(pos) - len + 1
    base[pos] <- len
    return(base)
  }
}


# calculate the median lengths of the variable loops 
ins.v <- split(insertions, insertions$Vloop)
lens <- unname(unlist(lapply(ins.v, function(x){median(x[,"Vlength"])})))
f <- estimateFreq(allseqs)
nucleotide <- function(x){
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
  nucl
}
simulateDNA <- function(p.enter, p.stay, rate, lambda){
  vlen <- lens[sample(1:5, 1)]
  
  seq <- sapply(1:vlen,nucleotide) 
  
  pois <- c()
  
  # this code runs with P(enter)
  if (runif(1) < p.enter){
    
    # get a poisson distributed value greater than 0 (0 is a waste of computation power)
    smpl <- 0
    while (smpl <= 0){
      smpl <- rpois(1,lambda=5)
    }
    
    if (runif(1) < 0.05){
      # add a non-multiple of 3 indel
    }else{

    }
  }
  
  
  # find all non-3 values
  non3 <- which(pois)
  
  # choose ~90% of these non-3 values to place in the "toRemove" vector
  toRemove <- which(sapply(1:length(non3), function(n){
    if(runif(1) < 0.05) return (F) else return(T)
  }))
  
  # the ~90% in the toRemove vector will be highlighted in the non3 vector,
  # which will be removed from the pois sample
  pois <- c(pois, smpl[-(non3[toRemove])])
  print(length(pois))
  
  
  # to estimate the ACTUAL rate of indels, you need to make failures occur after p.enter has been selected
  # algorithm:
    # probability of enter is chosen
    # draw a number from a poisson process
    # if the number is %% 3 == 0: 
      # keep it 100 percent
    # else if the number if %%3 != 0:
      # there's a low probability that it will be kept (penalty)
  
  # to estimate the OBSERVED rate of indels, it is a simpler process
  # algorithm:
    # probability of enter is chosen
    # draw a random number to determine whether a 
    # if number %% 3 == 0:
      # keep the number 
    # else if number %% 3 != 0:
      # while loop to redraw the number until 
  
  
  
  slips <- rep(0, nchar(seq)+1)
  
  
}

getSecond <- function(seq, p.enter, p.stay, lambda){
  
}

getTip <- function(oldtip, slip){
  nonzeros <- which(slip != 0)
  
  if (length(nonzeros) == 0){
    return(oldtip)
  }else{
    #print(oldtip)
    toCopy <- rep(T, nchar(oldtip))
    tip.chars <- strsplit(oldtip, "")[[1]]
    loc <- getSlipLocations(slip)[[1]]
    tab <- table(loc)
    for (n in 1:length(tab)){
      start <- as.numeric(names(tab)[n])
      stop <- start + (tab[[n]] - 1)
      # Case to catch when the 
      if (start > nchar(oldtip) || stop > nchar(oldtip)){
        start <- nchar(oldtip) - (stop - start)
        stop <- nchar(oldtip)
      }
      # re-adjustment of positions 
      # perform this action on all indices after the first
      if (n > 1 && !toCopy[start]){
        adjust <- min(which(toCopy[start:length(toCopy)])) - 1
        start <- start + adjust
        stop <- stop + adjust
      }
      #print(start)
      #print(stop)
      toCopy[start:stop] <- F
    }
    #tip.chars[toCopy] <- "-"
    return(paste0(tip.chars[toCopy],collapse=""))
  }
}

# used to go from c(0,0,0,3,0,0,0) to c(4,4,4)
getSlipLocations <- function(slip){
  nonzeros <- which(slip!=0)
  locations <- c()
  for (pos in nonzeros){
    locations <- c(locations, rep(pos, slip[pos]))
  }
  tab <- table(locations)
  return (list(loc=locations,len=length(slip)))
}

# used to go from c(4,4,4) to c(0,0,0,3,0,0,0) 
getSlipVector <- function(locs, length){
  vect <- rep(0,length)
  if (length(locs) == 0){
    return (vect)
  }else{
    tab <- table(locs)
    for (n in 1:length(tab)){
      vect[as.numeric(names(tab)[n])] <- unname(tab[n])
    }
    return(vect)
  }
}
delta <- function(rep=1,mean=0,sd=1.5){
  x <- rnorm(rep,mean=mean,sd=sd)
  if (x < 0){
    x <- abs(x)
    x <- ceiling(x)
    x <- -x
  }else{
    x <- ceiling(x)
  }
  x
}
position <- function(len){
  sample(len, 1)
}



path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

# FIX HEADERS
insertions$Header <- gsub("_\\d$","",insertions$Header)

# CASE: remove instances missing ancestor and tip 
insertions <- insertions[-c(which(insertions$Anc == "")),]

# CASE: remove instances with gaps in the ancestor but NO INSERTION
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]

# CASE: remove instances with insertion position 0
insertions <- insertions[-c(which(insertions$Pos==0)),]

# Restore all deletions found in tip sequences and adjust the POS values accordinaly
res <- as.data.frame(t(unname(mapply(restoreTipDel,insertions$Vseq, insertions$Anc, insertions$Seq, insertions$Pos))))
insertions$Vseq <- as.character(res[,1])
insertions$Pos <- as.numeric(as.character(res[,2]))

# CASE: restore all gaps from OTHER insertions
insertions$Anc <- mapply(removeOtherGaps, insertions$Anc, insertions$Vseq, insertions$Seq, insertions$Pos)

# SANITY CHECK: to make sure all seqs are equal
a <- nchar(insertions$Vseq) - nchar(insertions$Seq)
b <- nchar(gsub("-","",insertions$Anc))
sum(a!=b)==0
insertions[which(a!=b),]

# CASE: replace "R" nucleotides with the corresponding one found in the ancestor
cases <- which(grepl("[RYSWKMBDHVN]", insertions$Vseq))
for (idx in cases){
  toEdit <- insertions[idx, "Vseq"]
  pos <- gregexpr("[RYSWKMBDHVN]", toEdit)[[1]]
  toUse <- insertions[idx,"Anc"]
  insertions[idx, "Vseq"] <- paste0(substr(toEdit, 1, pos-1),substr(toUse,pos,pos), substr(toEdit, pos+1,nchar(toEdit)))
}

# generate slip list 
slip.list <- unname(mapply(createSlips, insertions$Anc, insertions$Seq, insertions$Pos))

# add the a/b replicate label to the headers
#insertions$Header <- unname(mapply(patLabel, insertions$Header, insertions$Pat))
names(slip.list) <- insertions$Header

# # for PRINTING the slip.list 
# for (elem in 1:length(slip.list)){
#   if(sum(slip.list[[elem]]) > 0){
#     print(slip.list[elem])
#   }
# }
# # C.-.-.QT.10R.-.-_289_1_b
# # SHUFFLING --- randomly shuffle the slip locations around 
# slip.list <- lapply(slip.list, function(x){
#   total <- sum(x)
#   locs <- sample(length(x), total, replace=T)
#   getSlipVector(locs, length(x))
# })

nucleotides <- c("A","C","G","T")
idx <- which(unname(lapply(slip.list,sum))>0)
# For use in the proposal function
changeSlip <- function(){
  # choose a sequence to edit
  seq <- sample(length(idx),1)
  
  # convert it to indices
  slip <- slip.list[[idx[seq]]]
  slip.idx <- getSlipLocations(slip)
  
  # choose a slip event to change
  toEdit <- position(length(slip.idx[[1]]))
  
  # this is to ensure that the proposed change is never outside the slip region
  proposal <- 0
  while(proposal <= 0 || proposal > length(slip)){
    proposal <- slip.idx[[1]][toEdit] + delta()
  }
  slip.idx[[1]][toEdit] <- proposal
  # save this globally so that the change gets fixed
  slip.list[[idx[seq]]] <<- getSlipVector(slip.idx[[1]],slip.idx[[2]])
}
nt <- c("A", "C", "G", "T")
require(stringr)
estimateFreq <- function(seqs){
  nt <- c("A", "C", "G", "T")
  output <- c()
  for (n in 1:length(nt)){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
}
require(expm)
# returns the F81 transition probability matrtix 
getMat <- function(rate, branch){
  
  nt <- c("A", "C", "G", "T")
  
  # generate the F81 rate matrix 
  mat <- matrix(rep(f, each=4), nrow=4, ncol=4,dimnames=list(nt,nt))
  mat <- mat * rate
  diag(mat) <- sapply(1:4, function(x) -(sum(rate * f[-x]))) # entirely equivalent to -rate*(sum(f[-x]))
  
  # multiply by branch length
  mat <- branch * mat
  
  # exponentiate and return
  tmat <- expm(mat)
  tmat
}
nt <- c("A", "C", "G", "T")
pairllh <- function(anc, newtip, rate, branch){
  tmat <- getMat(rate,branch)
  achars <- strsplit(anc, "")[[1]]
  tchars <- strsplit(newtip, "")[[1]]
  print(branch)
  result <- mapply(function(tchar,achar){
    # initializes a matrix with a 
    tip.llh <- matrix(as.numeric(tchar == nt), nrow=4,ncol=1,dimnames=list(nt))
    
    # finalize the calculation for tip likelihood
    # dot product
    llh <- tmat %*% tip.llh
    llh <- llh * f
    if (!achar %in% nt){
      print(achar)
    }
    # likelihood given the exact nucleotide (state) that we see in the ancestor
    final.llh <- llh[achar,]
    return(log(final.llh))
  }, achars, tchars)
  
  # if (class(result) == "list"){
  #   print(anc)
  #   print(newtip)
  # }
  sum(result)
}

nt <- c("A", "C", "G", "T")
allseqs <- c(insertions$Vseq, insertions$Anc)
f <- estimateFreq(allseqs)

branches <- insertions$length
anc.seqs <- gsub("-", "", insertions$Anc)
require(parallel)
seqllh <- function(rate){
  #print("Starting SeqLLH ... ")
  tip.seqs <- unname(mcmapply(getTip, insertions$Vseq, slip.list, mc.cores=16))
  #print(head(tip.seqs))
  rate <- rep(rate, length(anc.seqs))
  #print("Calculating pairwise ... ")
  total.llh <- unname(mcmapply(pairllh, anc.seqs, tip.seqs, rate, branches, mc.cores=16))
  #print(head(total.llh))
  sum(total.llh)
}


likelihood <- function(param){
  #print("Starting LLH...")
  p.enter <- param[1]
  p.stay  <- param[2]
  rate <- param[3]
  
  slips <- unname(unlist(slip.list))
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  if (any(is.na(slips))){
    llh <- log(0)

  }else{
    llh <-  x*log(1-p.enter) + y*log(p.enter) + y*log(1-p.stay) + z*log(p.stay)
  }

  llh + seqllh(rate)
}
# --------------------------------------------------

prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  rate <- param[3]
  
  prior.pe <- dlnorm(p.enter,meanlog=-8,sdlog=2, log=T)
  prior.ps <- dlnorm(p.stay,meanlog=-0.2,sdlog=0.1,log=T)
  prior.rate <- dlnorm(p.stay,meanlog=-7,sdlog=2,log=T)
  
  return(prior.pe + prior.ps + prior.rate)
}

posterior <- function(param){
  #print(prior(param))
  #print(likelihood(param))
  if (any(param > 1)){
    return(log(0))
  }else{
    return(prior(param) + likelihood(param))
  }
}

proposalFunction <- function(param){
  p.enter <- param[1]
  p.stay <- param[2]
  rate <- param[3]
  
  num <- runif(1)
  if (num > 0.9){
    if (num-0.9 < 1/30){
      p.enter <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
    }else if(num-0.9 > 2/30){
      p.stay <- rlnorm(1,meanlog=log(param[2]),sdlog=0.02)
    }else{
      rate <- rlnorm(1,meanlog=log(param[3]),sdlog=0.08)
    }
  }else{
    # perform a change on the sliplist 
    # this will NOT change your three parameters, but WILL change the likelihood
    changeSlip()
  }
  return(c(p.enter, p.stay, rate))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,3))
  chain[1,] <- startvalue
  
  for (i in 1:iterations){
    proposal <- proposalFunction(chain[i,])
    p.current <- posterior(chain[i,])
    p.next <- posterior(proposal)
    print(p.current)
    print(p.next)
    if(is.na(p.current) || is.na(p.next)){
      print("ERROR: Posterior could not be calculated")
      print(paste("Chain value:", chain[i,1], chain[i,2], chain[i,3], sep=" "))
      break
    }
    prop <- exp(p.next - p.current)
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- proposal
      
    # if the proportion is less than the random uniform sample, REJECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
    }
    print(paste("STATE",i,":", chain[i,1], chain[i,2], chain[i,3], sep=" "))
    
  }
  return(chain)
}

# RUN MCMC
startvalue <- c(0.001, 0.8, 0.001)
chain <- runMCMC(startvalue, 10000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- 0.1*(nrow(chain)-1)
acceptance <- 1 - mean(duplicated(chain[-(1:burnin),]))
print(paste0("Acceptance: ", acceptance))


med1 <- round(median(chain[-(1:burnin),1]),4)
med2 <- round(median(chain[-(1:burnin),2]),4)
med3 <- round(median(chain[-(1:burnin),3]),4)
par(mfrow=c(1,3), mar=c(5,5,4,1))
# PLOTTING 
hist(chain[-(1:burnin),1],nclass=30, main="Posterior of Enter", xlab="Prob(Enter)",ylab="Frequency",col="lightskyblue")
abline(v = med1, col='red',lwd=2)
text(0.000168, 70000, paste0("Median = ", med1))
text(0.000168, 60000, paste0("Acceptance = ", as.character(acceptance)))
hist(chain[-(1:burnin),2],nclass=30, main="Posterior of Stay", xlab="Prob(Stay)", ylab="Frequency",col="lightskyblue")
text(0.905, 70000, paste0("Median = ", med2))
abline(v = med2, col='red',lwd=2)
hist(chain[-(1:burnin),3],nclass=30, main="Posterior of Rate", xlab="Rate", ylab="Frequency",col="lightskyblue")
#text(0.00017, 400, paste0("Acceptance = ",as.character(acceptance)))
dev.off()

len <- length(chain[,1])
thinned <- seq(burnin, len, 100)

par(mfrow=c(1,2), mar=c(5,5,4,1))
plot(chain[thinned,1], type = "l", xlab="MCMC Steps" , ylab="Prob(Enter)",main = "Chain values of Enter")
abline(h = med1, col="red")
plot(chain[thinned,2], type = "l", xlab="MCMC Steps" , ylab="Prob(Stay)",main = "Chain values of Stay")
abline(h = med2, col="red")

# get slip locations 

# randomly sample a location on the slip locations 

# change that value by a rnorm derived number

# get slip vector again 

# generate the newtip sequence from this slip vector

# calculate the likelihood of the newtip sequence given the ancestor
# this requres the transition probability matrix 
# remember: matrix exponentiation 

# transition matrix 

# receive the newtip and the ancestor for input 
# first check that they are a) the same length, b) contain no gaps 
# modify the checkDiff function or the transition function to generate a matrix of transition probs
# start by getting counts at each location 


