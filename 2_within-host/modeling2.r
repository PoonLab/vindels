require(bbmle)
require(ape)
# used to fill in deletion gaps found in the tip sequences 
nt <- c("A", "C", "G", "T")
require(stringr)
source("~/vindels/2_within-host/utils.r")

estimateFreq <- function(seqs){
  nt <- c("A", "C", "G", "T")
  output <- c()
  for (n in 1:length(nt)){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
}
createSlips <- function(anc, ins, pos){
  # start out with a base vector containing nchar number of zeros 
  # remove the gap characters from the ancestral sequence 
  anc <- gsub("-","",anc)
  base <- rep(0,nchar(anc))  # removed the nchar(anc) + 1 because there cannot be a slip at the final position (there is nothing to skip over)
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
sum(a!=b)==0              # should be TRUE
insertions[which(a!=b),]  # should be nrow = 0

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

# # C.-.-.QT.10R.-.-_289_1_b
# # SHUFFLING --- randomly shuffle the slip locations around 
slip.list <- lapply(slip.list, function(x){
  total <- sum(x)
  locs <- sample(length(x), total, replace=T)
  getSlipVector(locs, length(x))
})


getTip <- function(oldtip, slip){
  nonzeros <- which(slip != 0)
  
  if (length(nonzeros) == 0){
    return(oldtip)
  }else{
    #print(oldtip)
    toCopy <- rep(T, nchar(oldtip))logfile <- 
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
collapseVect <- function(vect){
  nonzero <- which(vect != 0)
  
  for (i in length(nonzero):2){
    # check whether there is overlap in these slip positions, if so, amalgamate
    current <- nonzero[i]
    previous <- nonzero[i-1]
    
    # checks whether the previous nonzero is adjacent
    if ((current - 1) == previous){
      vect[previous] <- vect[previous] + vect[current]
      vect[current] <- 0
    }
  }
  return(vect)
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
    return(collapseVect(vect))
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

idx <- which(unname(lapply(slip.list,sum))>0)
# For use in the proposal function
changeSlip <- function(slip.list){
  # choose a sequence to edit
  seq <- sample(length(idx),1)
  
  # convert it to indices
  slip <- slip.list[[idx[seq]]]
  slip.idx <- getSlipLocations(slip)
  
  # choose a slip event to change
  toEdit <- sample(length(slip.idx[[1]]),1)
  
  # this is to ensure that the proposed change is never outside the slip region
  proposal <- 0
  while(proposal <= 0 || proposal > length(slip)){
    proposal <- slip.idx[[1]][toEdit] + delta()
  }
  # save the change to the slip list
  slip.idx[[1]][toEdit] <- proposal
  
  # save the whole list
  slip.list[[idx[seq]]] <- getSlipVector(slip.idx[[1]],slip.idx[[2]])
  slip.list
}


require(expm)
# returns the F81 transition probability matrtix 
getMat <- function(rate, branch){
  
  nt <- c("A", "C", "G", "T")
  
  # generate the F81 rate matrix 
  mat <- matrix(rep(f, each=4), nrow=4, ncol=4,dimnames=list(nt,nt))
  mat <- mat * rate
  diag(mat) <- sapply(1:4, function(x) -(rate*sum(f[-x]))) # equivalent to -rate*(sum(f[-x]))
  
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
# NORMAL DATA: 
f <- estimateFreq(c(insertions$Vseq, insertions$Anc))
branches <- insertions$length
anc.seqs <- gsub("-", "", insertions$Anc)

# SIMULATED DATA: 

require(parallel)
seqllh <- function(rate, slip.list){
  #print("Starting SeqLLH ... ")
  tip.seqs <- unname(mcmapply(getTip, insertions$Vseq, slip.list, mc.cores=16))
  #print(head(tip.seqs))
  rate <- rep(rate, length(anc.seqs))
  #print("Calculating pairwise ... ")
  total.llh <- unname(mcmapply(pairllh, anc.seqs, tip.seqs, rate, branches, mc.cores=16))
  #print(head(total.llh))
  sum(total.llh)
}


likelihood<- function(param, slip.list){
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
  llh + seqllh(rate, slip.list)
}


# --------------------------------------------------

prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  rate <- param[3]
  
  prior.pe <- dlnorm(p.enter,meanlog=-8,sdlog=0.5, log=T)
  prior.ps <-  dunif(p.stay, min=0.5, max=1.0, log=T) # dlnorm(p.stay,meanlog=-0.17,sdlog=0.05,log=T)
  prior.rate <- dunif(rate, min=0, max=0.01, log=T )#dlnorm(rate,meanlog=log(erate),sdlog=0.3,log=T)
  
  return(prior.pe + prior.ps + prior.rate)
}

posterior <- function(param, slip){
  #print(prior(param))
  #print(likelihood(param))
  if (any(param > 1)){
    return(log(0))
  }else{
    return(prior(param) + likelihood(param, slip))
  }
}

proposalFunction <- function(param, slip_current){
  p.enter <- param[1]
  p.stay <- param[2]
  rate <- param[3]
  
  num <- runif(1)
  s2p <- 0.85
  if (num > s2p){
    if (num-s2p < (1/3 *(1-s2p))){
      p.enter <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
    }else if(num-s2p > (2/3 *(1-s2p))){
      p.stay <- rlnorm(1,meanlog=log(param[2]),sdlog=0.02)
    }else{
      rate <- rlnorm(1,meanlog=log(param[3]),sdlog=0.08)
    }
    slip_proposed <- slip_current
  }else{
    # perform a change on the sliplist 
    # this will NOT change your three parameters, but WILL change the likelihood
    slip_proposed <- changeSlip(slip_current)
  }
  return(list(param=c(p.enter, p.stay, rate), slip=slip_proposed))
}


runMCMC <- function(startvalue, iterations, slip.list){
  chain <- array(dim = c(iterations+1,3))
  chain[1,] <- startvalue
  slip_current <- slip.list
  logfile <- file("~/PycharmProjects/hiv-withinhost/slip-model.csv", "w")
  write("p(Enter), p(Stay), Rate, Slip-changed, Accept", 
      file=logfile)
  for (i in 1:iterations){
    # calculate posterior of current position
    p.current <- posterior(chain[i,], slip_current)
    
    # generate proposal 
    proposal <- proposalFunction(chain[i,], slip_current)
    param_proposed <- proposal[[1]]
    slip_proposed <- proposal[[2]]
    
    # calculate posterior of the new proposal 
    p.next <- posterior(param_proposed, slip_proposed)
    
    print(p.current)
    print(p.next)
    # print(paste("PROPOSAL",i,":", param_proposed, sep=" "))
    # print(paste0("Sliplist change proposed: ", any(unname(unlist(slip_current))!=unname(unlist(slip_proposed)))))
    s.change <- any(unname(unlist(slip_current))!=unname(unlist(slip_proposed)))
    # to catch problematic posterior calculations 
    if(is.na(p.current) || is.na(p.next)){
      print("ERROR: Posterior could not be calculated")
      print(paste("Chain value:", chain[i,1], chain[i,2], chain[i,3], sep=" "))
      break
    }
    
    prop <- exp(p.next - p.current)
    print(prop)
    
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- param_proposed
      slip_current <- slip_proposed
      #print("Accept")
      accept <- T
    # if the proportion is less than the random uniform sample, REJECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
      #print("Reject")
      accept <- F
    }
    print(paste("STATE",i,":", chain[i,1], chain[i,2], chain[i,3], sep=" "))
    if (i %% 10 == 0){
      write(paste(c(chain[i,], as.numeric(s.change), as.numeric(accept)),collapse=",") , file=logfile, append=T)
    }
  }
  return(list(chain=chain, slip=slip_current))
  close(logfile)
}

# RUN MCMC
startvalue <- c(0.001, 0.8, 0.001)
chain <- runMCMC(startvalue, 50000, slip.list)



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


