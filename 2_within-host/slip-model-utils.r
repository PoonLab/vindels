# SLIP MODEL UTILS

createSlips <- function(anc, len, pos){
  # start out with a base vector containing nchar number of zeros 
  # remove the gap characters from the ancestral sequence 
  anc <- gsub("-","",anc)
  base <- rep(0,nchar(anc))  # removed the nchar(anc) + 1 because there cannot be a slip at the final position (there is nothing to skip over)
  # if there is no insertion, simply add a zero to complete the vector 
  if (len == 0 || is.na(len)){
    return (base)
    
    # if there is an insertion, add the slip count  at the appropriate position
  }else{
    pos <- as.numeric(pos) - len + 1
    base[pos] <- len
    return (base)
  }
}

estimateFreq <- function(seqs){
  require(stringr)
  nt <- c("A", "C", "G", "T")
  output <- c()
  for (n in 1:length(nt)){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
}

getSlipVector <- function(locs, length){
  # used to go from c(4,4,4) to c(0,0,0,3,0,0,0) 
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

getTip <- function(oldtip, slip){
  # used to determine the new tip sequence using the slip index
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

getSlipLocations <- function(slip){
  # used to go from c(0,0,0,3,0,0,0) to c(4,4,4)
  nonzeros <- which(slip!=0)
  locations <- c()
  for (pos in nonzeros){
    locations <- c(locations, rep(pos, slip[pos]))
  }
  tab <- table(locations)
  return (list(loc=locations,len=length(slip)))
}

delta <- function(sd=3){
  # chooses a normally distributed discrete value above 0
  x <- rnorm(1,mean=0,sd=sd)
  if (x < 0){
    x <- abs(x)
    x <- ceiling(x)
    x <- -x
  }else{
    x <- ceiling(x)
  }
  x
}

getMat <- function(rate, branch){
  # returns the F81 transition probability matrtix 
  require(expm)
  nt <- c("A", "C", "G", "T")
  
  # generate the F81 rate matrix 
  mat <- matrix(rep(f, each=4), nrow=4, ncol=4,dimnames=list(nt,nt))
  inv.freq <- sapply(1:4, function(x) 1-f[x])
  mat <- mat * rate
  #diag(mat) <- sapply(1:4, function(x) -(rate*sum(f[-x]))) # equivalent to -rate*(sum(f[-x]))
  diag(mat) <- sapply(1:4, function(x) -(rate*inv.freq[x]))
  # multiply by branch length
  mat <- branch * mat
  
  # exponentiate and return
  tmat <- expm(mat)
  tmat
}

changeSlip <- function(slip.list){
  slip.idx <- getSlipLocations(slip.list)
  
  tab <- table(slip.idx[[1]])
  
  # choose a slip event to change
  toEdit <- sample(length(slip.idx[[1]]),1)
  
  stack <- tab[as.character(slip.idx[[1]][toEdit])]
  if (stack > 1){   # have a set probability of moving the entire stack instead of individual events 
    
    all.idx <- which(slip.idx[[1]] == as.numeric(names(stack[1])))
    rnum <- runif(1)
    if (rnum < 0.9){
      # rewrite the toEdit variable with all indices in the stack
      toEdit <- all.idx 
    }else if(rnum > 0.9 && rnum < 0.95){
      if (runif(1) < 0.5){
        toEdit <- sample(all.idx, floor(length(all.idx) * 0.5))
      }else{
        toEdit <- sample(all.idx, ceiling(length(all.idx) * 0.5))
      }
    } 
  }
  
  # this is to ensure that the proposed change is never outside the slip region
  proposal <- slip.idx[[1]][toEdit]
  while(any(proposal) <= 0 || any(proposal) > length(slip.list)){
    proposal <- slip.idx[[1]][toEdit] + delta()
  }
  # save the change to the slip list
  slip.idx[[1]][toEdit] <- proposal
  
  # save the whole list
  getSlipVector(slip.idx[[1]],slip.idx[[2]])
}

pairllh <- function(anc, newtip, rate, branch){
  tmat <- getMat(rate,branch)
  achars <- strsplit(anc, "")[[1]]
  tchars <- strsplit(newtip, "")[[1]]
  #print(branch)
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
  
  sum(result)
}

seqllh <- function(rate, slip.list){
  # likelihood of the entire slip.list
  # HIGHEST TIME COMPLEXITY -- only calculated when : 
  # a) rate parameter is changed
  # b) before the MCMC starts for the first iteration
  
  require(parallel)
  #print("Starting SeqLLH ... ")
  tip.seqs <- unname(mcmapply(getTip, indels$tip, slip.list, mc.cores=16))
  #print(head(tip.seqs))
  rate <- rep(rate, length(anc.seqs))
  #print("Calculating pairwise ... ")
  unname(mcmapply(pairllh, anc.seqs, tip.seqs, rate, branches, mc.cores=16))
  #print(head(total.llh))
}