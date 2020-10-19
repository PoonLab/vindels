# SLIP MODEL UTILS
require(expm)
require(stringr)
require(parallel)

# Function to convert ancestral sequences into vectors of slip events
createSlips <- function(anc, len, pos){
  # start out with a base vector containing nchar number of zeros 
  base <- rep(0,nchar(anc)) 
  # remove the gap characters from the ancestral sequence 
  anc <- gsub("-","",anc)

  # if there is no insertion, return the vector of zeros 
  if (len == 0 || is.na(len)){
    return (base)
    
  # if there is an insertion, include the slip count at the appropriate position
  }else{
    pos <- as.numeric(pos) - (len-1)   # pos marks the end of the insertion; adjust by subtracting length-1 
    base[pos] <- len
    return (base)
  }
}

# Receives a vector of sequences and returns a vector of c(A,C,G,T) nucleotide proportions 
# found in the data
estimateFreq <- function(seqs){
  
  nt <- c("A", "C", "G", "T")
  
  output <- c()
  for (n in 1:4){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
}

# Function to convert vector of slip locations c(4,4,4,4,4,4)
# to slip vectors c(0,0,0,6,0,0)
getSlipVector <- function(locs, length){
  
  vect <- rep(0,length)
  if (length(locs) == 0){
    return (vect)
  }else{
    tab <- table(locs)
    vect <- replace(vect, as.numeric(names(tab)), as.numeric(tab))
    return(vect)
  }
}

# Function to convert vector of slip vectors c(0,0,0,6,0,0) 
# to slip locations c(4,4,4,4,4,4)
getSlipLocations <- function(slip){
  # used to go from c(0,0,0,3,0,0,0) to c(4,4,4)
  nonzeros <- which(slip!=0)
  if (length(nonzeros) == 0){
    return (list(loc=c(),len=length(slip)))
  }else{
    return (list(loc=rep(nonzeros, slip[nonzeros]),len=length(slip)))
  }
}

# Function to determine the new tip sequence from the old one using the slip vector
getTip <- function(oldtip, slip){
  
  nonzeros <- which(slip != 0)
  
  # no slip events
  if (length(nonzeros) == 0){
    return(oldtip)
  
  # one or more slip events
  }else{
    toCopy <- rep(T, nchar(oldtip))     # default vector of T values 
    tip.chars <- strsplit(oldtip, "")[[1]]   # nucleotide sequence (split up)
    loc <- getSlipLocations(slip)[[1]]     # returns the locations of all slip events 
    tab <- table(loc)
    
    # for every distinct slip event specified in the slip vector, 
    for (n in 1:length(tab)){
      start <- as.numeric(names(tab)[n])
      stop <- start + (tab[[n]] - 1)
      
      
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

# Function for drawing discrete  distribution 0
delta <- function(sd=2){
  
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
  
  # set up the matrix
  mat <- matrix(rep(f, each=4), nrow=4, ncol=4,dimnames=list(nt,nt))
  # save the inverse of each frequency
  inv.freq <- sapply(1:4, function(x) 1-f[x])
  # multiply by rate 
  mat <- mat * rate
  # adjust the diagonals so that every row sums to 0
  diag(mat) <- sapply(1:4, function(x) -(rate*inv.freq[x]))
  # multiply by branch length
  mat <- branch * mat
  
  # exponentiate and return
  expm(mat)
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
  proposal <- slip.idx[[1]][toEdit] + delta(0.5)
  while(any(proposal <= 0) || any(proposal > length(slip.list))){
    proposal <- slip.idx[[1]][toEdit] + delta(0.5)
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

  result <- mapply(function(tchar,achar){
    # initializes simple row vector with 1 for the given nucleotide, 0 for the others
    tip.llh <- matrix(as.numeric(tchar == nt), nrow=4,ncol=1,dimnames=list(nt))
    
    # finalize the calculation for tip likelihood
    # dot product
    llh <- tmat %*% tip.llh    # tip.llh = c(1,0,0,0) for A
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

# likelihood of the entire slip.list
# HIGHEST TIME COMPLEXITY -- only calculated when : 
# a) rate parameter is changed
# b) before the MCMC starts for the first iteration

seqllh <- function(rate, slip.list){

  #print("Starting SeqLLH ... ")
  tip.seqs <- unname(mcmapply(getTip, indels$tip, slip.list, mc.cores=16))
  #print(head(tip.seqs))
  rate <- rep(rate, length(anc.seqs))
  #print("Calculating pairwise ... ")
  unname(mcmapply(pairllh, anc.seqs, tip.seqs, rate, branches, mc.cores=16))
  #print(head(total.llh))
}

# this function uses the START positions of insertions, not the end positions
insert <- function(str, indel,pos ){
  #pos <- pos+1
  vect <- strsplit(str,"")[[1]]
  if (pos < 1){
    return (NA)
  }
  
  if (pos == 1){
    return (paste(c(indel, vect[1:length(vect)]),collapse=""))
    
    # adds insert after the vector (at the end )
  }else if (pos == length(vect)+1){
    return (paste(c(vect[1:length(vect)], indel),collapse=""))
    # involves slicing the before and after parts of the vector around the insert
    # used when pos = 2 : nchar-1
  }else{
    return (paste(c(vect[1:pos-1], indel, vect[pos:length(vect)]),collapse=""))
  }
}
delete <- function(str, indel, pos){
  vect <- rep(T, nchar(str))
  
  if (pos > (nchar(str))){
    return (NA)
  }
  
  start <- pos - nchar(indel) + 1
  end <- pos
  
  vect[start:end] <- F
  
  chars <- strsplit(str, "")[[1]]
  
  return(paste(chars[which(vect)],collapse=""))
}
