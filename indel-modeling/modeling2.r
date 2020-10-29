require(bbmle)
require(ape)
# used to fill in deletion gaps found in the tip sequences 
nt <- c("A", "C", "G", "T")
require(stringr)
source("~/vindels/2_within-host/utils.r")


createSlips <- function(anc, len, pos){
  # start out with a base vector containing nchar number of zeros 
  # remove the gap characters from the ancestral sequence 
  anc <- gsub("-","",anc)
  base <- rep(0,nchar(anc))  # removed the nchar(anc) + 1 because there cannot be a slip at the final position (there is nothing to skip over)
  # if there is no insertion, simply add a zero to complete the vector 
  if (len == 0){
    return (base)
    
    # if there is an insertion, add the slip count  at the appropriate position
  }else{
    pos <- as.numeric(pos) - len + 1
    base[pos] <- len
    return (base)
  }
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
    return(collapseVect(vect))
  }
}

collapseVect <- function(vect){
  nonzero <- which(vect != 0)
  
  if (length(nonzero) < 2){
    return (vect)
  }else{
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
}
# ------  DATA IMPORT AND MANIPULATION -----

path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/all/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

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
slip.list <- unname(mapply(createSlips, insertions$Anc, nchar(insertions$Seq), insertions$Pos))

# add the a/b replicate label to the headers
#insertions$Header <- unname(mapply(patLabel, insertions$Header, insertions$Pat))
names(slip.list) <- insertions$Header

# # C.-.-.QT.10R.-.-_289_1_b
# # SHUFFLING --- randomly shuffle the slip locations around 
slip.list <- lapply(slip.list, function(x){
  total <- sum(x)
  if (total == 0){
    return (x)
  }else{
    locs <- sample(length(x), total, replace=T)
    getSlipVector(locs, length(x))
  }
})


# ------- PREAMBLE FUNCTIONS -------
estimateFreq <- function(seqs){
  nt <- c("A", "C", "G", "T")
  output <- c()
  for (n in 1:length(nt)){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
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





delta <- function(rep=1,mean=0,sd=1.5){
  # chooses a normally distributed discrete value above 0
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
  rand <- sample(length(idx),1)
  
  # convert it to indices
  slip <- slip.list[[idx[rand]]]
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
  slip.list[[idx[rand]]] <- getSlipVector(slip.idx[[1]],slip.idx[[2]])
  return(list(slip=slip.list, idx=idx[rand]))
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

# NORMAL DATA: 
f <- estimateFreq(c(insertions$Vseq, insertions$Anc))
branches <- insertions$Date
anc.seqs <- gsub("-", "", insertions$Anc)

# likelihood of the entire slip.list
# only calculated when : 
  # a) rate parameter is changed
  # b) before the MCMC starts for the first iteration
require(parallel)

# ------ MCMC FUNCTIONS ------
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
  
  # if (class(result) == "list"){
  #   print(anc)
  #   print(newtip)
  # }
  sum(result)
}

seqllh <- function(rate, slip.list){
  #print("Starting SeqLLH ... ")
  tip.seqs <- unname(mcmapply(getTip, insertions$Vseq, slip.list, mc.cores=16))
  #print(head(tip.seqs))
  rate <- rep(rate, length(anc.seqs))
  #print("Calculating pairwise ... ")
  unname(mcmapply(pairllh, anc.seqs, tip.seqs, rate, branches, mc.cores=16))
  #print(head(total.llh))
}

likelihood<- function(param, slip.list, llh.list){
  #print("Starting LLH...")
  p.enter <- param[1]
  p.stay  <- param[2]
  
  # affine gap likelihood
  # utilizes both p.enter + p.stay
  slips <- unname(unlist(slip.list))
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  if (any(is.na(slips))){
    llh <- log(0)
  }else{
    llh <-  x*log(1-p.enter) + y*log(p.enter) + y*log(1-p.stay) + z*log(p.stay)
  }
  llh + sum(llh.list)
}

prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  rate <- param[3]
  
  prior.pe <- dlnorm(p.enter,meanlog=-8,sdlog=0.5, log=T)
  prior.ps <-  dlnorm(p.stay,meanlog=-0.17,sdlog=0.05,log=T)
  prior.rate <- dlnorm(rate,meanlog=log(0.0001), sdlog=0.5,log=T)
  
  return(prior.pe + prior.ps + prior.rate)
}

posterior <- function(param, slip, llh){
  #print(prior(param))
  #print(likelihood(param))
  if (any(param > 1)){
    return(log(0))
  }else{
    return(prior(param) + likelihood(param, slip, llh))
  }
}

proposalFunction <- function(param, slip_current, llh_current){
  p.enter <- param[1]
  p.stay <- param[2]
  rate <- param[3]
  
  num <- runif(1)
  s2p <- 0.90
  
  # CHANGE PARAMETERS 
  if (num > s2p){
    if (num-s2p < (1/3 *(1-s2p))){
      p.enter <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
      llh_proposed <- llh_current   # stays the same
    }else if(num-s2p > (2/3 *(1-s2p))){
      p.stay <- rlnorm(1,meanlog=log(param[2]),sdlog=0.02)
      llh_proposed <- llh_current   # stays the same
    }else{
      rate <- rlnorm(1,meanlog=log(param[3]),sdlog=0.08)
      llh_proposed <- seqllh(rate, slip_current)  # recalcuate using the new rate
    }
    slip_proposed <- slip_current    # stays the same
  
  # CHANGE SLIP
  }else{
    unpack <- changeSlip(slip_current)
    slip_proposed <- unpack[[1]]        # this is the new slip.list with a single change
    changed <- unpack[[2]]              # this is the location at which the slip.list was changed
    
    # recalculate SINGLE pairwise llh at position = changed
    llh_proposed <- llh_current
    new.seq <- getTip(insertions$Vseq[changed], slip.list[[changed]])
    llh_proposed[changed] <- pairllh(anc.seqs[changed], new.seq, rate, branches[changed])
  }
  return(list(param=c(p.enter, p.stay, rate), slip=slip_proposed, llh=llh_proposed))
}

runMCMC <- function(startvalue, iterations, slip.list){
  # timing
  start.time <- proc.time()
  
  # initialize the chain
  chain <- array(dim = c(iterations+1,3))
  chain[1,] <- startvalue
  
  # start the slip list and the llh list
  slip_current <- slip.list
  llh_current <- seqllh(startvalue[3], slip.list)
  
  # keep a logfile up to date
  logfile <- file("~/PycharmProjects/hiv-withinhost/slip-model.csv", "w")
  write("p(Enter), p(Stay), Rate, Slip-changed, Accept", file=logfile)
  
  for (i in 1:iterations){
    # calculate posterior of current position
    p.current <- posterior(chain[i,], slip_current, llh_current)
    
    # generate proposal 
    proposal <- proposalFunction(chain[i,], slip_current, llh_current)
  
    
    # calculate posterior of the new proposal (parameters, slip_proposed, llh_proposed)
    p.next <- posterior(proposal[[1]], proposal[[2]], proposal[[3]])
    
    #print(paste("Current:", p.current, "Next:", p.next, sep=" "))
    #print(paste(proposal[[1]], sep=" "))
    s.change <- any(unname(unlist(slip_current))!=unname(unlist(proposal[[2]])))
    #print(paste0("Sliplist change proposed: ", s.change))
    
    # to catch problematic posterior calculations 
    if(is.na(p.current) || is.na(p.next)){
      print("ERROR: Posterior could not be calculated")
      print(paste("Chain value:", chain[i,1], chain[i,2], chain[i,3], sep=" "))
      break
    }
    
    prop <- exp(p.next - p.current)
    #print(prop)
    
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- proposal[[1]]
      slip_current <- proposal[[2]]
      llh_current <- proposal[[3]]
      #print("Accept")
      accept <- T
    # if the proportion is less than the random uniform sample, REJECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
      #print("Reject")
      accept <- F
    }
    
    if (i %% 10 == 0){
      print(paste("STATE",i,":", chain[i,1], chain[i,2], chain[i,3], sep=" "))
      write(paste(c(chain[i,], as.numeric(s.change), as.numeric(accept), proc.time() - start.time), collapse=",") , file=logfile, append=T)
    }
  }
  return(list(chain=chain, slip=slip_current))
  close(logfile)
}

# RUN MCMC
startvalue <- c(0.001, 0.8, 0.001)
chain <- runMCMC(startvalue, 100000, slip.list)




