# slippage model functions
# CURRENT VERSION : #4 
# Addition of penalty parameter 
# this parameter will selectively remove non3 indel sequences from analysis based on a beta distribution  
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
  nt <- c("A", "C", "G", "T")
  output <- c()
  for (n in 1:length(nt)){
    counts <- sum(unname(sapply(seqs,function(x) str_count(x, nt[n]))))
    output[n] <- counts / sum(unname(sapply(seqs, nchar)))
  }
  output
}

setup <- function(tip, anc, len, pos, branches, shuffle){
  
  indels <<- data.frame(tip=tip, anc=anc, len=len, pos=pos, stringsAsFactors = F)
  
  branches <<- as.numeric(branches)
  anc.seqs <<- gsub("-", "", anc)
  
  f <<- estimateFreq(c(tip, anc))
  names(f) <<- nt
  
  # generate slip list 
  slip.list <<- unname(mapply(createSlips, anc, len, pos))
  
  # #SHUFFLING --- randomly shuffle the slip locations around
  # slip-model-3  --- location prior added 
  if (shuffle){
    slip.list <<- lapply(slip.list, function(x){
      total <- sum(x)
      if (total == 0){
        return (x)
      }else{
        locs <- getSlipLocations(x)[[1]]
        for (n in 1:length(locs)){
          proposed <- -1
          while(proposed < 0 || proposed > length(x)){
            proposed <- locs[n] + delta(3)
          }
          locs[n] <- proposed
        }
        getSlipVector(locs, length(x))
      }
    })
  }
  # needed for use in the CHANGESLIP function
  idx <<- which(unname(lapply(slip.list,sum))>0)
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
  #print(rate)
  #print(branch)
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

# For use in the proposal function
changeSlip <- function(slip.list){
  slip.idx <- getSlipLocations(slip.list)
  
  # choose a slip event to change
  toEdit <- sample(length(slip.idx[[1]]),1)
  
  # this is to ensure that the proposed change is never outside the slip region
  proposal <- 0
  while(proposal <= 0 || proposal > length(slip.list)){
    proposal <- slip.idx[[1]][toEdit] + delta()
  }
  # save the change to the slip list
  slip.idx[[1]][toEdit] <- proposal
  
  # save the whole list
  getSlipVector(slip.idx[[1]],slip.idx[[2]])
}

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

# likelihood of the entire slip.list
# only calculated when : 
# a) rate parameter is changed
# b) before the MCMC starts for the first iteration
seqllh <- function(rate, slip.list){
  require(parallel)
  #print("Starting SeqLLH ... ")
  tip.seqs <- unname(mcmapply(getTip, indels$tip, slip.list, mc.cores=16))
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
  
  # this calculates the observed probability of seeing non3 indels 
  totals <- unlist(lapply(slip.list, sum))
  totals <- totals[totals != 0]
  obs.fix <- sum(totals %% 3 != 0) / (sum(totals %% 3 == 0) / 0.286)
  
  llh.fix <- dnorm(obs.fix, mean=param[4], sd=param[5])
  
  # affine gap likelihood
  # utilizes both p.enter + p.stay
  slips <- unname(unlist(slip.list))
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  if (any(is.na(slips))){
    slipllh <- log(0)
  }else{
    slipllh <-  x*log(1-p.enter) + y*log(p.enter) + y*log(1-p.stay) + z*log(p.stay)
  }
  slipllh + sum(llh.list) + llh.fix
}

prior <- function(param){
  prior.pe <- dunif(param[1], min=1e-6, max=1e-2, log=T) 
  prior.ps <- dunif(param[2], min=0.4, max=0.9, log=T) 
  prior.rate <- dunif(param[3], min=1e-7, max=1e-3, log=T)
  prior.fix <- dunif(param[4], min=0.01, max=1, log=T)
  prior.fixsd <- dunif(param[5], min=0, max=5, log=T)
  #prior.pe <- dlnorm(p.enter,meanlog=log(0.00015),sdlog=1.2, log=T)
  #prior.ps <- dlnorm(p.stay,meanlog=log(0.75),sdlog=0.1,log=T)
  #prior.rate <- dlnorm(rate,meanlog=log(0.00001), sdlog=1.5,log=T)
  
  return(prior.pe + prior.ps + prior.rate + prior.fix + prior.fixsd)
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

  num <- runif(1)
  s2p <- 0.95
  
  # CHANGE PARAMETERS 
  if (num > s2p){
    num <- (num - s2p) / (1-s2p)
    if (num < 1/5){
      param[1] <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
      llh_proposed <- llh_current   # stays the same
    }else if(num > 2/10 && num < 5/10){
      param[2] <- rlnorm(1,meanlog=log(param[2]),sdlog=0.04)
      llh_proposed <- llh_current   # stays the same
    }else if(num > 5/10 && num < 8/10){
      param[3] <- rlnorm(1,meanlog=log(param[3]),sdlog=0.02)
      llh_proposed <- seqllh(param[3], slip_current)  # recalcuate using the new rate
    }else if(num > 8/10 && num < 9/10){
      param[4] <- rlnorm(1,meanlog=log(param[4]),sdlog=0.05)
      llh_proposed <- llh_current
    }else{
      param[5] <- rlnorm(1,meanlog=log(param[5]),sdlog=0.05)
      llh_proposed <- llh_current
    }
    slip_proposed <- slip_current    # stays the same
    # CHANGE SLIP
  }else{
    # choose a sequence to edit
    rand <- sample(length(idx),1)
    changed <- idx[rand]          # this is the location at which the slip.list was changed
    # convert it to indices
    slip_proposed <- slip_current
    slip_proposed[[changed]] <- changeSlip(slip_current[[changed]]) # this is the new slip.list with a single change
    
    # recalculate SINGLE pairwise llh at position = changed
    llh_proposed <- llh_current
    new.tip <- getTip(indels$tip[changed], slip_proposed[[changed]])
    llh_proposed[changed] <- pairllh(anc.seqs[changed], new.tip, param[3], branches[changed])
  }
  return(list(param=param, slip=slip_proposed, llh=llh_proposed))
}

# MODIFIED TO CONTINUE RUN 

runMCMC <- function(startvalue, iterations, runno){
  # timing
  start.time <- proc.time()
  
  # initialize the chain
  chain <- array(dim = c(iterations+1,5))
  chain[1,] <- startvalue
  
  # start the slip list and the llh list
  slip_current <- slip.list
  llh_current <- seqllh(startvalue[3], slip.list)
  
  # keep a logfile up to date
  logfile <- file(paste0("~/PycharmProjects/hiv-withinhost/slip-model-", 
                         runno,#substr(gsub("[\\ :-]","",Sys.time()), 9, 12),
                         ".csv"), "w")
  write("p(Enter), p(Stay), Rate, Fixation, Fix-SD, Slip-changed, Accept, Time", file=logfile)
  
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
    
    if (i %% 100 == 0){
      print(paste(c("STATE",i,":", chain[i,]), collapse=" "))
      write(paste(c(chain[i,], as.numeric(s.change), as.numeric(accept), (proc.time() - start.time)[[3]]), collapse=",") , file=logfile, append=T)
    }
    if (i %% 100000 == 0){
      slip_current <<- slip_current
      llh_current <<- llh_current
    }
  }
  return(list(chain=chain, slip=slip_current))
  close(logfile)
}
