# SLIP MODEL 
# CURRENT VERSION : #4-1
# 4-1 :  Addition of fixation parameter
# This parameter will selectively remove non3 indel sequences 
# from analysis using a beta distribution prior
source("~/vindels/2_within-host/slip-model-utils.r")

setup <- function(tip, anc, len, pos, branches, shuffle){
  
  indels <<- data.frame(tip=tip, anc=anc, len=len, pos=pos, stringsAsFactors = F)
  
  branches <<- as.numeric(branches)
  anc.seqs <<- gsub("-", "", anc)
  
  f <<- estimateFreq(c(tip, anc))
  names(f) <<- nt
  
  # list of slip vectors -- one for each sequence -- determining how sequence will be read
  slip.list <<- unname(mapply(createSlips, anc, len, pos))
  
  # ----SHUFFLING --- 
  # randomly shuffle the slip locations around
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
            proposed <- locs[n] + delta(2)
          }
          locs[n] <- proposed
        }
        getSlipVector(locs, length(x))
      }
    })
  }
  
  # locations of all insertion-containing sequences 
  idx <<- which(unname(lapply(slip.list,sum))>0)
  
  # finds all insertion lengths in the slip.list
  lens <- unlist(lapply(slip.list, sum))
  lens <- lens[lens != 0]
  
  # all multiple-of-3 insertion lengths
  a <<- sum(lens %% 3 == 0)
  # all non-3 insertion lengths
  b <<- sum(lens %% 3 != 0)
  
  # REMOVE:  est.total <<- a / 0.29
  # REMOVE: est.fix <<- b / (est.total - a)
  non3 <- 1:100
  non3 <<- non3[-which(non3%%3==0)]
}

likelihood<- function(param, slip.list, llh.list){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  # ---- Fixation Parameter ----
  # this determines the likelihood of the fixation parameter 
  # mean determined by calculated 'est.fix' value
  # param[4] is the proposed 'fix' value
  
  est.non3 <- round((a / (1 - sum(dgeom(non3-1, (1-p.stay))))) - a)
   
  llh.fix <- dbinom(b, est.non3, param[4], log=T)    # mean=est.fix, sd=0.005, log=T)
  
  adj.enter <- p.enter * ((a+b) / (a + est.non3))
  
  # ---- Affine Likelihood ---- 
  # utilizes both adj.enter + p.stay
  slips <- unname(unlist(slip.list))
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  if (any(is.na(slips))){
    aff.llh <- log(0)
  }else{
    aff.llh <-  x*log(1-adj.enter) + y*log(adj.enter) + y*log(1-p.stay) + z*log(p.stay)
  }
  
  aff.llh + sum(llh.list) + llh.fix
}

prior <- function(param){
  prior.pe <- dunif(param[1], min=1e-7, max=1e-2, log=T) 
  prior.ps <- dunif(param[2], min=0.4, max=0.95, log=T) 
  prior.rate <- dunif(param[3], min=1e-7, max=1e-3, log=T)
  prior.fix <- dbeta(param[4], shape1=3, shape2=20, log=T)
  #prior.fixsd <- dexp(param[5], rate=1, log=T)
  #prior.pe <- dlnorm(p.enter,meanlog=log(0.00015),sdlog=1.2, log=T)
  #prior.ps <- dlnorm(p.stay,meanlog=log(0.75),sdlog=0.1,log=T)
  #prior.rate <- dlnorm(rate,meanlog=log(0.00001), sdlog=1.5,log=T)
  
  return(prior.pe + prior.ps + prior.rate + prior.fix) # + prior.fixsd)
}

posterior <- function(param, slip, llh){
  if (any(param > 1)){
    return(log(0))
  }else{
    return(c(prior(param) , likelihood(param, slip, llh)))
  }
}

proposalFunction <- function(param, slip_current, llh_current){

  num <- runif(1)
  s2p <- 0.97
  
  # CHANGE PARAMETERS 
  if (num > s2p){
    num <- (num - s2p) / (1-s2p)
    if (num < 1/5){
      param[1] <- rlnorm(1,meanlog=log(param[1]),sdlog=0.1)
      llh_proposed <- llh_current   # stays the same
    }else if(num > 2/10 && num < 5/10){
      param[2] <- rlnorm(1,meanlog=log(param[2]),sdlog=0.02)
      llh_proposed <- llh_current   # stays the same
    }else if(num > 5/10 && num < 8/10){
      param[3] <- rlnorm(1,meanlog=log(param[3]),sdlog=0.08)
      llh_proposed <- seqllh(param[3], slip_current)  # recalcuate using the new rate
    }else{
      param[4] <- rlnorm(1,meanlog=log(param[4]),sdlog=0.02)
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

runMCMC <- function(startvalue, iterations, runno, notes){
  # timing
  start.time <- proc.time()
  
  # initialize the chain
  chain <- array(dim = c(iterations+1,4))
  chain[1,] <- startvalue
  
  # start the slip list and the llh list
  slip_current <- slip.list
  llh_current <- seqllh(startvalue[3], slip.list)
  
  # keep a logfile up to date
  logfile <- file(paste0("~/PycharmProjects/hiv-withinhost/15_modeling/slip-", 
                         as.character(runno),#substr(gsub("[\\ :-]","",Sys.time()), 9, 12),
                         ".csv"), "w")
  notes <- gsub("^", "#",notes)
  notes <- gsub("\n", "\n#", notes)
  write(notes, file=logfile)
  write("p(Enter), p(Stay), Rate, Fix, Likelihood, Posterior, Slip-changed, Accept, Time", file=logfile, append=T)
  
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
    
    prop <- exp(sum(p.next) - sum(p.current))
    #print(prop)
    
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- proposal[[1]]
      slip_current <- proposal[[2]]
      llh_current <- proposal[[3]]
      llh <- p.next[2]
      #print("Accept")
      accept <- T
      # if the proportion is less than the random uniform sample, REJECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
      #print("Reject")
      accept <- F
      llh <- p.current[2]
    }
    
    if (i %% 100 == 0){
      print(paste(c("STATE",i,":", chain[i,], sum(p.current)), collapse=" "))
      write(paste(c(chain[i,], llh, sum(p.current), as.numeric(s.change), as.numeric(accept), (proc.time() - start.time)[[3]]), collapse=",") , file=logfile, append=T)
    }
    if (i %% 50000 == 0){
      slip_current <<- slip_current
      llh_current <<- llh_current
      slip <- unlist(lapply(slip_current, function(x){
        paste(x, collapse="")
      }))
      writeLines(slip, con=paste0("~/PycharmProjects/hiv-withinhost/15_modeling/list-",as.character(runno),"-",i,'.csv'))
      writeLines(as.character(llh_current), con=paste0("~/PycharmProjects/hiv-withinhost/15_modeling/llh-",as.character(runno),"-",i,'.csv'))
    }
  }
  return(list(chain=chain, slip=slip_current))
  close(logfile)
}
