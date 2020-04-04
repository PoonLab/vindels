require(bbmle)
require(ape)
# used to fill in deletion gaps found in the tip sequences 
nt <- c("A", "C", "G", "T")
require(stringr)
source("~/vindels/2_within-host/utils.r")
source("~/vindels/2_within-host/slippage-model.r")


# ------  DATA IMPORT AND MANIPULATION -----

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



nt <- c("A", "C", "G", "T")

# NORMAL DATA: 
f <- estimateFreq(c(insertions$Vseq, insertions$Anc))
branches <- insertions$Date
anc.seqs <- gsub("-", "", insertions$Anc)

# likelihood of the entire slip.list
# only calculated when : 
  # a) rate parameter is changed
  # b) before the MCMC starts for the first iteration

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




