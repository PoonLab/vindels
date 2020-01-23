require(bbmle)
# used to fill in deletion gaps found in the tip sequences 

source("~/vindels/2_within-host/utils.r")

removeDeletions <- function(vseq, anc){
  if(!grepl("-",vseq)){
    return(vseq)
  }else{
    tip.chars <- strsplit(vseq, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(tip.chars=="-")
    
    for (c in idx){
      tip.chars[c] <- anc.chars[c]
    }
    test <- paste0(tip.chars,collapse="")
    test
  }
}
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

newtip <- function(oldtip, slip){
  nonzeros <- which(slip != 0)
  
  if (length(nonzeros) == 0){
    return(oldtip)
  }else{
    
    newtip <- c()
    tip.chars <- strsplit(oldtip, "")[[1]]
    new.chars <- c()
    n <- 0
    pos <- 0
    while (pos < nchar(oldtip)){
      pos <- pos + 1
      n <- n+1
      #print(pos)
      count <- slip[n]
      #print(count)
      if (count != 0){
        pos <- pos + (count-1)
      }else{
        new.chars <- c(new.chars, tip.chars[pos])
      }
    }
    return(paste0(new.chars,collapse=""))
  }
}

checkDiff <- function(seq1, seq2){
  if (seq1 == seq2){
    return(integer(0))
  }else if (nchar(seq1) != nchar(seq2)){
    return(NA)
  }
  
  seq1 <- str_split(seq1, "")[[1]]
  seq2 <- str_split(seq2, "")[[1]]
  
  chars <- rbind(seq1, seq2)
  which(chars[1,]!=chars[2,])
}

path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]
insertions <- insertions[-c(which(insertions$Pos==0)),]
insertions$Vseq <- unname(mapply(removeDeletions,insertions$Vseq, insertions$Anc))

slip.list <- unname(mapply(createSlips, insertions$Anc, insertions$Seq, insertions$Pos))
insertions$Header <- unname(mapply(patLabel, insertions$Header, insertions$Pat))
names(slip.list) <- insertions$Header

for (elem in 1:length(slip.list)){
  if(sum(slip.list[[elem]]) > 0){
    print(slip.list[elem])
  }
}

getSlipLocations <- function(slips){
  if (sum(slips)== 0){
    return (c(integer(0),length(slips)))
  }else{
    nonzeros <- which(slips!=0)
    locations <- c()
    for (n in nonzeros){
      locations <- c(locations, rep(n, slips[n]))
    }
    return (c(locations,length(slips)))
  }
}

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
ceiling(rnorm(1,mean=0,sd=2))

likelihood <- function(slip){
  if (slip <= 0){
    return(0)
  }else if (slip > 1){
    return(0)
  }
  sum(dgeom(counts,prob=slip, log=T))
}

prior <- function(slip){
  prior <- dunif(slip, log = T)
  
  return(prior)
}


posterior <- function(slip){
  prior(slip) + likelihood(slip)
}

proposalFunction <- function(slip){
  return(rnorm(1,mean=slip, sd=0.01))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,1))
  chain[1,] <- startvalue
  
  
  for (i in 1:iterations){
    
    proposal <- proposalFunction(chain[i,])
    prop <- exp(posterior(proposal) - posterior(chain[i,]))
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
    if (runif(1) < prop) {
      chain[i+1,] <- proposal
      
      # if the proportion is less than the random uniform sample, REJCECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
    }
    if (i %% 100 == 0){
      print(paste0("STATE ",i,": ", chain[i,1]))
    }
  }
  return(chain)
  
}

# RUN MCMC
startvalue <- 0.5
chain <- runMCMC(startvalue, 2000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- 200
acceptance <- 1 - mean(duplicated(chain[-(1:burnin),]))
print(paste0("Acceptance: ", acceptance))
