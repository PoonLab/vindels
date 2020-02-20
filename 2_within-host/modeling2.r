require(bbmle)
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


# used to go from c(0,0,0,3,0,0,0) to c(4,4,4)
getSlipLocations <- function(slip){
  nonzeros <- which(slip!=0)
  locations <- c()
  for (pos in nonzeros){
    locations <- c(locations, rep(pos, slip[pos]))
  }
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


changeSlip <- function(slip.list=slip.list){
  idx <- which(unname(lapply(slip.list,sum))>0)
  seq <- sample(length(idx),1)
  slip <- slip.list[[idx[seq]]]
  slip.idx <- getSlipLocations(slip)
  
  # choose a slip event to change
  toEdit <- position(length(slip.idx[[1]]))
  slip.idx[[1]][toEdit] <- slip.idx[[1]][toEdit] + delta()
  slip.list[[idx[seq]]] <<- getSlipVector(slip.idx[[1]],slip.idx[[2]])
}

path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

# PROBLEMATIC CASE: remove instances with gaps in the ancestor but NO INSERTION
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]

# CASE: remove instances with insertion position 0
insertions <- insertions[-c(which(insertions$Pos==0)),]

# Restore all deletions found in tip sequences and adjust the POS values accordinaly
res <- as.data.frame(t(unname(mapply(restoreDel,insertions$Vseq, insertions$Anc, insertions$Seq, insertions$Pos))))
insertions$Vseq <- as.character(res[,1])
insertions$Pos <- as.numeric(as.character(res[,2]))

# generate slip list 
slip.list <- unname(mapply(createSlips, insertions$Anc, insertions$Seq, insertions$Pos))

# add the a/b replicate label to the headers
insertions$Header <- unname(mapply(patLabel, insertions$Header, insertions$Pat))
names(slip.list) <- insertions$Header

# for printing the slip.list 
for (elem in 1:length(slip.list)){
  if(sum(slip.list[[elem]]) > 0){
    print(slip.list[elem])
  }
}
# C.-.-.QT.10R.-.-_289_1_b
# randomly shuffle the slip locations around 
slip.list <- lapply(slip.list, function(x){
  total <- sum(x)
  locs <- sample(length(x), total, replace=T)
  getSlipVector(locs, length(x))
})

nucleotides <- c("A","C","G","T")


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


# transitionCounts <- function(seq1, seq2){
#   len <- nchar(seq1)
#   nt <- c("A", "C", "G", "T")
#   counts <- matrix(0, nrow=4, ncol=4,dimnames=list(nt,nt))
#   
#   #chars1 <- strsplit(seq1, "")[[1]]
#   #chars2 <- strsplit(seq2, "")[[1]]
#   if (seq1 != ""){
#     for (i in 1:(len-1)){
#       x <- substr(seq1, i, i)
#       y <- substr(seq2, i ,i)
#       counts[x,y] <- counts[x,y] + 1
#     }
#   }
#   counts
# }
nt <- c("A", "C", "G", "T")
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

getMat <- function(rate){
  nt <- c("A", "C", "G", "T")
  mat <- matrix(rep(f, each=4), nrow=4, ncol=4,dimnames=list(nt,nt))
  mat <- mat * rate
  diag(mat) <- sapply(1:4, function(x) -(sum(rate * f[-x]))) 
  mat
}


getTMat <- function(mat, branch){
  mat <- branch * mat
  tmat <- expm(mat)
  tmat
}

allseqs <- c(insertions$Vseq, insertions$Anc)
f <- estimateFreq(allseqs)
pairllh <- function(seq1, seq2, rate, branch){
  
  tmat <- getTMat(getMat(rate), branch)
  
}



likelihood <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  # Log likelihood of each tip/anc pair
  #(1 - p.slip)^x * p.slip^y * (1-p.stay)^y * p.stay^z
  llh <-  x*log(1-p.enter) + y*log(p.enter) + y*log(1-p.stay) + z*log(p.stay)
  llh
  
}
# --------------------------------------------------

prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  prior.pe <- dlnorm(p.enter,meanlog=-10,sdlog=2, log=T)
  prior.ps <- dlnorm(p.stay,meanlog=-0.15,sdlog=0.05,log=T)
  
  return(prior.pe + prior.ps)
}

posterior <- function(param){
  #print(prior(param))
  #print(likelihood(param))
  prior(param) + likelihood(param)
}

proposalFunction <- function(param){
  p.enter <- param[1]
  p.stay <- param[2]
  
  num <- runif(1)
  if (num > 0.9){
    p.enter <- rlnorm(1,meanlog=param[1],sdlog=0.1)
    p.stay <- rlnorm(1,meanlog=param[2],sdlog=0.01)
  }else{
    # perform a change on the sliplist 
    changeSlip()
    # this will NOT change your two parameters, but WILL change the likelihood
  }
  
  return(c(p.enter,p.stay))
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
