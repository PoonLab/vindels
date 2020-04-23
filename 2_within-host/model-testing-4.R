# model testing 
# relies upon modeling2.r
# used to test the 

# SIMULATE DNA 
# ------------------
# SIMULATE DNA SEQUENCES 
source('~/vindels/2_within-host/utils.r')
source("~/vindels/2_within-host/slippage-model-4-1.r")
# calculate the median lengths of the variable loops
# ins.v <- split(insertions, insertions$Vloop)
# lens <- unname(unlist(lapply(ins.v, function(x){median(x[,"Vlength"])})))
# f <- estimateFreq(c(insertions$Vseq, insertions$Anc))
# names(f) <- nt
# rm (ins.v)
# rm(insertions)
lens <- c(75, 126, 105,  90 , 33)
f <- c(0.4158261, 0.1649874, 0.1761463, 0.2423612 )
names(f) <- nt

# Generates a sequence that adheres to the nucleotide frequencies of the data 
genSeq <- function(len){
  seq <- c()
  for (n in 1:len){
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
    seq[n] <- nucl
  }
  paste(seq, collapse="")
}

nt <- c("A", "C", "G", "T")

# used to test for the prior probability of P.ENTER
# res <- sapply(1:100, function(x){sum(sapply(1:29250, function(x){sum(runif(120) < 0.000109)}))})
# sum(res!=0)
#lens <- c(10,20,20,30,30)
simPair <- function(p.enter, p.stay, rate){
  vlen <- lens[sample(1:5, 1)]
  anc <- genSeq(vlen)
  
  # pick a branch rate value from a distribution 
  branch <- rlnorm(1, meanlog=4, sdlog=1.05)
  
  # ----SUBSTITUTIONS--------
  # add in substitutions based on the rate value 
  
  anc.chars <- strsplit(anc, "")[[1]]
  tmat <- getMat(rate, branch)
  
  # generate the tip sequence based on the transition rate matrix 
  tip <- paste(sapply(1:nchar(anc), function(n){
    rand <- runif(1)
    probs <- tmat[anc.chars[n],]
    
    if (rand < probs[1]){
      "A"
    }else if (rand > probs[1] && rand < probs[2]){
      "C"
    }else if (rand > (probs[1]+probs[2]) && rand < probs[3]){
      "G"
    }else{
      "T"
    }
  }), collapse="")
  
  # ------ INDELS -------
  # determine the number of insertions that occur 
  count <- sum(runif(vlen) < p.enter)
  
  if (count > 0){
    
    for (n in 1:count){
      len <- 0
      exit <- 0
      while(exit < p.stay ){
        len <- len + 1
        exit <- runif(1)
      }
      
      # add the length and choose a random location for the slip event 
      idx <- sample(1:vlen,1)
      
      # generate the sequence to insert / delete 
      indel <- genSeq(len)
      
      # generate the tip and ancestor sequence by adding / removing sequence 
      tip <- insert(tip, indel, idx)
      anc <- insert(anc, rep("-",len), idx)
    }
  }
  #print("finished indels")
  return(list(tip=tip,anc=anc,branch=branch))
  # to estimate the ACTUAL rate of indels, you need to make failures occur after p.enter has been selected
  # algorithm:
  # probability of enter is chosen
  # draw a number from a poisson process
  # if the number is %% 3 == 0: 
  # keep it 100 percent
  # else if the number if %%3 != 0:
  # there's a low probability that it will be kept (penalty)
  
}

# vec <- c()
# for (i in 1:50){
# SIMULATE TIP + ANCESTOR SEQUENCES
all.seqs <- sapply(1:25000, function(n){
  #print(n)
  pair <- simPair(0.00016, 0.75, 0.00001)
  # VALUE 1 = Tip, VALUE 2 = Ancestor, VALUE 3 = Branch length
  return(c(pair[[1]], pair[[2]], pair[[3]]))
})



insertions <- as.data.frame(t(all.seqs), stringsAsFactors = F)
colnames(insertions) <- c("tip", "anc", "branch")

# Generate the length and position columns
data <- t(unname(sapply(insertions$anc, function(x){
  # returns c(length, position) of insertion events
  gaps <- gregexpr("-",x)[[1]]
  if (length(gaps) == 1 && gaps == -1){
    return (c(NA, NA))
  }else{
    return (c(length(gaps), max(gaps)))
  }
})))

insertions$len <- data[,1]
insertions$pos <- data[,2]

# vec <- c()
# for (i in 1:1000){
#   x <- round(rgeom(1000, 0.17))
#   x <- x[x>0]
#   vec[i] <- sum(x%%3==0) / 1000
# }
# median(vec)   # this number is used to calibrate the MCMC
# rm(vec)

# VERSION 4:
# generate a vector of T/F values to remove non3 values from the data
toInclude <- sapply(insertions$len, function(len){
  if (is.na(len)){
    T
  }else if (len %% 3 != 0){
    runif(1) < 0.12
  }else{
    T
  }
})

# filter out insertions based on the fixation parameter
new.ins <- insertions[toInclude,]
rm(insertions)
setup(new.ins$tip, new.ins$anc, new.ins$len, new.ins$pos, new.ins$branch, F)

# ----- algorithm testing ---

# all.len <- insertions[!is.na(insertions$len), 'len']
# filt.len <- new.ins[!is.na(new.ins$len), 'len']
# sum(all.len %% 3 == 0)
# sum(all.len %% 3 != 0)
# x <- sum(filt.len %% 3 == 0)
# vec[i] <- sum(filt.len %% 3 != 0) / ((x / 0.276) - x)
#}


# RUN MCMC
startvalue <- c(0.01, 0.55, 0.000001, 0.3)
chain <- runMCMC(startvalue, 200000, '10', 'fix3')


# fix3 : (0.00016, 0.75, 0.00001)
# loc-prior2 : (0.00016, 0.75, 0.00001)

# --- print out the whole slip list  ----
indels$slip <- lapply(lapply(slip_current, function(x){
  getSlipLocations(x)[[1]]}), function(slip){
    if(is.null(slip)){
      NA
    }else{
      paste(slip, collapse=',')
    }
})

input <- readLines("~/PycharmProjects/hiv-withinhost/15_modeling/list-10.csv")
input <- unname(sapply(input, function(x){
  unname(sapply(strsplit(x, "")[[1]], as.numeric))
}))

# ----- For checking -----
csv <- read.csv("~/PycharmProjects/hiv-withinhost/slip-model-v4-1.csv", stringsAsFactors = F, skip=1, header=F)
colnames(csv) <- c('p.enter', 'p.stay', "rate" ,'slip.changed', 'accept', 'time')

# -----IDEA TO IMPLEMENT ------
# to estimate the ACTUAL rate of indels, you need to make failures occur after p.enter has been selected
# algorithm:
# probability of enter is chosen
# draw a number from a poisson process
# if the number is %% 3 == 0: 
# keep it 100 percent
# else if the number if %%3 != 0:
# there's a low probability that it will be kept (penalty)
