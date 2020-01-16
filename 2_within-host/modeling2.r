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

insert <- function(vect, pos, ins){
  if (pos == 1){
    return (c(ins, vect[2:length(vect)]))
    
    # adds insert after the vector (at the end )
  }else if (pos == length(vect)){
    return (c(vect[1:length(vect)-1], ins))
    # involves slicing the before and after parts of the vector around the insert
    # used when pos = 2 : nchar-1 
  }else{
    return (c(vect[1:pos-1], ins, vect[pos:length(vect)]))
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


insertions <- read.csv("~/PycharmProjects/hiv-withinhost/10_nucleotide/ins-sep-all.csv",row.names=1, stringsAsFactors = F)
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]
insertions$Vseq <- unname(mapply(removeDeletions,insertions$Vseq, insertions$Anc))

slip.list <- unname(mapply(createSlips, insertions$Anc, insertions$Seq, insertions$Pos))
insertions$Header <- unname(mapply(patLabel, insertions$Header, insertions$Pat))
names(slip.list) <- insertions$Header

for (elem in 1:length(slip.list)){
  if(sum(slip.list[[elem]]) > 0){
    print(slip.list[elem])
  }
}

insertions <- insertions[-c(which(insertions$Pos==0)),]

# #slips <- matrix(rep(rep(0,121),10), nrow=10,ncol=121)
# slips <- c(rep(0,121000))
# slips[sample(121000,100)] <- sample(3,100,replace = T) * 3
slips <- nchar(insertions[insertions$Count!=0, "Seq"])
slips <- slips[-167]
slips <- c(slips, rep(0, sum(nchar(insertions$Vseq))))

# randomly shuffle all the entries
rnd <- sample(length(slips),length(slips))
slips <- slips[rnd]





# EXPERIMENTAL -------------------------------
# slip = data frame of slip events for EACH SEQUENCE 

# GOAL: 
# make a likelihood function that will compute the probability of a tip/anc pair 
# it will do the following: 
# a) count the number of substitution differences 

llh.pair <- function(tip, anc, slips){
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  
}

# attempted to find the log - likelihood of an affine model 
likelihood2 <- function(slip, param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  x <- sum(slip == 0)
  y <- sum(slip != 0)
  z <- sum(slip[which(slip!=0)] - 1)
  
  llh <- c()
  # Log likelihood of each tip/anc pair
  #(1 - p.slip)^x * p.slip^y * (1-p.stay)^y * p.stay^z
  llh <- 2 * y * x * z * log(1-p.enter) * log(p.enter) * log(1-p.stay) * log(p.stay)
  sum(llh)
  
}
# --------------------------------------------------

likelihood <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  x <- sum(slips == 0)
  y <- sum(slips != 0)
  z <- sum(slips[which(slips!=0)] - 1)
  
  (1-p.enter)^x*(p.enter)^y*(1-p.stay)^y*(p.stay)^z
}



prior <- function(param){
  p.enter <- param[1]
  p.stay  <- param[2]
  
  prior.pe <- dlnorm(p.enter,meanlog=-40,sdlog=5)
  prior.ps <- dlnorm(p.stay,meanlog=-0.3,sdlog=0.15)
  
  return(prior.pe + prior.ps)
}


posterior <- function(param){
  print(prior(param))
  print(likelihood(param))
  prior(param) + likelihood(param)
}

proposalFunction <- function(param){
  p.enter <- rlnorm(1,meanlog=param[1],sdlog=0.1)
  p.stay <- rlnorm(1,meanlog=param[2],sdlog=0.01)
  return(c(p.enter,p.stay))
}


runMCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,2))
  chain[1,] <- startvalue
  
  
  for (i in 1:iterations){
    
    proposal <- proposalFunction(chain[i,])
    print(proposal)
    print(posterior(proposal))
    print(posterior(chain[i,]))
    prop <- exp(posterior(proposal) - posterior(chain[i,]))
    print(prop)
    # if the proportion exceeds the random uniform sample, ACCEPT the proposed value
  *  if (runif(1) < prop) {
      chain[i+1,] <- proposal
      
      # if the proportion is less than the random uniform sample, REJCECT the proposed value stick with current 
    } else {
      chain[i+1,] <- chain[i,]
    }
    if (i %% 100 == 0){
      print(paste0("STATE ",i,": ", chain[i,]))
    }
  }
  return(chain)
  
}

# RUN MCMC
startvalue <- c(0.0000000000001,0.6)
chain2 <- runMCMC(startvalue, 5000)


# sets the burnin size, removes all rows from the chain that are associated with the burnin 
burnin <- 10000
acceptance <- 1 - mean(duplicated(chain2[-(1:burnin),]))

# PLOTTING 
hist(chain[-(1:burnin),1],nclass=30, main="Posterior of x", xlab="True value = red line" )
abline(v = median(chain[-(1:burnin),1]))

plot(chain[-(1:burnin),1], type = "l", xlab="MCMC Steps" , main = "Chain values of x")
abline(h = median(chain[-(1:burnin),1]), col="red")