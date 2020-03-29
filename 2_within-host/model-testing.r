# model testing 
# relies upon modeling2.r


# SIMULATE DNA 
# ------------------
# SIMULATE DNA SEQUENCES 

# this calculates the mean substitution rate across all the observed sequences 
# sum(all substitutions) / sum(all nucl sites)
diffs <- sum(unname(mapply(function(x,y) {length(checkDiff(x,y))}, 
                           tip.seqs, anc.seqs))) / sum(nchar(tip.seqs))


# calculate the median lengths of the variable loops 
ins.v <- split(insertions, insertions$Vloop)
lens <- unname(unlist(lapply(ins.v, function(x){median(x[,"Vlength"])})))
rm (ins.v)
# Generates a sequence that adheres to the nucleotide frequencies of the data 
f <- estimateFreq(allseqs)
names(f) <- nt
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
pairllh <- function(anc, newtip, rate, branch){
  tmat <- getMat(rate,branch)
  achars <- strsplit(anc, "")[[1]]
  tchars <- strsplit(newtip, "")[[1]]
  print(branch)
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

# used to test for the prior probability of P.ENTER
# res <- sapply(1:100, function(x){sum(sapply(1:29250, function(x){sum(runif(120) < 0.000109)}))})
# sum(res!=0)

estimateSubs <- function(tip, anc){
  res <- unname(mapply(function(x,y){length(checkDiff(x,y))}, tip, anc))
  lens <- unname(sapply(tip, nchar))
  sum(res) / sum(lens)
}

erate <- estimateSubs(insertions$Vseq, insertions$Anc)

simPair <- function(p.enter, p.stay, rate){
  vlen <- lens[sample(1:5, 1)]
  anc <- genSeq(vlen)
  
  # make a copy of the ancestor to edit
  tip <- anc
  
  # SUBSTITUTIONS
  # add in substitutions based on the rate value 
  subs <- which(sapply(1:vlen, function(x){runif(1) < rate}))
  
  if (length(subs) > 0){
    tip.chars <- str_split(tip, "")[[1]]
    anc.chars <- str_split(anc, "")[[1]]
    
    for (i in subs){
      # choose whether to mutate the tip sequence or ancestor (assumed 50% chance)
      if (runif(1) < 0.5){
        toChange <- tip.chars
        tipBool  <- T
      }else{
        toChange <- anc.chars
        tipBool <- F
      }
      
      # store the nucleotide at the chosen location
      n.idx <- which(names(f)==toChange[i])
      # calculate the probabilities of changing to each of the others nucleotides 
      t.probs <- f[-n.idx] / (1- f[[n.idx]])
      
      # based on the probabilities, change to the new nucleotide
      rnum <- runif(1)
      if (rnum < t.probs[[1]]){
        out.nt <-  names(t.probs)[1]
      }else if(rnum > t.probs[[1]] && rnum < (t.probs[[1]]+t.probs[[2]])){
        out.nt <- names(t.probs)[2]
      }else{
        out.nt <- names(t.probs)[3]
      }
      
      if (tipBool){
        tip.chars[i] <- out.nt
      }else{
        anc.chars[i] <- out.nt
      }
      
    }
    tip <- paste0(tip.chars,collapse="")
    anc <- paste0(anc.chars,collapse="")
  }
  
  # INDELS 
  # determine the number of insertionsd that occur 
  count <- sum(runif(vlen) < p.enter)
  if (count > 0){
    for (n in 1:count){
      if (runif(1) < 0.18){
        # add a non-multiple of 3 indel
        
        len <- 0
        # used to count how many SUBSEQUENT nucleotides will be added AFTER the first one
        while (len %% 3 == 0){
          # algorithm to compute an appropriate length
          len <- 0
          exit <- 1
          while(exit > p.stay){
            len <- len + 1
            exit <- runif(1)
          }
        }
        
      }else{
        # add a multiple of 3 indel 
        len <- 1
        # keep trying until an indel of length > 0, and multiple of 3 is found 
        while (len %% 3 != 0 || len == 0){
          len <- 0
          exit <- 1
          # grows an indel size
          while(exit > p.stay){
            len <- len + 1 
            exit <- runif(1)
          }
        }
      }
      # add the length and choose a random location for the slip event 
      idx <- sample(1:length(vlen)+1,1)
      
      # generate the sequence to insert / delete 
      indel <- genSeq(len)
      
      # generate the tip and ancestor sequence by adding / removing sequence 
      tip <- insert(tip, indel, idx)
      anc <- insert(anc, rep("-",len), idx)
    }
  }
  
  
  
  
  return(list(tip=tip,anc=anc))
  # to estimate the ACTUAL rate of indels, you need to make failures occur after p.enter has been selected
  # algorithm:
  # probability of enter is chosen
  # draw a number from a poisson process
  # if the number is %% 3 == 0: 
  # keep it 100 percent
  # else if the number if %%3 != 0:
  # there's a low probability that it will be kept (penalty)
  
  # to estimate the OBSERVED rate of indels, it is a simpler process
  # algorithm:
  # probability of enter is chosen
  # draw a random number to determine whether a 3 or non3 should be added
  # if runif(1) < 0.05:
  # add a non3 number 
  # else:
  # add a multiple of three 
}

# -----------------