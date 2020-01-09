
transitionCounts <- function(seq){
  len <- nchar(seq)
  nt <- c("A", "C", "G", "T", "X")
  counts <- matrix(0, nrow=5, ncol=5)
  
  rownames(counts) <- nt
  colnames(counts) <- nt
  #print(seq)
  if (seq != ""){
    for (i in 1:(len-1)){
      x <- substr(seq, i, i)
      y <- substr(seq, i+1 ,i+1)
      counts[x,y] <- counts[x,y] + 1
    }
  }
  counts
}


checkDiff <- function(seq1, seq2){
  if (seq1 == seq2){
    return(NULL)
  }else if (nchar(seq1) != nchar(seq2)){
    return(NA)
  }
  
  seq1 <- str_split(seq1, "")[[1]]
  seq2 <- str_split(seq2, "")[[1]]
  
  chars <- rbind(seq1, seq2)
  which(chars[1,]!=chars[2,])
}


# used to detect similar/identical sequences found adjacent to insertions 
flankCheck <- function(indel,pos,vseq,wobble, offset=0){
  len <- nchar(indel)
  pos <- as.numeric(pos)
  
  beforeBool <- F
  afterBool <- F
  
  beforeIdx <- NaN
  afterIdx <- NaN
  
  beforeDiff <- NaN
  afterDiff <- NaN
  
  beforeSeq <- ""
  afterSeq <- ""
  
  # BEFORE (5')
  # find the best matching position from 0 to offset 
  lowest <- 1000
  bestIdx <- -1
  
  for (idx in 0:offset){
    # subtract length to get the start of the sequence 
    # needs to be enough nucleotides to check
    if ((pos - len - idx) >= 0){
      before <- substr(vseq, pos-len-idx+1, pos-idx)
      #print(before)
      diffs <- checkDiff(indel, before)
      if (length(diffs) < lowest){
        lowest <- length(diffs)
        bestIdx <- idx
      }
    }
  }
  
  if (bestIdx != -1 & lowest <= round(wobble*len)){
    beforeBool <- T
    beforeIdx <- bestIdx
    beforeDiff <- lowest
    beforeSeq <- substr(vseq,pos-len-bestIdx+1, pos-bestIdx)
  }
  
  # AFTER 
  # find the best matching position from 0 to offset 
  lowest <- 1000
  bestIdx <- -1
  for (idx in 0:offset){
    if ((pos + len + idx) <= nchar(vseq)){
      # then the PRECEDING position can be checked
      after <- substr(vseq, pos+len+idx+1, pos+2*len+idx)
      #print(after)
      diffs <- checkDiff(indel, after)
      if (length(diffs) < lowest){
        lowest <- length(diffs)
        bestIdx <- idx
      }
    }
  }
  
  if (bestIdx != -1 & lowest <= round(wobble*len)){
    afterBool <- T
    afterIdx <- bestIdx
    afterDiff <- lowest
    afterSeq <- substr(vseq,pos+bestIdx+1, pos+len+bestIdx)
  }
  
  c(indel, vseq, as.logical(beforeBool),  as.numeric(beforeIdx),  as.numeric(beforeDiff), beforeSeq, as.logical(afterBool), as.numeric(afterIdx), as.numeric(afterDiff),afterSeq)
}


flankProps <- function(indel, pos, vseq){
  len <- nchar(indel)
  pos <- as.numeric(pos)
  
  if ((pos - len - idx) >= 0){
    before <- substr(vseq, pos-len-idx+1, pos-idx)
    #print(before)
    diffs <- checkDiff(indel, before)
    if (length(diffs) < lowest){
      lowest <- length(diffs)
      bestIdx <- idx
    }
  }
  
}

# adds an "X" character to signify the location of an insertion 
addX <- function(seq,pos){
  if (!is.na(pos)){
    paste0(substr(seq,1,pos),"X",substr(seq,pos+1,nchar(seq)))
  }else{
    seq
  }
}


labels <- function(header, patient, vloop){
  letter <- strsplit(patient, "-")[[1]][2]
  paste0(header,"_", letter, "_", vloop)
}


require(stringi)
simulateData <- function(seq){
  
  newseq <- seq
  
  pos <- 1
  
  while (pos <= nchar(seq)) {
    char <- substr(seq, pos,pos)
    p.slip <- sample(20000,1)
    
    #Slip event
    if (p.slip==1) {
      # perform the slip event using slip() and insert it into the sequence
      newseq <- paste0(substr(seq,0,pos), slip(seq,pos,0.09), substr(seq,pos+1,nchar(seq)))
    }
    pos <- pos + 1
    
  }
  newseq
}

slip <- function(seq, pos, p.exit){
  exit <- F 
  result <- ""
  
  count <- 0
  
  
  while (!exit & count < pos){
    result <- paste0(result, substr(seq, pos-count, pos-count))
    count <- count + 1
    # generate a random number, if its less than the probability given, quit 
    rand <- runif(1, min=0,max=1)
    if (rand < p.exit){
      exit <- T
    }
  }
  
  result <- stri_reverse(result)
  result
}