require(Biostrings)
# for changing headers from ACCNO_DATE format to ACCNO
getSubtype <- function(header){
  newheader <- strsplit(as.character(header),"\\.")[[1]][1]
  newheader
}

patLabel <- function(header, pat){
  label <- strsplit(pat, "-")[[1]][2]
  paste0(header,"_",label)
}
getLoop <- function(header, vloop){
  paste0(header,"_",vloop)
}

# used for handling entire columns of NA values
removeNA <- function(input, repl=""){
  if (all(is.na(input))){
    input <- repl
  }
  input
}

# retrieves the accno field from the full scale header 
getAccno <- function(input){
  accno <- strsplit(input, "\\.")[[1]][5]
  accno
}

bootstrap <- function(vector, replicates){
  
}


# specifically handles fields containing a comma
# and copies their data so they can be split into individual rows in a data frame
splitRows <- function(row){
  row <- data.frame(t(row),stringsAsFactors = F)
  seqs <- str_split(row[1,6], ",")[[1]]
  pos <- str_split(row[1,7],",")[[1]]
  len <- length(seqs)
  data.frame(row[rep(1,len),1:5], Seq=seqs, Pos=pos, row[rep(1,len),8:10])
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


# uses the start position of the given vloop (retrieved from var.pos; separate file)
# to compute the indel position within gp120
addPos <- function(pos, header, vloop){
  if (pos == ""){
    return(NA)
  }
  var.pos <- read.csv(paste0(path,"3RegionSequences/variable/", strsplit(filename, "-")[[1]][1], ".csv"), stringsAsFactors = F)
  var.pos <- var.pos[,-c(2,5,8,11,14)]
  
  pos <- as.numeric(pos)
  vloop <- as.character(vloop)
  arg2 <- as.numeric(var.pos[var.pos$header==gsub("_\\d*$","",header), paste0('start.',vloop)])
  pos + arg2
}

transitionCounts <- function(seq){
  len <- nchar(seq)
  nt <- c("A", "C", "G", "T")
  counts <- matrix(0, nrow=5, ncol=5,dimnames=list(nt,nt))
  
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

# this function uses the START positions of insertions, not the end positions
insert <- function(str, indel,pos ){
  #pos <- pos+1
  vect <- strsplit(str,"")[[1]]
  if (pos < 1){
    return (NA)
  }
  
  if (pos == 1){
    return (paste(c(indel, vect[1:length(vect)]),collapse=""))

    # adds insert after the vector (at the end )
  }else if (pos == length(vect)+1){
    return (paste(c(vect[1:length(vect)], indel),collapse=""))
    # involves slicing the before and after parts of the vector around the insert
    # used when pos = 2 : nchar-1
  }else{
    return (paste(c(vect[1:pos-1], indel, vect[pos:length(vect)]),collapse=""))
  }
}
delete <- function(str, indel, pos){
  vect <- rep(T, nchar(str))
  
  if (pos > (nchar(str) - nchar(indel) + 1)){
    return (NA)
  }
  
  start <- pos
  end <- pos + nchar(indel) - 1
  
  vect[start:end] <- F
  
  chars <- strsplit(str, "")[[1]]
  
  return(paste(chars[which(vect)],collapse=""))
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

# 14_flanking 
flankCheck <- function(indel,pos,vseq,wobble=1/6, offset=0){
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
      before <- substr(vseq, pos-2*len-idx+1, pos-idx-len)
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
      after <- substr(vseq, pos+idx+1, pos+len+idx)
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

# 14_flanking
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

# modeling 
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

# 13_1 nglycs and modeling
restoreDel <- function(tip, anc, indel, pos){
  # this is used to restore any gaps in tip sequences containing insertions 
  # Reasoning:
  # I need to restore the tip sequence to its ORIGINAL STATE where no deletions have occurred
  # want to focus on investigating specific insertions at one time 
  # Any and all gaps in the tip sequence are deletions and need to be restored
  
  if(!grepl("-",tip) || indel == "" || is.na(pos)){
    return(c(tip,pos))
  }else{
    tip.chars <- strsplit(tip, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(tip.chars=="-")
    
    # perform a readjustment of the position (for insertions only 
    if (any(idx < pos)){
      pos <- as.character(pos + sum(idx < pos))
    }
    
    tip.chars[idx] <- anc.chars[idx]
    tip <- paste0(tip.chars,collapse="")
    return(c(tip, pos))
  }
}

# 13_1 nglycs and modeling
restoreIns <- function(tip, anc, indel){
  # this is used to restore any gaps in ancestor sequences containing deletions 
  # Reasoning:
  # I need to restore the ancestor sequence to its ORIGINAL STATE where no insertions have occurred
  # want to focus on investigating specific insertions at one time
  # Any and all gaps in the tip sequence are deletions and need to be restored
  
  if(!grepl("-",anc)){
    return(anc)
  }else{
    tip.chars <- strsplit(tip, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(anc.chars=="-")
    
    anc.chars[idx] <- tip.chars[idx]
    anc <- paste0(anc.chars,collapse="")
    return(anc)
  }
}




# general
csvcount <- function(input,delim=","){
  if (is.na(input)){
    return(0)
  }
  commas <- str_count(input, delim)
  if (commas > 0){
    result <- commas + 1  
  }else if(is.nan(input) || input == ""){
    result <- 0
  }else{
    result <- 1
  }
  result
}

# general
# used for extracting condensed CSV information 
extractInfo <- function(input){
  if (length(input)==1 && input == ""){
    return(c("",""))
  }else{
    insertions <- strsplit(input, ":")
  }
  seq <- c()
  pos <- c()
  
  for (ins in insertions[[1]]){
    fields <- strsplit(ins, "-")
    seq <- c(seq, fields[[1]][1])
    pos <- c(pos, as.numeric(fields[[1]][2]) )
  }
  return(c(paste(seq,collapse=","), paste(pos,collapse=",")))
}


# OLD ----------------
# used to compute ancestor sequences using the tip sequence  
insOriginal <- function(indel, pos, vseq){
  if (indel == ""){
    return(vseq)
  }
  pos <- as.character(pos)
  seqs <- strsplit(indel, ",")[[1]]
  idxs <- as.numeric(strsplit(pos, ",")[[1]])
  # iterate through the sequences and positions
  for (i in 1:length(seqs)){
    len <- nchar(seqs[i])
    idx <- idxs[i]
    
    # cut out insertion : substring before and up to start of insertion, substring from end of insertion until the end 
    vseq <- paste0(substr(vseq, 0, idx-len) , substr(vseq, idx+1, nchar(vseq)))
    idxs[(i+1):length(idxs)] <- idxs[(i+1):length(idxs)] - len
  }
  vseq
}

delOriginal <- function(indel, pos, vseq){
  if (indel == ""){
    return(vseq)
  }
  
  seqs <- strsplit(indel, ",")[[1]]
  idxs <- as.numeric(strsplit(pos, ",")[[1]])
  # iterate through the sequences and positions
  for (i in 1:length(seqs)){
    len <- nchar(seqs[i])
    idx <- idxs[i]
    
    # add back deletion : substring before and up to the point of deletion, deletion sequence, substring after point of deletion until end 
    vseq <- paste0(substr(vseq, 0, idx) , seqs[i], substr(vseq, idx+1, nchar(vseq)))
    idxs[(i+1):length(idxs)] <- idxs[(i+1):length(idxs)] + len
  }
  vseq
}


# 13 -- glycosylation analysis 
translate <- function(dna) {
  require(ape)
  
  if (nchar(dna) %% 3 != 0) {
    return(NA)
  }
  dnabin <- as.DNAbin(DNAString(dna))
  aabin <- trans(dnabin)[[1]]
  aaseq <- paste(as.character(aabin),collapse="")
  aaseq
}

# takes in an amino acid sequence and returns the locations of all Nglycs
extractGlycs <- function(aaseq){
  result <- gregexpr("N[^P][ST][^P]", aaseq)[[1]]  # used for 0 indexing these position values for analysis in python 
  paste(result, collapse=",")
}

countGlycs <- function(field){
  if ("," %in% field){
    return(str_count(field, ",") + 1)
  }else if (field == -1){
    return(0)
  }else{
    return(1)
  }
}
