require(stringr)
require(ape)
findAncestor <- function(header, tree){
  # this expression will return results for NODES ONLY
  # second column provides the CAPTURED TIP LABELS from within the node label
  header <- str_replace(header, "_[0-9]$", "")
  tips <- str_match_all(header,"([^\\)\\(,\n:]+):")[[1]][,2]
  if (length(tips) == 0){
    # no colons; this means its a TIP 
    # the index in the tre$tip.label vector is the final result
    index <- match(header, tree$tip.label)
  }else{
    # retreive all descendants of every node and tip in the tree
    desc <- Descendants(tree)
    
    # find the numeric labels of all extracted tips 
    matches <- match(tips, tree$tip.label)
    
    # find the SINGLE node in the descendants list that contains the exact same subset of tips
    index <- which(sapply(desc, function(x){ifelse(length(x) == length(matches) && all(x==matches),T,F)}))
  }
  if (length(index)!=1){
    return(paste0("PROBLEM:",as.character(index)))
  }
  return(index)
}

bootstrap <- function(sampl, n=1000){
  replicate(n, mean(sample(sampl, replace=T)))
}

getSubtype <- function(header){
  newheader <- strsplit(as.character(header),"\\.")[[1]][1]
  newheader
}
consecutive <- function(vect){
  cons <- 0
  vals <- c()
  for (i in 2:length(vect)){
    if (vect[i] == (vect[i-1]+1)){
      cons = cons + 1
    }else{
      vals <- c(vals, cons)
      cons <- 0
    }
  }
}


getPat <- function(header, pat){
  label <- strsplit(pat, "-")[[1]][2]
  paste0(header,"_",label)
}
extractPat <- function(header){
  x <- strsplit(header,"_")[[1]]
  c(x[1], x[2])
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


# specifically handles fields containing a comma
# and copies their data so they can be split into individual rows in a data frame
splitRows <- function(row,col.idx){
  row <- data.frame(t(row),stringsAsFactors = F)
  seqs <- str_split(row[1,col.idx[1]], ",")[[1]]
  pos <- str_split(row[1,col.idx[2]],",")[[1]]
  len <- length(seqs)
  data.frame(row[rep(1,len),1:(min(col.idx)-1)], 
             indel=seqs, 
             pos=pos, 
             row[rep(1,len),(max(col.idx)+1):ncol(row)])
}


# adds an "X" character to signify the location of an insertion 
addX <- function(seq,pos){
  if (!is.na(pos)){
    paste0(substr(seq,1,pos),"X",substr(seq,pos+1,nchar(seq)))
  }else{
    seq
  }
}


labels <- function(header, patient){
  letter <- str_match(patient , "-([ab])")[1,2]
  paste(header,letter,sep="_")
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
  
  # Insert sequence immediately before the "pos"th position
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
  
  if (pos > (nchar(str))){
    print("DELETE PROBLEM")
    return (NA)
  }

  start <- pos - nchar(indel) + 1
  end <- pos
  #if(class(start) != "numeric") print("DELETE PROBLEM");print(class(start))
  if (start < 1 ){
    print("DELETE -- INCORRECT START")
  }
  vect[start:end] <- F
  
  chars <- strsplit(str, "")[[1]]
  
  return(paste(chars[which(vect)],collapse=""))
}

translate <- function(dna) {
  if (!is.na(dna)) {
    if (nchar(dna) %% 3 != 0) {
      extra <- nchar(dna) %% 3
      dna <- substr(dna,1, nchar(dna)-extra)
    }
    
    dnabin <- as.DNAbin(DNAString(dna))
    aabin <- trans(dnabin)[[1]]
    aaseq <- paste(as.character(aabin),collapse="")
    aaseq <- gsub("X", "-", aaseq)
    return(truncate(aaseq))
  }else{
    print("START")
    print(dna)
    print(class(dna))
    return("")
  }
}



checkDiff <- function(seq1, seq2){
  require(stringr)
  if (seq1 == seq2){
    return(integer(0))
  }else if (nchar(seq1) != nchar(seq2)){
    return(NA)
  }
  
  seq1 <- str_split(seq1, "")[[1]]
  seq2 <- str_split(seq2, "")[[1]]
  
  which(seq1 != seq2)
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
restoreOtherIndels <- function(anc, tip, indel, pos){
  # find the location of all gap characters in tip/anc
  gaps <- gregexpr("-",anc)[[1]]
  if (length(gaps) > 1 || gaps != -1){
    # create a vector of positions to be copied over
    idx <- rep(F, nchar(tip))
    idx[gaps] <- T
    
    # retrieve the boundaries of the indel
    end <- as.numeric(pos)
    start <- as.numeric(pos) - nchar(indel) + 1
    
    # make sure it ignores the indel sequence
    toIgnore <- start:end
    idx[toIgnore] <- F
    
    # fill in all other insertions / deletions that are not the one of interest ***
    anc.chars <- strsplit(anc, "")[[1]]
    tip.chars <- strsplit(tip, "")[[1]]
    anc.chars[idx] <- tip.chars[idx]
    anc <- paste0(anc.chars, collapse="")
    return(anc)
    
  }
  return(anc)
}

restoreInsAnc2 <- function(anc, indel, pos){
  # find the location of all gap characters in tip/anc
  
  # create a vector of positions to be copied over
  idx <- rep(T, nchar(anc))
  
  # retrieve the boundaries of the indel
  end <- as.numeric(pos)
  start <- as.numeric(pos) - nchar(indel) + 1
  
  # make sure it ignores the indel sequence
  toIgnore <- start:end
  idx[toIgnore] <- F
  
  # fill in all other insertions / deletions that are not the one of interest ***
  anc.chars <- strsplit(anc, "")[[1]]
  anc <- paste0(anc.chars[idx], collapse="")
  return(anc)
  
}

restoreAllGaps <- function(tip, anc){
  # this is used to restore any deletion gaps in tip sequences containing insertions 
  # Reasoning:
  # I need to restore the tip sequence to its ORIGINAL STATE where no deletions have occurred
  # want to focus on investigating specific insertions at one time 
  # Any and all gaps in the tip sequence are deletions and need to be restored
  
  if(!grepl("-",tip)){
    return(tip)
  }else{
    tip.chars <- strsplit(tip, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(tip.chars=="-")
    
    # perform a readjustment of the position (for insertions only 
    # if (!is.na(pos)){
    #   pos <- gregexpr("[ACTG]", tip)[[1]][pos]
    # }
    tip.chars[idx] <- anc.chars[idx]
    tip <- paste0(tip.chars,collapse="")
    return(tip)
  }
  # this is used to restore any insertion gaps in ancestor sequences containing deletions
  # Reasoning:
  # I need to restore the ancestor sequence to its ORIGINAL STATE where no insertions have occurred
  # want to focus on investigating specific insertions at one time
  # Any and all gaps in the tip sequence are deletions and need to be restored
  
}

# for 14_flanking ONLY 
restoreTipDel <- function(tip, anc, indel, pos){
  # this is used to restore any deletion gaps in tip sequences containing insertions 
  # Reasoning:
  # I need to restore the tip sequence to its ORIGINAL STATE where no deletions have occurred
  # want to focus on investigating specific insertions at one time 
  # Any and all gaps in the tip sequence are deletions and need to be restored
  
  if(!grepl("-",tip)){
    return(c(tip,pos))
  }else{
    tip.chars <- strsplit(tip, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(tip.chars=="-")
    
    # perform a readjustment of the position (for insertions only 
    if (!is.na(pos)){
      pos <- gregexpr("[ACTG]", tip)[[1]][pos]
    }
    tip.chars[idx] <- anc.chars[idx]
    tip <- paste0(tip.chars,collapse="")
    return(c(tip, pos))
  }
}


# general
csvcount <- function(input,delim=","){
  if (is.na(input) || input == "-1"){
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
extractInfo <- function(input, delim=":"){
  if (is.na(input)){
    return(c(NA,NA))
  }else{
    insertions <- strsplit(input, delim)
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
truncate <- function(aaseq){
  if (grepl("\\*",aaseq)){
    idx <- gregexpr("\\*",aaseq)[[1]]-1
    aaseq <- substr(aaseq, 1, idx)
  }
  aaseq
}


# takes in an amino acid sequence and returns the locations of all Nglycs
extractGlycs <- function(aaseq, vect=F){
  aaseq <- truncate(aaseq)
  result <- gregexpr("N[^P][ST][^P]", aaseq)[[1]]  # used for 0 indexing these position values for analysis in python 
  if(!vect){
    paste(result, collapse=",")
  }else{
    result
  }
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
