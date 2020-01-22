require(ape)
require(stringr)
require(Biostrings)

source("~/GitHub/vindels/2_within-host/utils.r")

removeDeletionsTip <- function(vseq, anc){
  if(!grepl("-",vseq)){
    return(vseq)
  }else{
    tip.chars <- strsplit(vseq, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(tip.chars=="-")
    
    for (c in idx){
      tip.chars[c] <- anc.chars[c]
    }
    res <- paste0(tip.chars,collapse="")
    res
  }
}

removeOtherGaps <- function(anc, tip, indel, pos){
  # find the location of all gap characters in tip/anc
  gaps <- gregexpr("-",anc)[[1]]
  
  # create a vector of positions to be copied over
  idx <- rep(T, nchar(tip))
  idx[gaps] <- F
  
  # retrieve the boundaries of the indel
  end <- as.numeric(pos)
  start <- as.numeric(pos) - nchar(indel) + 1
  
  # make sure it ignores the indel sequence
  toIgnore <- start:end
  idx[toIgnore] <- T

  anc.chars <- strsplit(anc, "")[[1]]
  anc <- paste0(anc.chars[idx], collapse="")
  
  return (anc)
}


restoreInsGaps <- function(vseq, anc, indel, pos){
  if(!grepl("-",vseq)){
    return(c(vseq,pos))
  }else{
    tip.chars <- strsplit(vseq, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(tip.chars=="-")
    
    
    # perform a readjustment of the position (for insertions only )
    if (any(idx < pos)){
      pos <- as.character(pos + sum(idx < pos))
    }
    
    tip.chars[idx] <- anc.chars[idx]
    tip <- paste0(tip.chars,collapse="")
    return(c(tip, pos))
  }
}

restoreDelGaps <- function(vseq, anc, indel){
  if(!grepl("-",anc)){
    return(anc)
  }else{
    tip.chars <- strsplit(vseq, "")[[1]]
    anc.chars <- strsplit(anc, "")[[1]]
    idx <- which(anc.chars=="-")
    
    anc.chars[idx] <- tip.chars[idx]
    anc <- paste0(anc.chars,collapse="")
    return(anc)
  }
}



randomizationTest <- function(indel, seq){
  
  # PREPARATION 
  
  # remove all deletions in the tip sequences 
  # generate the ancestor sequence ONLY containing a single insertion
    # remove the dashes from the proper insertion
    # restore all other insertions in the sequence using 
  
  
  
  # INSERTIONS 
  
  
  # calc the size of the insertion 
  len <-  nchar(indel)
  # determine the locations of all N-glyc sites in the ancestral sequence 
  
  # create 100x random numbers (sample) between 1 and nchar(anc) - len(insertion) 
  # this will be the start point of the test insertion 
  
  # for every number in this random sample 
  
    # generate the result sequence ; use substring to add the insertion sequence into the ancestor 
  
    # recalculate the number of N-glyc sites 
    # adjust the location of all N-glyc locations falling AFTER the random position to check their similarity 
  
    # report this as a positive or a negative result zzzz
  
  
  
}


#PycharmProjects/hiv-withinhost/
path <- "~/PycharmProjects/hiv-withinhost/"
path <- "~/Lio/"
ins <- read.csv(paste0(path, "13_nglycs/ins-sep.csv"),  sep="\t", stringsAsFactors = F)
del <- read.csv(paste0(path,"13_nglycs/del-sep.csv"), sep="\t", stringsAsFactors = F)
ins <- ins[,-c(3,4)]
del <- del[,-c(3,4)]

ins$Vpos <- NULL
del$Vpos <- NULL

# apply an adjust to the deletion locations to make them the same as insertions 
del$Pos <- as.numeric(del$Pos) + nchar(del$Seq)

# Insertions : fill in gaps found in the tip sequences 
res <- as.data.frame(t(unname(mapply(restoreInsGaps,ins$Vseq, ins$Anc, ins$Seq,ins$Pos))))
ins$Vseq <- as.character(res[,1])
ins$Pos <- as.numeric(as.character(res[,2]))


# Deletions : fill in gaps found in the ancestral sequences 
del$Anc <- unname(mapply(restoreDelGaps, del$Anc, del$Vseq, del$Seq))

# Insertions : 
ins$Anc <- unname(mapply(removeOtherGaps, ins$Anc,ins$Vseq, ins$Seq, ins$Pos))
del$Vseq <- 
#headers <- c("accno", "vloop", "indel", "pos", "tip","anc", "patient")
#colnames(ins) <- headers
#colnames(del) <- headers

new.ins <- ins
new.del <- del

# new.ins$anc <- unname(mapply(insAlign, ins$seq, ins$pos, ins$anc, ins$tip))
# new.ins$tip <- ins$tip
# 
# new.del$anc <- del$anc
# new.del$tip <- unname(mapply(delAlign, del$seq, del$pos, del$anc, del$tip))

new.ins$anc.glycs <- unlist(sapply(sapply(ins$Anc, translate), extractGlycs))
new.del$anc.glycs <- unlist(sapply(sapply(del$Anc, translate), extractGlycs))

new.ins$tip.glycs <- unlist(sapply(sapply(ins$Tip, translate), extractGlycs))
new.del$tip.glycs <- unlist(sapply(sapply(del$Tip, translate), extractGlycs))


write.table(new.ins, paste0(path,"13_nglycs/ins-edit.csv"), sep="\t", quote=F, row.names=F)
write.table(new.del, paste0(path,"13_nglycs/del-edit.csv"), sep="\t", quote=F, row.names=F)



ins$aaseq <- NULL
del$aaseq <- NULL

ins$original <- mapply(insOriginal, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq)
del$original <- mapply(delOriginal, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq)
# for both: 
# determine the start and stop of all nglycs 
# collect them in one column separated by "-", comma separated


# insertions 
# determine what the original sequence is 
# start = pos - nchar(insertion sequence)
# end = pos 

# original seq: insertion removed  
# original <- paste0(substring(vloop, 0, start), substring(vloop, end, nchar(vloop)))


# deletions 
# start = pos - nchar(insertion sequence)
# end = pos 

# original seq: deletion added back in 
# original <- paste0(substring(vloop,0,start), del, substring(vloop,end,nchar(vloop)))

# --------------- OUTDATED

# insAlign <- function(indels, pos, anc, seq){
#   i.list <- str_split(indels, ",")[[1]]
#   p.list <- str_split(pos, ",")[[1]]
#   
#   for (idx in 1:length(i.list)){
#     
#     len <- nchar(i.list[idx])
#     ix <- i.list[idx]
#     px <- as.numeric(p.list[idx])
#     
#     anc <- paste0(substr(anc, 0, px-len), paste(rep("-", len),collapse=""), substr(anc,px-len+1, nchar(anc)))
#   }
#   
#   anc
# }


# 
# delAlign <- function(indels, pos, anc, seq){
#   i.list <- str_split(indels, ",")[[1]]
#   p.list <- str_split(pos, ",")[[1]]
#   p.list <- as.numeric(p.list)
#   
#   for (idx in 1:length(i.list)){
#     len <- nchar(i.list[idx])
#     ix <- i.list[idx]
#     px <- p.list[idx]
#     
#     seq <- paste0(substr(seq, 0, px), paste(rep("-", len),collapse=""), substr(seq,px+1, nchar(seq)))
#     p.list[(idx+1):length(p.list)] <- p.list[(idx+1):length(p.list)] + len
#   }
#   
#   seq
# }

