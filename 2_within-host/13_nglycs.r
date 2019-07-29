require(ape)
require(stringr)
require(Biostrings)
removeNA <- function(input, repl=""){
  if (is.na(input)){
    input <- repl
  }
  input
}

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

extractGlycs <- function(aaseq){
  result <- gregexpr("N[^P][ST][^P]", aaseq)[[1]] * 3 - 2
  result
}

countGlycs <- function(field){
  if ("," %in% field){
    return(str_count(field, ",") + 1)
  }else if (field == ""){
    return(0)
  }else{
    return(1)
  }
}

insAlign <- function(indels, pos, ancestor, seq){
  i.list <- str_split(indels, ",")[[1]]
  p.list <- str_split(pos, ",")[[1]]
  
  for (idx in 1:length(i.list)){
    
    len <- nchar(i.list[idx])
    ix <- i.list[idx]
    px <- as.numeric(p.list[idx])
    
    ancestor <- paste0(substr(ancestor, 0, px-len), paste(rep("-", len),collapse=""), substr(ancestor,px-len+1, nchar(ancestor)))
  }
  
  ancestor
}

delAlign <- function(indels, pos, ancestor, seq){
  i.list <- str_split(indels, ",")[[1]]
  p.list <- str_split(pos, ",")[[1]]
  p.list <- as.numeric(p.list)
  
  for (idx in 1:length(i.list)){
    len <- nchar(i.list[idx])
    ix <- i.list[idx]
    px <- p.list[idx]
    
    seq <- paste0(substr(seq, 0, px), paste(rep("-", len),collapse=""), substr(seq,px+1, nchar(seq)))
    p.list[(idx+1):length(p.list)] <- p.list[(idx+1):length(p.list)] + len
  }
  
  seq
}



#CAAGGGATGGAGGAAAAAACAATACGGAGACATTCAGACCT
#PycharmProjects/hiv-withinhost/
path <- "~/PycharmProjects/hiv-withinhost/"
ins <- read.csv(paste0(path, "13_nglycs/ins.csv"),  sep="\t", stringsAsFactors = F)
del <- read.csv(paste0(path,"13_nglycs/del.csv"), sep="\t", stringsAsFactors = F)

headers <- c("accno", "vloop", "seq", "pos", "ancestor", "tipseq")
colnames(ins) <- headers
colnames(del) <- headers

#ins$Seq <- sapply(ins$Seq, translate)
#del$Seq <- sapply(del$Seq, translate)

new.ins <- ins[,c(1:4)]
new.del <- del[,c(1:4)]

new.ins$ancestor <- unname(mapply(insAlign, ins$seq, ins$pos, ins$ancestor, ins$tipseq))
new.ins$tipseq <- ins$tipseq

new.del$ancestor <- del$ancestor
new.del$tipseq <- unname(mapply(delAlign, del$seq, del$pos, del$ancestor, del$tipseq))

new.ins$anc.glycs <- sapply(sapply(ins$ancestor, translate), extractGlycs)
new.del$anc.glycs <- sapply(sapply(del$ancestor, translate), extractGlycs)

new.ins$tip.glycs <- sapply(sapply(ins$tipseq, translate), extractGlycs)
new.del$tip.glycs <- sapply(sapply(del$tipseq, translate), extractGlycs)

checkGlycs <- function(list1, list2){
  
}

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
