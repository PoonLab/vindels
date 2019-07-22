require(ape)
require(stringr)
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
  dnabin <- as.DNAbin.DNAString(dna)
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
CAAGGGATGGAGGAAAAAACAATACGGAGACATTCAGACCT
#PycharmProjects/hiv-withinhost/
ins <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/ins.csv", row.names=1, stringsAsFactors = F)
del <- read.csv("~/PycharmProjects/hiv-withinhost/13_nglycs/del.csv", row.names=1, stringsAsFactors = F)

#ins$Pos <- ins$Pos + 1
#del$Pos <- del$Pos + 1


#ins$Seq <- sapply(ins$Seq, translate)
#del$Seq <- sapply(del$Seq, translate)

ins$aaseq <- sapply(ins$Vseq, translate)
del$aaseq <- sapply(del$Vseq, translate)

ins$glycs <- sapply(ins$aaseq, extractGlycs)
del$glycs <- sapply(del$aaseq, extractGlycs)

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

