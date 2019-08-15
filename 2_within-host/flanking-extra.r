
m.before <- matrix(nrow=2, ncol=10)
m.after <- matrix(nrow=2,ncol=10)
for (wobble in 0:1){
  
  for (offset in 0:9){
    flanking <- unname(mapply(insCheck, indel=ins$Seq, pos=ins$Pos, vseq=ins$Vseq, wobble=wobble, offset=offset))
    flanking <- as.data.frame(t(flanking))
    flanking <- cbind( ins[,c(1,7)], len=nchar(ins$Seq), flanking)
    colnames(flanking) <- c("accno","pos", "len", "indel", "vseq","before.bool", "before.offset", "before.diff", "before.seq","after.bool", "after.offset", "after.diff",  "after.seq")    
    flanking$before.bool <- as.logical(flanking$before.bool)
    flanking$after.bool <- as.logical(flanking$after.bool)
    
    
    before.prop <- sum(flanking$before.bool) / nrow(flanking)
    after.prop <- sum(flanking$after.bool) / nrow(flanking)
    
    m.before[wobble+1, offset+1] <- before.prop
    m.after[wobble+1,offset+1] <- after.prop
  }
}

rownames(m.before) <- 0:9
colnames(m.before) <- 0:9
rownames(m.after) <- 0:9
colnames(m.after) <- 0:9


# K S TEST
# ------------------------------------------------------------
# path <- '~/PycharmProjects/hiv-withinhost/'
# ins <- read.csv(paste0(path,"10_nucleotide/ins.csv"), stringsAsFactors = F, row.names = 1)
# 
# rownames(ins) <- 1:nrow(ins)
# 
# 
# 
# 
# # randomly rearrange the same nucleotides found in the vseq 
# 
# null_dist <- c()
# for (i in 1:nrow(ins)){
#   seq <- str_split(ins$Seq[i], ",")[[1]]
#   pos <- str_split(ins$Pos[i], ",")[[1]]
#   vseq <- ins$ins.unchanged[i]
#   letters <- str_split(vseq, "")[[1]]
#   
#   for (s in 1:length(seq)){
#     if (nchar(seq[s]) > 3){
#     px <- as.numeric(pos[s]) - nchar(seq[s]) +1
#     for (j in 1:1000){
#       rseq <- paste(letters[sample(1:nchar(vseq))], collapse = "")
#       res <- gregexpr(seq[s],rseq)[[1]]
#       if (res != -1){
#         res <- res - px
#         null_dist <- c(null_dist, res)
#       }
#     }
#     }
#   }
# }
# 
# test_dist <- c()
# findMatch <- function(indel, pos, vseq){
#   seq <- str_split(indel, ",")[[1]]
#   pos <- str_split(pos, ",")[[1]]
#   
#   res <- c()
#   for (s in 1:length(seq)){
#     if (nchar(seq[s]) >= 3){
#       px <- as.numeric(pos[s]) - nchar(seq[s]) +1
#       #print(seq[s])
#       #print(vseq)
#       match <- gregexpr(seq[s],vseq)[[1]]
#       if (length(match) > 1 | match != -1){
#         match <- match - px
#         res <- c(res, match)
#       }
#     }
#   }
#   res
# }
# 
# testDist <- unname(unlist(mapply(findMatch, ins$Seq, ins$Pos, ins$Vseq)))
# 
# hist(testDist, breaks=seq(min(testDist)-0.5,max(testDist)+0.5), col='red', xlab="Length (nt)")


