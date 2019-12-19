ATcontent <- function(indel, pos, vseq, sample_length){
  pos <- as.numeric(pos)
  nucleotides <- c("A", "C","G","T")
  
  counts <- function(indel,pos,vseq){
    len <- nchar(indel)
    nucleotides <- c("A", "C","G","T")
    
    if (pos - sample_length >= 0){
      upstream <- substr(vseq, pos-sample_length, pos)
      #print(upstream)
      ucounts <- c()
      for (n in 1:4){
        ucounts[n] <- str_count(upstream,nucleotides[n])
      }
    }else{
      ucounts <- c(NA,NA,NA,NA)
    }
    
    if ((pos + len + sample_length) <= nchar(vseq)){
      downstream <- substr(vseq, pos+len+1, pos+len+sample_length)
      dcounts <- c()
      for (n in 1:4){
        dcounts[n] <- str_count(downstream,nucleotides[n])
      }
      
    }else{
      dcounts <- c(NA,NA,NA,NA)
    }
    c(ucounts,dcounts)
  }
  nucl <- unname(mapply(counts, indel, pos, vseq))
  up <- as.data.frame(t(nucl))[,1:4]
  print(nrow(up[is.na(up[,1]),]))
  print(nrow(up))
  up <- up[!is.na(up[,1]),]
  colnames(up) <- nucleotides
  
  down <- as.data.frame(t(nucl))[,5:8]
  print(nrow(down[is.na(down[,1]),]))
  print(nrow(down))
  down <- down[!is.na(down[,1]),]
  colnames(down) <- nucleotides
  
  data.frame(upstream=colSums(up)/(nrow(up)*sample_length),downstream=colSums(down)/(nrow(down)*sample_length))
  
}

flanking <- read.csv("~/PycharmProjects/hiv-withinhost/14_flanking/flanking.csv",row.names=1,stringsAsFactors = F)


# RANDOMIZATION TEST 
# -----------------------------------------
sampleString <- function(vloop, pos){
  idx <- sample(1:(nchar(vloop)-14),100, replace=TRUE)
  strings <- sapply(idx, function(x){substr(vloop, x, x+14)})
  props <- list()
  nucl <- c("A","C","G","T")
  for (n in 1:4){
    props[[n]] <- unname(sapply(strings, function(x){str_count(x, nucl[n])/nchar(x)}))
  }
  props
}

iSample <- list(numeric(nrow(flanking)*100),numeric(nrow(flanking)*100),numeric(nrow(flanking)*100),numeric(nrow(flanking)*100))

# generates the randomly sampled substrings for each indel
for (row in 1:nrow(flanking)){
  itemp <- sampleString(flanking[row,"vseq"])
  for (i in 1:4){
    iSample[[i]][((row-1)*100+1):(row*100)] <- itemp[[i]]
  }
}

# compares the observed proportion to the overall distribution of each nucleotide 
usign <- c()
dsign <- c()
u.pct <-c()
d.pct <-c()
for (i in 1:4){
  idist <- iSample[[i]]
  u.prop <- at.props[i,1]
  d.prop <- at.props[i,2]
  
  # find the percentile of the upstream and downstream AT content
  percentile <- ecdf(idist)
  u.pct[i] <- percentile(u.prop)
  d.pct[i] <- percentile(d.prop)
  
  # check for significance
  iQT <- quantile(idist, probs=c(0.025,0.975))
  
  # highlight significant differences 
  if (u.prop < iQT[[1]]){
    usign[i] <- "lower"
  }else if(u.prop > iQT[[2]]){
    usign[i] <- "higher"
  }else{
    usign[i] <- ""
  }
  
  # highlight significant differences 
  if (d.prop < iQT[[1]]){
    dsign[i] <- "lower"
  }else if(d.prop > iQT[[2]]){
    dsign[i] <- "higher"
  }else{
    dsign[i] <- ""
  }
  
}
print(u.pct)
ins.nt$sign <- isign
del.nt$sign <- dsign