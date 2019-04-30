require(bbmle)
require(stringr)
require(ape)
poisll <- function (rate, counts, times){
  prob <- (exp(-rate*times))*((rate*times)^(counts))/(factorial(counts))    # (e^-rt)*((rt)^k)/(k!)
  sum(log(prob))
}


csvcount <- function(input){
  commas <- str_count(input, ",")
  if (commas > 0){
    result <- commas + 1  
  }else if(input == ""){
    result <- 0
  }else{
    result <- 1
  }
}

extractInfo <- function(input){
  if (length(input)==1 && input == ""){
    return(c("",""))
  }else{
    insertions <- str_split(input, ",")
  }
  seq <- c()
  pos <- c()
  
  for (ins in insertions[[1]]){
     fields <- str_split(ins, "-")
     seq <- c(seq, fields[[1]][1])
     pos <- c(pos, as.numeric(fields[[1]][2]))
  }
  return(c(paste(seq,collapse=","), paste(pos,collapse=",")))
}

cutHeader <- function(header){
  newheader <- str_split(as.character(header),"_")[[1]][1]
  newheader
}


#TESTING ----------------------
#csvfile <- read.csv(file, sep="\t", stringsAsFactors = F)
#output <- sapply(csvfile$Ins, extractInfo)



# INSERTION PARSING ----------
folder <- Sys.glob("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/insertions/*.csv")
treedir <- "/home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/prelim/"
data.df <- data.frame()
vr.df <- list()
total.df <- data.frame()
count <- 0
for (file in folder){
  print(file)
  basename <- str_split(basename(file),"\\.")[[1]][1]
  count <- count + 1
  csvfile <- read.csv(file,sep="\t", stringsAsFactors = F)
  if (all(is.na(csvfile$Ins))){
    csvfile$Ins <- ""
  }
  csvfile$Count <- sapply(csvfile$Ins, csvcount) 
  
  tre <- read.tree(paste0(treedir, basename, ".tree.sample"))
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]   #branches will match exactly with the tre$tip.label order
  csvfile$Date <- branches[match(csvfile$Accno, tre$tip.label)]
  csvfile$Accno <- unname(sapply(csvfile$Accno, cutHeader))
  
  output <- sapply(csvfile$Ins, extractInfo)
  output <- unname(output)
  output <- t(output)
  output <- as.data.frame(output)
  output$V1 <- as.character(output$V1)
  output$V2 <- as.character(output$V2)
  csvfile <- cbind(csvfile, output)
  csvfile$Ins <- NULL
  
  colnames(csvfile) <- c("Accno", "Vloop", "Count","Date", "Seq", "Pos")
  total.df <- rbind(total.df, csvfile)
  temp.df <- split(csvfile, csvfile$Vloop)
  
  
  for (i in 1:5){
    vr.df[i][[1]] <- rbind(vr.df[i][[1]], temp.df[i][[1]])
  }
}
require(BSDA)
# RATE ANALYSIS -------------
rates <- c()
vlengths <- c(84,156,105,90,33)
all.df <- data.frame()
for (vloop in 1:5){
  current <- vr.df[[vloop]]
  finalized <- current[current$Date < 325 & current$Count < 2,]
  #print(nrow(current) - nrow(finalized))
  
  fit <- glm(finalized$Count ~ 1, offset=log(finalized$Date), family="poisson")
  rate <- exp(coef(fit)[[1]])*365/vlengths[vloop]
  rates <- c(rates, rate)
  print(summary(fit))
  #print(1 - (fit$deviance/fit$null.deviance))
  #all.df <- rbind(all.df, finalized)
  #EDA(residuals(fit))
  
  #to.remove <- c(order(residuals(fit), decreasing = T)[1:52])
  #finalized2 <- finalized[-to.remove,]
  #fit2 <- glm(finalized2$Count ~ finalized2$Date, family="poisson")
  #EDA(residuals(fit2))
}


ggplot(all.df, aes(x=Date, y=Count, group=Count)) + geom_density_ridges(colour="white", fill="blue", scale=1, bandwidth=5)

# Get raw insertion counts 
sum(vr.df[[1]]$Count)
sum(vr.df[[5]]$Count)


