require(bbmle)
require(stringr)
require(ape)

csvcount <- function(input){
  commas <- str_count(input, ",")
  if (commas > 0){
    result <- commas + 1  
  }else if(input == ""){
    result <- 0
  }else{
    result <- 1
  }
  result
}

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
     pos <- c(pos, as.numeric(fields[[1]][2]))
  }
  return(c(paste(seq,collapse=","), paste(pos,collapse=",")))
}

cutHeader <- function(header){
  newheader <- str_split(as.character(header),"_")[[1]][1]
  newheader
}


# INSERTION PARSING ----------
#Sys.glob("~/Lio/deletions/*.csv")
ifolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/insertions/*.csv")
dfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/deletions/*.csv")
treedir <- "~/PycharmProjects/hiv-withinhost/7SampleTrees/prelim/"
data.df <- data.frame()
ins.vr <- list()
del.vr <- list()
iTotal <- data.frame()
dTotal <- data.frame()
count <- 0
for (file in 1:length(ifolder)){
  print(file)
  basename <- strsplit(basename(ifolder[file]),"\\.")[[1]][1]
  count <- count + 1
  iCSV <- read.csv(ifolder[file], stringsAsFactors = F)
  dCSV <- read.csv(dfolder[file], stringsAsFactors = F)
  if (all(is.na(iCSV$Ins))){
    iCSV$Ins <- ""
  }
  
  if (all(is.na(dCSV$Del))){
    dCSV$Del <- ""
  }
  iCSV$Count <- sapply(iCSV$Ins, csvcount) 
  dCSV$Count <- sapply(dCSV$Del, csvcount)
  
  tre <- read.tree(paste0(treedir, basename, ".tree.sample"))
  #tre$tip.label <- unname(sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][1]}))
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]   #branches will match exactly with the tre$tip.label order
  iCSV$Date <- branches[match(iCSV$Accno, tre$tip.label)]
  dCSV$Date <- branches[match(dCSV$Accno, tre$tip.label)]
  
  iCSV$Accno <- unname(sapply(iCSV$Accno, cutHeader))
  dCSV$Accno <- unname(sapply(dCSV$Accno, cutHeader))
  
  insInfo <- sapply(iCSV$Ins, extractInfo)
  insInfo <- unname(insInfo)
  insInfo <- t(insInfo)
  insInfo <- as.data.frame(insInfo)
  insInfo$V1 <- as.character(insInfo$V1)
  insInfo$V2 <- as.character(insInfo$V2)
  iCSV <- cbind(iCSV, insInfo)
  iCSV$Ins <- NULL
  
  delInfo <- sapply(dCSV$Del, extractInfo)
  delInfo <- unname(delInfo)
  delInfo <- t(delInfo)
  delInfo <- as.data.frame(delInfo)
  delInfo$V1 <- as.character(delInfo$V1)
  delInfo$V2 <- as.character(delInfo$V2)
  dCSV <- cbind(dCSV, delInfo)
  dCSV$Del <- NULL
  
  colnames(iCSV) <- c("Accno","Vloop", "Vlength", "Count","Date", "Seq", "Pos")
  colnames(dCSV) <- c("Accno", "Vloop", "Vlength", "Count","Date", "Seq", "Pos")
  
  iTotal <- rbind(iTotal, iCSV)
  dTotal <- rbind(dTotal, dCSV)
  
  iTemp <- split(iCSV, iCSV$Vloop)
  dTemp <- split(dCSV, dCSV$Vloop)
  
  
  for (i in 1:5){
    ins.vr[i][[1]] <- rbind(ins.vr[i][[1]], iTemp[i][[1]])
    del.vr[i][[1]] <- rbind(del.vr[i][[1]], dTemp[i][[1]])
  }
}



for (i in 1:5){
  write.csv(ins.vr[i][[1]], paste0("~/Lio/V",i,".csv"))
  write.csv(ins.vr[i][[1]], paste0("~/Lio/V",i,".csv"))
}

require(BSDA)
# RATE ANALYSIS -------------
irates <- c()
drates <- c()
vlengths <- c(84,156,105,90,33)
all.df <- data.frame()
for (vloop in 1:5){
  iData <- ins.vr[[vloop]]
  dData <- del.vr[[vloop]]
  
  iFinal <- iData[iData$Date < 325 & iData$Count < 2,]
  dFinal <- dData[dData$Date < 325 & dData$Count < 2,]
  #print(nrow(current) - nrow(iFinal))
  
  ifit <- glm(iFinal$Count ~ 1, offset=log(iFinal$Date), family="poisson")
  irate <- exp(coef(ifit)[[1]])*365/vlengths[vloop]
  irates <- c(irates, irate)
  print(summary(ifit))
  
  dfit <- glm(dFinal$Count ~ 1, offset=log(dFinal$Date), family="poisson")
  drate <- exp(coef(dfit)[[1]])*365/vlengths[vloop]
  drates <- c(drates, drate)
  print(summary(dfit))
  
  #print(1 - (fit$deviance/fit$null.deviance))
  #all.df <- rbind(all.df, iFinal)
  #EDA(residuals(fit))
  
  #to.remove <- c(order(residuals(fit), decreasing = T)[1:52])
  #iFinal2 <- iFinal[-to.remove,]
  #fit2 <- glm(iFinal2$Count ~ iFinal2$Date, family="poisson")
  #EDA(residuals(fit2))
}

indelrates <- data.frame(VLoop=c("V1","V2","V3","V4","V5"), Rate=rates)

require(ggplot2)

plot <- ggplot(indelrates, aes(x=VLoop, 
                            y=AdjRate,
                            width=1)) + geom_bar(colour="black",
                                                 stat="identity", 
                                                 fill="dodgerblue",
                                                 position="dodge", 
                                                 show.legend=F) 

plot <- plot + labs(x="Variable Loop", 
            y=expression(paste("Indel Rate (Events/Nt/Year x ", 10^-3, ")", sep = "")))+scale_fill_manual(values=colors2)+scale_y_continuous(expand = c(0, 0),
                                                                                                                                             limits = c(0, 1.5))+theme(panel.grid.major.y = element_line(color="black",size=0.3),
                                                                                                                                                                     panel.grid.major.x = element_blank(),
                                                                                                                                                                     panel.grid.minor.y = element_blank(),
                                                                                                                                                                     panel.grid.minor.x = element_blank(),
                                                                                                                                                                     panel.spacing=unit(1, "mm"),
                                                                                                                                                                     #panel.background=element_rect(fill="gray88",colour="white",size=0),
                                                                                                                                                                     plot.margin =margin(t = 20, r = 20, b = 20, l = 8, unit = "pt"),
                                                                                                                                                                     axis.line = element_line(colour = "black"), 
                                                                                                                                                                     axis.title.y=element_text(size=20,margin=margin(t = 0, r = 15, b = 0, l = 0)),
                                                                                                                                                                     axis.title.x=element_text(size=20,margin=margin(t = 15, r = 0, b = 0, l = 0)),
                                                                                                                                                                     strip.text.x = element_text(size=16),
                                                                                                                                                                     axis.text=element_text(size=14),
                                                                                                                                                                     legend.position="none")


ggplot()
ggplot(all.df, aes(x=Date, y=Count, group=Count)) + geom_density_ridges(colour="white", fill="blue", scale=1, bandwidth=5)

# Get raw insertion counts 
sum(vr.df[[1]]$Count)
sum(vr.df[[5]]$Count)


