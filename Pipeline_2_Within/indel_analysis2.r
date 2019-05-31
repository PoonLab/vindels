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
  newheader <- str_split(as.character(header),"\\.")[[1]][1]
  newheader
}

removeNA <- function(input){
  if (all(is.na(input))){
    input <- ""
  }
  input
}

getAccno <- function(input){
  accno <- strsplit(input, "\\.")[[1]][5]
  accno
}


# INSERTION PARSING ----------
ifolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/ins_fix/*.csv")
dfolder <- Sys.glob("~/PycharmProjects/hiv-withinhost/9Indels/del_fix/*.csv")
treedir <- "~/PycharmProjects/hiv-withinhost/7SampleTrees/prelim/"
data.df <- data.frame()
ins.vr <- data.frame()
del.vr <- data.frame()
iTotal <- data.frame()
dTotal <- data.frame()
count <- 0
sequences <- list()
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
  #iCSV$Ins <- removeNA(iCSV$Ins)
  #dCSV$Del <- removeNA(dCSV$Del)
  iAccno <- unname(sapply(iCSV$Accno, getAccno))
  dAccno <- unname(sapply(dCSV$Accno, getAccno))
  
  
  iCSV$Subtype <- unname(sapply(iCSV$Accno, cutHeader))
  dCSV$Subtype <- unname(sapply(dCSV$Accno, cutHeader))
  
  
  # store the sequences from these two data frames for nucleotide analysis
  sequences$ins <- as.character(iCSV$Seq)
  sequences$del <- as.character(dCSV$Seq)
  
  # remove them as they arent needed for this analysis
  dCSV$Seq <- NULL
  iCSV$Seq <- NULL
  
  # creates the counts column
  iCSV$Count <- sapply(iCSV$Ins, csvcount) 
  dCSV$Count <- sapply(dCSV$Del, csvcount)
  
  tre <- read.tree(paste0(treedir, basename, ".tree.sample"))
  #tre$tip.label <- unname(sapply(tre$tip.label, function(x){strsplit(x, "_")[[1]][1]}))
  
  # adjusts the tre tip labels to match 
  tre$tip.label <- unname(sapply(tre$tip.label, function(x){strsplit(x,"_")[[1]][1]}))
  
  branches <- tre$edge.length[tre$edge[,2] <=Ntip(tre)]   #branches will match exactly with the tre$tip.label order
  iCSV$Date <- branches[match(iAccno, tre$tip.label)]
  dCSV$Date <- branches[match(dAccno, tre$tip.label)]
  

  
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

  
  colnames(iCSV) <- c("Accno","Vloop", "Vlength","Subtype", "Count","Date", "Seq", "Pos")
  colnames(dCSV) <- c("Accno", "Vloop", "Vlength","Subtype", "Count","Date", "Seq", "Pos")
  
  iTotal <- rbind(iTotal, cbind(iCSV, sequences[[1]]))
  dTotal <- rbind(dTotal, cbind(dCSV, sequences[[2]]))
  
  # iTemp <- split(iCSV, iCSV$Vloop)
  # dTemp <- split(dCSV, dCSV$Vloop)

  # for (i in 1:5){
  #   ins.vr[i][[1]] <- rbind(ins.vr[i][[1]], iTemp[i][[1]])
  #   del.vr[i][[1]] <- rbind(del.vr[i][[1]], dTemp[i][[1]])
  # }
}
iTotal[,8] <- as.character(iTotal[,8])
dTotal[,8] <- as.character(dTotal[,8])

# for (i in 1:5){
#   write.csv(iTotal[iTotal$Vloop==i,c(1,2,6,8)], paste0("~/PycharmProjects/hiv-withinhost/10_nucleotides/ins/Ins_V",i,".csv"))
#   write.csv(dTotal[dTotal$Vloop==i,c(1,2,6,8)], paste0("~/PycharmProjects/hiv-withinhost/10_nucleotides/del/Del_V",i,".csv"))
# 
# }


require(BSDA)
# RATE ANALYSIS -------------
irates <- c()
drates <- c()
vlengths <- c(84,156,105,90,33)
all.df <- data.frame()
for (vloop in 1:5){
  iData <- iTotal[iTotal$Vloop==vloop,]
  dData <- dTotal[dTotal$Vloop==vloop,]
  
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



vloops <- c("V1","V2","V3","V4","V5")

insrates <- data.frame(VLoop=vloops, iRate=irates, AdjRate=irates*10^3)
delrates <- data.frame(VLoop=vloops, dRate=drates, AdjRate=drates*10^3)

indels <- cbind(insrates, delrates[,c(2,3)])

require(Rmisc)
require(ggplot2)

g1 <- ggplot(insrates, aes(x=VLoop, y=AdjRate,width=1)) + 
  geom_bar(colour="black", stat="identity",fill="dodgerblue",position="dodge",show.legend=F) + 
  labs(x="Variable Loop", 
       y=expression(paste("       Insertion Rate \n(Events/Nt/Year x  ", 10^-3, ")", sep = "")))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 5))+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 42, r = 10, b = 8, l = 18, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 3, b = 0, l = 12)),
        axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position="none")+ geom_text(aes(y=0.3,x=3 ),
                                           label="N/A", 
                                           size=6)


g2 <- ggplot(delrates, aes(x=VLoop, y=AdjRate,width=1)) + 
  geom_bar(colour="black", stat="identity",fill="firebrick1",position="dodge",show.legend=F) + 
  labs(x="Variable Loop", 
       y="Deletion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_reverse(lim=c(5,0))+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = -10, r = 10, b = 8, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=18,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=18,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14),
        legend.position="none") + geom_text(aes(y=0.3,x=3 ),
                                            label="N/A", 
                                            size=6)
  
multiplot(g1,g2)



# SUBTYPE STRATIFICATION 
# ----------------------------------------------
ins <- list()
ins$B <- iTotal[iTotal$Subtype=="B",]
ins$C <- iTotal[iTotal$Subtype=="C",]
del <- list()
del$B <- dTotal[dTotal$Subtype=="B",]
del$C <- dTotal[dTotal$Subtype=="C",]

sub_irates <- data.frame()
sub_drates <- data.frame()

for (i in 1:2){
  
  iData <- ins[[i]]
  dData <- del[[i]]
  irates <- c()
  drates <- c()
  
  st <- c("B", "C")
  
  for (vloop in 1:5){
    iFinal <- iData[iData$Vloop==vloop & iData$Date < 325 & iData$Count < 2,]
    dFinal <- dData[dData$Vloop==vloop & dData$Date < 325 & dData$Count < 2,]
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
  sub_irates <- rbind(sub_irates, data.frame(rate=irates, adjrate=(irates*10^3),subtype=rep(st[[i]],5), VLoop=vloops))
  sub_drates <- rbind(sub_drates, data.frame(rate=drates, adjrate=(drates*10^3),subtype=rep(st[[i]],5), VLoop=vloops))
}


# SUBTYPE PLOT 
g1 <- ggplot(sub_irates, aes(x=VLoop, y=adjrate,fill=subtype)) + 
  geom_bar(stat="identity",position="dodge") + 
  labs(title="Insertions",x="Variable Loop", 
       y="Insertion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 2))+
  scale_fill_brewer(palette="Set1")+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 12, r = 12, b = 15, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=20,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=20,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=25,hjust=0.5),
        legend.position="right", legend.text=element_text(size=18), legend.title=element_text(size=20)) + geom_text(aes(y=0.1,x=3 ),
                                            label="N/A", 
                                            size=6)
g1


# SUBTYPE PLOT 
g2 <- ggplot(sub_drates, aes(x=VLoop, y=adjrate,fill=subtype)) + 
  geom_bar(stat="identity",position="dodge") + 
  labs(title="Deletions",x="Variable Loop", 
       y="Deletion Rate")+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 6))+
  scale_fill_brewer(palette="Set1")+
  theme(panel.grid.major.y = element_line(color="black",size=0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing=unit(1, "mm"),
        #panel.background=element_rect(fill="gray88",colour="white",size=0),
        plot.margin =margin(t = 12, r = 12, b = 15, l = 24, unit = "pt"),
        axis.line = element_line(colour = "black"), 
        axis.title.y=element_text(size=20,margin=margin(t = 0, r = 11, b = 0, l = 6)),
        axis.title.x=element_text(size=20,margin=margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size=16),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=25,hjust=0.5),
        legend.position="right", legend.text=element_text(size=18), legend.title=element_text(size=20)) + geom_text(aes(y=0.1,x=3 ),
                                                                                                                    label="N/A", 
                                                                                                                    size=6)
g2

# LENGTH ANALYSIS
# --------------------------------------------------------

categorize <- function(count){
  if (count <= 3){
    "1-3"
  }else if(count > 3 & count <= 6){
    "4-6"
  }else{
    "7+"
  }
}

iLength <- iTotal[,c(2,4,7)]
dLength <- dTotal[,c(2,4,7)]

iLength$Len <- sapply(iLength$Seq, nchar)
dLength$Len <- sapply(dLength$Seq, nchar)

iLength <- iLength[iLength$Len!=0,]
dLength <- dLength[dLength$Len!=0,]


iLength$Bin <- sapply(iLength$Len,categorize)
dLength$Bin <- sapply(dLength$Len,categorize)

iLength$Bin <- factor(iLength$Bin,levels=c("7+","4-6","1-3"))
dLength$Bin <- factor(dLength$Bin,levels=c("7+","4-6","1-3"))


colnames(iLength) <- c("Variable Loop", "Subtype", "Seq", "Len", "Insertion Length (nt)")

mosaic(~Vloop + Bin, data=iLength,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=20),
                            gp_varnames=gpar(fontsize=26),
                            set_varnames = c(Vloop="Variable Loop", 
                                             Bin="Insertion Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(Vloop=c("V1","V2","V4","V5")))


mosaic(~Vloop + Bin, data=dLength,
       shade=T, main=NULL, direction="v",
       spacing=spacing_equal(sp = unit(0.7, "lines")),
       residuals_type="Pearson",
       margins=c(2,2,6,2),
       labeling_args = list(tl_labels = c(F,T), 
                            tl_varnames=c(F,T),
                            gp_labels=gpar(fontsize=20),
                            gp_varnames=gpar(fontsize=26),
                            set_varnames = c(Vloop="Variable Loop", 
                                             Bin="Deletion Length (nt)"),
                            offset_labels=c(0,0,0,0),rot_labels=c(0,0,0,0), just_labels=c("center","center","center","center")),
       legend=legend_resbased(fontsize = 20, fontfamily = "",
                              x = unit(0.5, "lines"), y = unit(2,"lines"),
                              height = unit(0.8, "npc"),
                              width = unit(1, "lines"), range=c(-10,10)),
       set_labels=list(Vloop=c("V1","V2","V4","V5")))



ggplot()
ggplot(all.df, aes(x=Date, y=Count, group=Count)) + geom_density_ridges(colour="white", fill="blue", scale=1, bandwidth=5)

# Get raw insertion counts 
sum(vr.df[[1]]$Count)
sum(vr.df[[5]]$Count)


