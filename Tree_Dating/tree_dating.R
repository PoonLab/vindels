dfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/Dates_edit/",full.names=FALSE)
tfolder <- list.files(path="~/PycharmProjects/hiv-evolution-master/7_1Trees/", full.names=FALSE)
setwd("~/PycharmProjects/hiv-evolution-master/8_Dated_Trees")

names <- c()

for (i in 88:length(tfolder)){
  tre <- read.tree(paste("~/PycharmProjects/hiv-evolution-master/7_1Trees",tfolder[i],sep="/"))
  csv <- read.csv(paste("~/PycharmProjects/hiv-evolution-master/Dates_edit",dfolder[i],sep="/"), header=FALSE)
  filename <- strsplit(tfolder[i], "\\+")[[1]][1]
  print(filename)
  
  tre2 <- multi2di(tre) 
  
  #loading a temporary data frame of mid dates
  names(csv) <- c('accno', 'date')
  temp <- NULL
  for (x in 1:nrow(csv)){
    dtx <- get.range(csv$date[x])
    middate <- 1970 + as.double(as.Date(mid.date(dtx[1],dtx[2])))/365.25
    temp <- rbind(temp, data.frame(accno=csv$accno[x],date=middate))
  }
  #rearranges the temp file to be in the same order as the tree
  index <- match(tre2$tip.label, temp$accno)
  temp <- temp[index,]
  
  #loading the vector of dates
  vect <- c()
  for (n in 1:nrow(temp)){
    vect[n] <- temp$date[n]
  }
  
  #rooting and date adjustments
  rtdtree <- handle.error(tre2, vect)
  print(rtdtree)
  if (is.null(rtdtree)){
    next
  }
  
  mu <- handle.warning(rtdtree, vect)
  print(mu)
  if (is.null(mu)){
    next
  }
  else {
    node.date <- estimate.dates(rtdtree, vect, mu)
    rtdtree$edge.length <- node.date[rtdtree$edge[, 2]] - node.date[rtdtree$edge[, 1]]
    write.tree(rtdtree, file=paste0(filename, "+.tree"))
    
  }
  
}







#--------------------------------------------TEST AREA
# t to transpose matrix result
# temp <- t(sapply(csv1$date, function(x) get.range(x)))
# csv1$min.date <- as.Date(temp[,1], origin='1970-01-01')
# csv1$max.date <- as.Date(temp[,2], origin='1970-01-01')



#vect is a vector that I am loading with csv 

# fn <- which(grepl("/A1_", tfolder))
# tre1 <- read.tree(tfolder[fn])
# index <- match(tre1$tip.label, temp$accno)
# temp <- temp[index,]


#rearranges the temp file to be in the same order as the tree

# names(vec) <- stss[,1]
# vec <- stss[,2]