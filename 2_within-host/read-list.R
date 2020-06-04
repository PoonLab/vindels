# slip list examination 

csv <- readLines("~/PycharmProjects/hiv-withinhost/15_modeling/list-21-full-stacks-100000.csv")

id <- which(unname(sapply(csv, function(x){
  res <- strsplit(x, "")[[1]]
  sum(sapply(res, as.numeric)) > 0
})))

csv[idx]

slip <- csv[which(sapply(csv, sum) > 0)]
