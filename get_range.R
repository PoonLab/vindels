get.range <- function(x) {
  x = as.character(x)
  items <- strsplit(x, '-')[[1]]
  year <- items[1]
  if (length(items) == 1) {
    ## year only
    low <- as.Date(paste(year, '01', '01', sep='-'))
    high <- as.Date(paste(year, '12', '31', sep='-'))
  } 
  else if (length(items) == 2) {
    ## year and month
    month <- items[2]
    
    # determine number of days in this month
    days <- n.days(year, month)
    low <- as.Date(paste(year, month, '01', sep='-'))
    high <- as.Date(paste(year, month, days, sep='-'))
  } 
  else {
    # year, month, day
    days <- items[3]
    month <- items[2]
    low <- as.Date(paste(year, month, days, sep='-'))
    high <- low
  }
  return (c(low, high))
}
