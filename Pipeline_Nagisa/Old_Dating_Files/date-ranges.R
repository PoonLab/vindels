dates <- read.table('NAME OF FILE')

# write function to generate new columns in data frame for start and end dates of range
# assume ISO-format dates in column $isodate

n.days <- function(year, month) {
  # @param year: string (number) or integer for full year, e.g., '1997'
  # @param month: string (number) or integer for month
  
  start <- as.Date(paste(year, month, '01', sep='-'))
  if (month == '12') {
    # increment to start of next year
    end <- as.Date(paste(as.integer(year)+1, '01', '01', sep='-'))
  } else {
    end <- as.Date(paste(year, as.integer(month)+1, '01', sep='-'))
  }
  return(as.integer(difftime(end, start)))
}

get.range <- function(x) {
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
    low <- as.Date(paste(year, month, days, sep='-'))
    high <- low
  }
  return (c(low, high))
}

