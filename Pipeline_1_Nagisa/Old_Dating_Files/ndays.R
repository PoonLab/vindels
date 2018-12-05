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
