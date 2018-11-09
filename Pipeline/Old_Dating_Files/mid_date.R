mid.date <- function(lo, hi){
  mid <- as.Date(lo) + as.integer(as.Date(hi)-as.Date(lo))/2
  return (mid)
}
