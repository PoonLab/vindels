handle.error <- function(tre, vect){
  rtdtree <- tryCatch(
    {rtdtree <- rtt(tre, vect)} ,
    error = function(c){
      message("Faulty rtt error, still running")
      return (NULL)
    }
  )
  return (rtdtree)
}

