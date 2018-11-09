handle.warning <- function(rtdtree, vect) {
  mu <- tryCatch(
    {mu <- estimate.mu(rtdtree,vect)},
    warning = function(c){
      message("A warning was thrown, still running")
      return (NULL)
    })
  return(mu)
}
