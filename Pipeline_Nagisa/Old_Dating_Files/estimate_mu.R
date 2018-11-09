estimate.mu <- function (t, node.dates, p.tol = 0.05) 
{
  g <- glm(node.depth.edgelength(t)[1:length(node.dates)] ~ 
             node.dates, na.action = na.omit)
  null.g <- glm(node.depth.edgelength(t)[1:length(node.dates)] ~ 
                  1, na.action = na.omit)
  if ((1 - pchisq(AIC(null.g) - AIC(g) + 2, df = 1)) > p.tol) {
    warning(paste("Cannot reject null hypothesis (p=", (1 - 
            pchisq(AIC(null.g) - AIC(g), df=1)), ")"))
  }
  coef(g)[[2]]
}