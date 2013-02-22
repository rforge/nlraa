## This is just a logistic growth function
logis <- function(x, Asym, xmid, scal){

  ans <- Asym / (1 + exp((xmid - x)/scal))

  return(ans)

}
