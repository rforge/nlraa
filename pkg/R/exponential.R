#########################################################################################################
## Equation 1. First order exponential (see text for description) 

exp1 <- function(time, Yo, rate){
  ans <- Yo * exp(-1*rate*time)        # Decrease over time, e.g. mass loss due to decomposition process
  return(ans)
#  Y <- Yo *(1-exp(-rate*time))        # Increase over time, e.g. cumulative CO2 production during decomposition process 
 }

exp2 <- function(time, Yo, rate){
  ans <- Yo *(1-exp(-rate*time))        # Increase over time, e.g. cumulative CO2 production during decomposition process
  return(ans)
}


