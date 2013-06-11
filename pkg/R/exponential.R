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

########################################################################################################
## Equation 2. Parallel or double first order exponential (see text for description) 

exp3 <- function(time, Yo, rate1, rate2, fraction){
  ans <- Yo * (fraction* exp(-rate1*time) + (1-fraction)*exp(-rate2*time))       # Decreasing pattern
  return(ans)
}

exp4 <- function(time, Yo, rate1, rate2, fraction){
   ans <- Yo- (Yo * (fraction* exp(-rate1*time) + (1-fraction)*exp(-rate2*time)))  # Increasing  pattern
   return(ans)
}





