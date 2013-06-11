## Eq. 1.2 N distribution within canopies (Jonshon et al., 2010 AoB)
## nb is the min N for photosynthesis
## kn is the N exctinction coefficient
## c emprical paramater, when c=0 then the there is no decay (uniform distribution)
## c emprical paramater, when c=1 then the decay is exponential
## c emprical paramater, when c>1 then the decay takes an inverse S shape

Nlai <- function(lai, Yo, nb, kn, c){
  ans <- Yo-(Yo-nb)*(1-exp(-kn*LAI))^c
  return(ans)
}


