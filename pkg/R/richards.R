#####################################################################################################################
## Equation 15. Richards equation (Richards, 1959) J Exp Bot 10: 290-300
## 4-parameter  
## tmax is the time when growth rate is maximum and RGR = theta/(1+v),
## thetha determines the curvature, Yo is the max Y
## v parameter determines the asymetrical growth. 
## If v=1, then Richards = logistic. If v=0, the Richards = Gompertz (see below Eq. 16)  
## Richards equation is an extension of the logistic to deal with asymmetrical growth.
## Applied to describe water stress index, crop N-uptake, seed germination, crop growth, LAI development, etc.
## Remark at time = 0, Y is not zero 
## Richeards equation has often been criticized because the shape parameter v has no obvious biological interpretation


rich  <- function(time, tmax, theta, v, Yo){
  ans  <- Yo/((1 + v*exp(-1*theta*(time - tmax)))^(1/v))
  return(ans)
# Y  <- Yo/((1 + v*exp(theta*(time - tmax)))^(1/v))       # by deleting the "-1" the shape is opposite
}

#####################################################################################################################
## Equation 16. Gompertz equation (Gompertz, 1825) 
## 3-parameter  
## tmax is the time when crop growth rate is maximum and Y = Yo/exp
## thetha determines the curvature, Yo is the max Y
## An alternative way to Richards equation to generate asymetric growth with less parameters
## Applied to describe water stress index, crop N-uptake, seed germination, crop growth, LAI development, etc.
## Remark at time = 0, Y is not zero


gompertz  <- function(time, tmax, theta, Yo){
  ans  <- Yo*exp(-1*exp(-1*theta*(time - tmax)))
  return(ans)

}

