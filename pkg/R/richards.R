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
  Y  <- Yo*exp(-1*exp(-1*theta*(time - tmax)))  

}

## Set parameter values and plot the relationship
time <- seq(0, 150,5)
ans1 <- Y(time, tmax = 60, theta = 0.09, Yo=100)
ans2 <- Y(time, tmax = 60, theta = 0.07, Yo=100)
ans3 <- Y(time, tmax = 60, theta = 0.05, Yo=100)
ans4 <- Y(time, tmax = 60, theta = 0.03, Yo=100)
xyplot(ans1 + ans2 + ans3 + ans4 ~ time, type="l", auto.key=TRUE, ylab = "text", xlab = "time")
