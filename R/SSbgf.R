## Implementing the beta growth function from (Yin et al 2003)

bgfInit <- function(mCall, LHS, data){

  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a bgf")
  }
  w.max <- max(xy[,"y"])
  t.e <- NLSstClosestX(xy, w.max)
  t.m <- t.e / 2
  value <- c(w.max, t.e, t.m)
  names(value) <- mCall[c("w.max","t.e","t.m")]
  value

}


bgf <- function(time, w.max, t.e, t.m){

  .expr1 <- t.e / (t.e - t.m)
  .expr2 <- (time/t.e)^.expr1
  .expr3 <- (1 + (t.e - time)/(t.e - t.m))
  .value <- w.max * .expr3 * .expr2

  ## Derivative with respect to t.e
  .exp1 <- ((time/t.e)^(t.e/(t.e - t.m))) * ((t.e-time)/(t.e-t.m) + 1)
  .exp2 <- (log(time/t.e)*((1/(t.e-t.m) - (t.e/(t.e-t.m)^2) - (1/(t.e - t.m)))))*w.max
  .exp3 <- (time/t.e)^(t.e/(t.e-t.m))
  .exp4 <- w.max * ((1/(t.e-t.m)) - ((t.e - time)/(t.e-t.m)^2))
  .exp5 <- .exp1 * .exp2 + .exp3 * .exp4 

  ## Derivative with respect to t.m
  .ex1 <- t.e * (time/t.e)^((t.e/(t.e - t.m))) * log(time/t.e) * ((t.e - time)/(t.e - t.m) + 1) * w.max
  .ex2 <- (t.e - time) * w.max * (time/t.e)^(t.e/(t.e-t.m))
  .ex3 <- (t.e - t.m)^2
  .ex4 <- .ex1 / .ex3 + .ex2 / .ex3
  
  .actualArgs <- as.list(match.call()[c("w.max", "t.e", "t.m")])

##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("w.max", 
                                                          "t.e", "t.m")))
    .grad[, "w.max"] <- .expr3 * .expr2
    .grad[, "t.e"] <- .exp5
    .grad[, "t.m"] <- .ex4 
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
    .value
}

SSbgf <- selfStart(bgf, initial = bgfInit, c("w.max", "t.e", "t.m"))

## Beta growth initial growth

bgf2 <- function(time, w.max, w.b, t.e, t.m, t.b){

  .expr1 <- (t.e - t.b) / (t.e - t.m) 
#  .expr11 <- pmax(c(0, time - t.b))
  .expr11 <- (time - t.b) 
  .expr2 <- .expr11/(t.e-t.b)
  .expr3 <- .expr2 ^ (.expr1) 
  .expr4 <- 1 + (t.e - time)/(t.e - t.m)
  .value <- w.b + (w.max - w.b) * .expr4 * .expr3

  .value[is.nan(.value)] <- 0
  
##   ## Derivative with respect to t.e
##   .exp1 <- ((time/t.e)^(t.e/(t.e - t.m))) * ((t.e-time)/(t.e-t.m) + 1)
##   .exp2 <- (log(time/t.e)*((1/(t.e-t.m) - (t.e/(t.e-t.m)^2) - (1/(t.e - t.m)))))*w.max
##   .exp3 <- (time/t.e)^(t.e/(t.e-t.m))
##   .exp4 <- w.max * ((1/(t.e-t.m)) - ((t.e - time)/(t.e-t.m)^2))
##   .exp5 <- .exp1 * .exp2 + .exp3 * .exp4 

##   ## Derivative with respect to t.m
##   .ex1 <- t.e * (time/t.e)^((t.e/(t.e - t.m))) * log(time/t.e) * ((t.e - time)/(t.e - t.m) + 1) * w.max
##   .ex2 <- (t.e - time) * w.max * (time/t.e)^(t.e/(t.e-t.m))
##   .ex3 <- (t.e - t.m)^2
##   .ex4 <- .ex1 / .ex3 + .ex2 / .ex3
  
##   .actualArgs <- as.list(match.call()[c("w.max", "t.e", "t.m")])

## ##  Gradient
##   if (all(unlist(lapply(.actualArgs, is.name)))) {
##     .grad <- array(0, c(length(.value), 3L), list(NULL, c("w.max", 
##                                                           "t.e", "t.m")))
##     .grad[, "w.max"] <- .expr3 * .expr2
##     .grad[, "t.e"] <- .exp5
##     .grad[, "t.m"] <- .ex4 
##     dimnames(.grad) <- list(NULL, .actualArgs)
##     attr(.value, "gradient") <- .grad
##   }
    .value
}
