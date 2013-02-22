## Calculating the mean bias from the 0 1 line

mbias <- function(obs,sim){

  stopifnot(length(obs) == length(sim))
  stopifnot(any(is.na(obs)!=TRUE) || any(is.na(sim)!=TRUE))

  x <- mean(obs - sim)

  x

}

rmse <- function(obs,sim){

  stopifnot(length(obs) == length(sim))
  stopifnot(any(is.na(obs)!=TRUE) || any(is.na(sim)!=TRUE))

  n <- length(obs)

  mse <- 1/n * sum(I(obs - sim)^2)

  list(mse=mse, rmse = sqrt(mse))

}
