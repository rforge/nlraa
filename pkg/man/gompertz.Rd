\name{gompertz}
\alias{gompertz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gompertz equation (Gompertz, 1825) 
  %%  ~~function to do ... ~~
}
\description{
  Gompertz equation (Gompertz, 1825). This is equation 2.3 in Table 1. 
}
\usage{
gompertz(time, tmax, theta, Yo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
a numeric vector of values at which to evaluate the model.
}
\item{tmax}{
  tmax is the time when crop growth rate is maximum and Y = Yo/exp
}
\item{theta}{
  theta determines the curvature
}
\item{Yo}{
   Yo is the max Y
}
}
\details{
An alternative way to Richards equation to generate asymetric growth with less parameters
Applied to describe water stress index, crop N-uptake, seed germination, crop growth, LAI development, etc.
Remark at time = 0, Y is not zero
}
\value{
  a numeric vector of length equal to the inputs
}
\references{
    Nonlinear Regression Models and Applications in Agricultural Research.
  Sotirios V. Archontoulis and
  Fernando E. Miguez. Agronomy Journal. doi: 10.2134/agronj2012.0506 
  %% ~put references to the literature/web site here ~
}
\author{
  Fernando E. Miguez
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{SSgompertz}}
}
\examples{
## Set parameter values and plot the relationship
time <- seq(0, 150,5)
ans1 <- gompertz(time, tmax = 60, theta = 0.09, Yo=100)
ans2 <- gompertz(time, tmax = 60, theta = 0.07, Yo=100)
ans3 <- gompertz(time, tmax = 60, theta = 0.05, Yo=100)
ans4 <- gompertz(time, tmax = 60, theta = 0.03, Yo=100)
xyplot(ans1 + ans2 + ans3 + ans4 ~ time, type="l", auto.key=TRUE, ylab = "text", xlab = "time")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
