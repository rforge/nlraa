\name{rich}
\alias{rich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Richards equation
}
\description{
  Richards equation (Richards, 1959) J Exp Bot 10: 290-300. This is
  equation 2.2 in Table 1 in the paper referenced below (Archontoulis
  and Miguez).
}
\usage{
rich(time, tmax, theta, v, Yo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
a numeric vector of values at which to evaluate the model.
}
\item{tmax}{
  tmax is the time when growth rate is maximum
}

\item{theta}{
 thetha determines the curvature. The relative growth rate is 
  RGR = theta/(1+v).
}
\item{v}{
  v parameter determines the asymetrical growth. 
%%     ~~Describe \code{v} here~~
}
\item{Yo}{
  Yo is the max Y
}
}
\details{
Richards equation is an extension of the logistic to deal with asymmetrical growth.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  a numeric vector of length equal to the inputs
}
\references{
  Richards equation (Richards, 1959) J Exp Bot 10: 290-300
  
  Nonlinear Regression Models and Applications in Agricultural Research.
  Sotirios V. Archontoulis and
  Fernando E. Miguez. Agronomy Journal. doi: 10.2134/agronj2012.0506 
}
\author{
  Fernando Miguez
%%  ~~who you are~~
}
\note{
Applied to describe water stress index, crop N-uptake, seed germination, crop growth, LAI development, etc.
Remark at time = 0, Y is not zero 
Richards equation has often been criticized because the shape parameter v has no obvious biological interpretation
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
There used to be an R package called richards but it is no longer maintained.
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(lattice)
time <- seq(0, 150,5)
ans1 <- rich(time, tmax = 60, theta = 0.07, v=1.3, Yo=100)
ans2 <- rich(time, tmax = 60, theta = 0.07, v=0.6, Yo=100)
ans3 <- rich(time, tmax = 60, theta = 0.07, v=0.3, Yo=100)
ans4 <- rich(time, tmax = 60, theta = 0.07, v=1.0, Yo=100)
xyplot(ans1 + ans2 + ans3 + ans4 ~ time,
       type="l", auto.key=TRUE, ylab = "text", xlab = "time")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
