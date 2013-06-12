\name{weibull}
\alias{weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Weibull equation (Weibull, 1951) 

}
\description{
   Weibull equation (Weibull, 1951). This is equation 2.4 in Table 1 of
   the paper referenced below.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weibull(time, a, b, Yo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
    a numeric vector of values at which to evaluate the model.
}
  \item{a}{
 a and b paramaters define the shape of the response
}
\item{b}{
   a and b paramaters define the shape of the response
}
\item{Yo}{
  Yo is the max value for Y and when time is zero, Y is also zero
}
}
\details{
Applied to describe water stress index, crop N-uptake, seed germination, crop growth, LAI development, etc.
Remark: a and b paramaters have not direct biological interperation
and the unit of paramater a depends on paramater b and therefore is
difficult to provide initial estimates for non-linear regression
analysis. Usually the initial paramaters were provided by trial and
error. But see SSweibull.
}
\value{
  a numeric vector of length equal to the inputs

}
\references{
  Nonlinear Regression Models and Applications in Agricultural Research.
  Sotirios V. Archontoulis and
  Fernando E. Miguez. Agronomy Journal. doi: 10.2134/agronj2012.0506 
}
\author{
  Fernando E. Miguez
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{SSweibull}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Set parameter values and plot the relationship
time <- seq(0, 200,5)
ans1 <- weibull(time, a = 0.00025, b = 0.5, Yo=100)
ans2 <- weibull(time, a = 0.00025, b = 1.5, Yo=100)
ans3 <- weibull(time, a = 0.00025, b = 2, Yo=100)
ans4 <- weibull(time, a = 0.00025, b = 3, Yo=100)
xyplot(ans1  + ans2 + ans3 + ans4 ~ time, type="l", auto.key=TRUE, ylab = "text", xlab = "time")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
