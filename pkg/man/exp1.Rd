\name{exp1}
\alias{exp1}
\alias{exp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simple exponential functions
%%  ~~function to do ... ~~
}
\description{
  First order exponential decrease (exp1) or increase (exp2). These are
  equations 1.1 and 1.2 in Table 1 in the paper referenced below.
}
\usage{
exp1(time, Yo, rate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
    a numeric vector of values at which to evaluate the model.
}
\item{Yo}{
     a numeric parameter representing the asymptote.
}
\item{rate}{
  rate of decrease
%%     ~~Describe \code{rate} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(lattice)
## Set parameter values and plot the relationship
time <- seq(0, 100,5)
ans1 <- exp1(time, Yo=1, rate=0.01)
ans2 <- exp1(time, Yo=1, rate=0.05)
ans3 <- exp1(time, Yo=1, rate=-0.01)

xyplot(ans1 + ans2 + ans3 ~ time, type="l", auto.key=TRUE,
       ylab = "OM decomposition", xlab = "time")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
