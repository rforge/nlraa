\name{exp3}
\alias{exp3}
\alias{exp4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Parallel or double first order exponential
}
\description{
  Parallel or double first order exponential
}
\usage{
exp3(time, Yo, rate1, rate2, fraction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
        a numeric vector of values at which to evaluate the model.
}
  \item{Yo}{
     a numeric parameter representing the asymptote.
}
\item{rate1}{
  first rate
}
\item{rate2}{
  second rate
}
\item{fraction}{
  Fraction (0-1) which determines the size of each pool
%%     ~~Describe \code{fraction} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    a numeric vector of length equal to the inputs

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Fernando E. Miguez
  Sotiris Archontoulis
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
## Set parameter values and plot the relationship
time <- seq(0, 100,5)
ans1 <- exp3(time, Yo=1, rate1=0.01, rate2=0.05, fraction=0.3)
ans2 <- exp3(time, Yo=1, rate1=0.01, rate2=0.05, fraction=0.1)
ans3 <- exp3(time, Yo=1, rate1=0.01, rate2=0.05, fraction=0.9)
xyplot(ans1 + ans2 + ans3 ~ time, type="l", auto.key=TRUE, ylab = "text", xlab = "time")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
