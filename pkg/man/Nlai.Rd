\name{Nlai}
\alias{Nlai}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Leaf N distribution within canopies
  %% ~~function to do ... ~~
}

\description{
  Function used to describe
  leaf N distribution within canopies (Jonshon et al., 2010 AoB). This
  is equation 1.9 in Table S.1 in the paper referenced below.
}

\usage{
Nlai(lai, Yo, nb, kn, c)
}

\arguments{
  \item{lai}{
   Cummulative leaf area index. Starting from the top. Should be zero
  right above the canopy. 
}
\item{Yo}{
  Value of leaf N concentration at lai = 0
}
  \item{nb}{
  Minimum N for photosynthesis
}
\item{kn}{
  extinction coefficient
}
\item{c}{
 emprical paramater: when c=0 then the there is no decay (uniform distribution),
 when c=1 then the decay is exponential,
 when c>1 then the decay takes an inverse S shape

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
lai  <- seq(1, 8, 0.5)
ans1 <- Nlai(lai, Yo=2, nb=0.4, kn=0.5, c=1)
ans2 <- Nlai(lai, Yo=2, nb=0.4, kn=0.5, c=4)
ans3 <- Nlai(lai, Yo=2, nb=0.4, kn=0.5, c=8)
ans4 <- Nlai(lai, Yo=2, nb=0.4, kn=0.5, c=12)
xyplot(ans1 + ans2 + ans3 + ans4~ LAI, type="l", auto.key=TRUE, ylab = "N distribution", xlab = "LAI")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

