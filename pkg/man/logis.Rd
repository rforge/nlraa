\name{logis}
\alias{logis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic growth function
}
\description{
The logistic growth function corresponding to the SSlogis
selfStart. This is Eq 2.1 in Table 1 in the paper referenced below.
}
\usage{
logis(x, Asym, xmid, scal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector of values at which to evaluate the model.
}
\item{Asym}{
   a numeric parameter representing the asymptote.
}
  \item{xmid}{
 a numeric parameter representing the ‘x’ value at the
          inflection point of the curve.  The value of ‘SSlogis’ will
          be ‘Asym/2’ at ‘xmid’.
}
  \item{scal}{
 a numeric scale parameter on the ‘input’ axis.
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
  \code{\link{SSlogis}}
}
\examples{
logis(3, 20, 5, 1)
}
\keyword{ models }
