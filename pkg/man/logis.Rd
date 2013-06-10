\name{logis}
\alias{logis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic growth function
}
\description{
Teh logistic growth function corresponding to the SSlogis selfStart
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
  \code{\link{SSlogis}}
}
\examples{
logis(3, 20, 5, 1)
}
\keyword{ models }
