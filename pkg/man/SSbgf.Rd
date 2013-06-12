\name{SSbgf}
\alias{SSbgf}
\alias{bgf}
\alias{bgf2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
selfStart for the beta growth function
}
\description{
Implements some self starting rules for the beta growth function as
desribed by Yin et. al (2003). This is equation 2.5 in the paper
referenced below (Archontoulis and Miguez).
}
\usage{
SSbgf(time, w.max, t.e, t.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
    a numeric vector of values at which to evaluate the model. This is
    typically expressed in some measure of time.
    
  }
  \item{w.max}{ A numeric parameter representing the maximum weight (or
    biomass).  }

  \item{t.e}{ A numeric parameter representing the time at
    which maximum weight is attained (w.max).  }

  \item{t.m}{ A numeric
parameter representing the time at which the inflection point of the
curve is attained. The maximum growth rate ocurrs at this time.  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a numeric vector of the same length as the input
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Yin, X., J. Goudriaan, E.A. Lantinga, J. Vos, and J.H.J. Spiertz. 2003b. A flexible sigmoid function of determinate growth. Annals of Botany 91: 361–371. (with
erratum in Annals of Botany 91: 753, 2003).

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
}
\examples{
## See the vignette for an example
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
