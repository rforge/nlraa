\name{sm}
\alias{sm}
\docType{data}
\title{
 Biomass accumulation for Maize, Fiber Sorghum and Sweet Sorghum collected by Sotiris Arhcountoulis
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(sm)}
\format{
  A data frame with 236 observations on the following 5 variables.
  \describe{
    \item{\code{DOY}}{a numeric vector, day of the year}
    \item{\code{Block}}{a numeric vector, block structure of the design}
    \item{\code{Input}}{a numeric vector, level of agricultural inputs. 1 is low and 2 is high}
    \item{\code{Crop}}{a factor with levels \code{F} \code{M} \code{S}. F is for fiber sorghum, M is for Maize and S is for sweet sorghum}
    \item{\code{Yield}}{a numeric vector, biomass yield in Mg/ha}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(sm)
## See vignette for examples
}
\keyword{datasets}
