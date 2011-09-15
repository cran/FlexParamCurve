\name{posneg_data}
\docType{data}
\alias{posneg_data}
\title{Simulated growth of black-browed albatrosses}
\description{
  The \code{posneg_data} data frame has 1300 rows and 3 columns of records of
  the simulated masses for black-browed albatross chicks between 0 and 166 
  days of age.
}
\usage{data(posneg_data)}
\format{
  This object of class \code{c("nfnGroupedData", "nfGroupedData",
    "groupedData", "data.frame")} containing the following columns:
  \describe{
    \item{mass}{
      a numeric vector of chick masses (g).
    }
    \item{age}{
      a numeric vector of chick ages (days).
    }
    \item{id}{
      an ordered factor indicating unique id of each simulated individual,
      i.e. which data belongs to which individual.
    }
  }
}
\details{   
  No published parameter estimates with associated variability are available for
  positive-negative growth curves. These data were simulated using
  an 8-parameter positive-negative Richards curve (\code{\link{SSposnegRichards}}
  (model 1)), using parameters drawn from normal distributions with the following
  means (standard deviations):
  \preformatted{
  Asym=4300 (180)
  K=0.06 (0.01)
  Infl=23 (0.4)
  M=0.1 (0.05)
  RAsym=1433.3 (540) #1/3 of Asym, more variable
  Rk=0.108 (0.03) #1.8 times faster recession, more variable
  Ri=Infl+87.259 (1.7) # more variable but linked to Infl
  RM=M (0.15) #more variable}
  These values were chosen through comparison of growth curves with Huin and Prince (2000) Fig 2
  and variability observed between individual chicks of little penguins in a 10 year dataset
  (Chiaradia and Nisbet unpublished data). Each simulated individual had 13 measurements
  stratified through the development period, with 1-13 day random differences in timing 
  of each measurement.
  This data object has methods for \code{\link{nlme}} grouped-data classes. 
}
\source{
  Huin, N. & Prince, P.A. (2000) 
  Chick growth in albatrosses: curve fitting with a twist. 
  Journal of Avian Biology, 31, 418-425.
}
\examples{
require(stats); require(graphics)
data(posneg_data)
modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
plot(mass ~ age, data = posneg_data, subset = id == "1",
     xlab = "Chick age (day)", las = 1,
     ylab = "Chick mass (g)",
     main = "posneg_data data and fitted curve (Chick #1 only)")
fm1 <- nls(mass ~ SSposnegRichards(age,Asym=Asym,K=K,Infl=Infl,M=1,RAsym=RAsym,
           Rk=Rk,Ri=Ri,RM=1,modno=22),
           data = posneg_data, subset = id == "1")
age <- seq(0, 166, length.out = 101)
lines(age, predict(fm1, list(age = age)))
}
\keyword{datasets}
