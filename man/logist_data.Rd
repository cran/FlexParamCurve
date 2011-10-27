\name{logist_data}
\docType{data}
\alias{logist_data}
\title{Simulated growth of whiskered terns}
\description{
  The \code{logist_data} data frame has 1100 rows and 3 columns of records of
  the simulated masses for whiskered tern chicks between 0 and 21 
  days of age.
}
\usage{data(logist_data)}
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
  an 3-parameter positive-negative Richards curve (\code{\link{SSposnegRichards}}
  (model 20)), using parameters drawn from normal distributions with the following
  means (standard deviations):
  \preformatted{
  Asym=92.35 (15.65)
  K=0.06 (0.138)
  Infl=0.294 (1.72)}
  These values were taken from Pallisson et al. (2008) for 75 chicks reported.
  Each simulated individual had 11 measurements
  stratified through the development period, with 1-2 day random differences in timing
  of each measurement.
  This data object has methods for \code{\link{nlme}} grouped-data classes. 
}
\source{
  Paillisson, J.-M., Latraube, F. & Reeber, S. (2008) 
  Assessing growth and age of Whiskered Tern 
  \eqn{Chlidonias hybrida} chicks using biometrics. Ardea, 96, 271-277.
}
\examples{
require(stats); require(graphics)
data(logist_data)
modpar(logist_data$age, logist_data$mass, force4par = TRUE) #create pnmodelparams for fixed parameters
plot(mass ~ age, data = logist_data, subset = id == "0.002",
     xlab = "Chick age (day)", las = 1,
     ylab = "Chick mass (g)",
     main = "logist_data and fitted curve (Chick #2 only)")
change.pnparameters(M=1) # set curve to logistic (M=1) in subsequent fit
fm1 <- nls(mass ~ SSposnegRichards(age,Asym=Asym,K=K,Infl=Infl,
           modno=20),
           data = logist_data, subset = id == "0.002")
age <- seq(0, 166, length.out = 101)
lines(age, predict(fm1, list(age = age)))
}
\keyword{datasets}
