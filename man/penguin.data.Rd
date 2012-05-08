\name{penguin.data}
\docType{data}
\alias{penguin.data}
\title{Field data on growth of little penguins \eqn{Eudyptula minor}}
\description{
  The \code{penguin.data} data frame has 2244 rows and 11 columns of records of
  the measured masses for little penguin chicks between 13 and 74 
  days of age collected at Philip Island, Victoria, Australia in 2000 and 2002
  (see Chiaradia & Nisbet 2006).
}
\usage{data(penguin.data)}
\format{
  This object of class \code{c("nfnGroupedData", "nfGroupedData",
    "groupedData", "data.frame")} containing the following columns:
  \describe{
    \item{site}{
      Three character factor for the site (only one site in dataset).
    }
    \item{year}{
      A factor specifying the year of measurement.
    }
    \item{bandid}{
      an ordered factor indicating unique id of each individual: the union of
      the laying date of the nest relative to the colony and the band combination
    }
    \item{siteyear}{
      A factor specifying levels of year for different sites (only one site in dataset).
    }
    \item{weight}{
      a numeric vector of chick masses (g).
    }
    \item{ckage}{
      a numeric vector of chick ages (days).
    }
    \item{Jdate}{
      a numeric vector of first egg-laying date of the nest(days), relative to the mean
      laying date for all nests in that year.
    }
    \item{nest}{
      A factor of unique codes that identify each nest.
    }
    \item{ck}{
      A factor of hatching order for each chick (A = first hatched, B = second hatched).
    }
    \item{outcome}{
      A factor of codes for fate of each chick (F = fledged; only fledged chicks included).
    }
    \item{clutch}{
      A factor of size of clutch/brood that each chick comes from (either 1- or 2-chick brood).
    }
  }
}
\details{   
  Data were collected as outlined in Chiaradia & Nisbet (2006). Penguin chicks are generally considered
  to exhibit a double-Gompertz growth form.Please contact Andre Chiaradia (a.chiaradia@penguins.org.au) 
  for use in collaborations.
}
\source{
  Chiaradia, A. & Nisbet, I.C.T. (2006) 
  Plasticity in parental provisioning and chick growth in Little Penguins \eqn{Eudyptula minor}
  in years of high and low breeding success. Ardea, 94, 257-270. 
}
\examples{
require(stats); require(graphics)
data(penguin.data)
\donttest{
modpar(penguin.data$ckage, penguin.data$weight, pn.options = "myoptions")
plot(weight ~ ckage, data = penguin.data, subset = bandid == penguin.data$bandid[1],
     xlab = "Chick age (day)", las = 1,
     ylab = "Chick mass (g)",
     main = "penguin.data and fitted curve (Chick #307 only)")
fm1 <- nls(weight ~ SSposnegRichards(ckage,Asym=Asym,K=K,Infl=Infl, RAsym=RAsym,modno=31, pn.options= "myoptions"),
           data = penguin.data, subset = bandid == penguin.data$bandid[1])
ckage <- seq(0, 74, length.out = 101)
lines(ckage, predict(fm1, list(ckage = ckage)))
}
}
\keyword{datasets}
