\name{tern.data}
\docType{data}
\alias{tern.data}
\title{Field data on growth of common terns \eqn{Sterna hirundo}}
\description{
  The \code{tern.data} data frame has 1164 rows and 12 columns of records of
  the measured masses for common tern chicks between 1 and 30 
  days of age collected at at Grays Beach, MA, in 1973 (Nisbet 1975) and
  Monomoy, MA, in 1975 (Nisbet et al. 1978)
}
\usage{tern.data}
\format{
  This object of class \code{c("nfnGroupedData", "nfGroupedData",
    "groupedData", "data.frame")} containing the following columns:
  \describe{
    \item{site}{
      Four character factor for the two sites (MYMA: Monomoy Island, MA; GBCT: Grays Beach, CT).
    }
    \item{year}{
      A factor specifying the year of measurement: 1973 or 1976.
    }
    \item{bandid}{
      an ordered factor indicating unique id of each individual: the union of
      the laying date of the nest relative to the colony and the band combination
    }
    \item{siteyear}{
      A factor specifying levels of year for different sites (different years at each site).
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
      A factor of hatching order for each chick (A = first hatched, B = second hatched
      C = third hatched).
    }
    \item{outcome}{
      A factor of codes for fate of each chick (F = fledged; only fledged chicks included).
    }
    \item{eggmass}{
      A numeric vector of the mass of the egg (from which the chick hatched) at laying.
    }
    \item{clutch}{
      A factor of size of clutch/brood that each chick comes from (either 1- or 2-chick brood).
    }
  }
}
\details{   
  Data were collected as outlined in Nisbet (1975)[Grays Beach, MA, 1973] and
  Nisbet et al.(1978) [Monomoy, MA, 1975].
  Please contact Ian Nisbet <icnisbet@verizon.net> for use in collaborations.
}
\source{
  Nisbet, I.C.T. (1975) Selective effects of predation in a tern colony. Condor, 77, 221-226.
  Nisbet, I.C.T., Wilson, K.J. & Broad, W.A. (1978) Common Terns raise young after death
  of their mates. Condor, 80, 106-109.
}
\examples{
require(stats); require(graphics)
#view data
tern.data
#create pnmodelparams for fixed parameters
modpar(tern.data$ckage, tern.data$weight, force4par = TRUE, pn.options = "myoptions")
plot(weight ~ ckage, data = tern.data, subset = bandid == tern.data$bandid[1],
     xlab = "Chick age (day)", las = 1,
     ylab = "Chick mass (g)",
     main = "tern.data and fitted curve (Chick #156 only)")
fm1 <- nls(weight ~ SSposnegRichards(ckage,Asym=Asym,K=K,Infl=Infl,modno=32,
           pn.options= "myoptions"),
           data = tern.data, subset = bandid == tern.data$bandid[1])
ckage <- seq(0, 30, length.out = 101)
lines(ckage, predict(fm1, list(ckage = ckage)))
}
\keyword{datasets}
