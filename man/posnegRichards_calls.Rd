\name{posnegRichards_calls}
\docType{data}
\alias{posnegRichards_calls}
\title{List of calls for fitting 36 SSposnegRichards models in \code{\link{nlsList}}}
\description{
  The \code{posnegRichards_calls} list has two components of 20 and 17 rows and 1 column, 
  respectively, called 'Variable M models' (models with 4 estimable positive parameters) 
  and "Fixed M models" (models with 3 estimable positive parameters, i.e. M is fixed
  to value in \eqn{pnmodelparams}.
}
\usage{data(posnegRichards_calls)}
\format{
  This object of class \code{\link{list}} containing the components:
  \describe{
    \item{Examples of 36 calls}{
      a list of 36 possible FlexParamCurve models, including a custom model (17)
      and two models used by internally by \code{\link{modpar}} (18 and 19)
    }
  }
}
\details{   
  A list object to provide users with examples of how to fit 36 different \code{\link{nlsList}}
  models using the selfStart function \code{\link{SSposnegRichards}}.
}
\examples{
# extract the call for fitting a nls model with 8-parameter positive-begative Richards curve (model 1)
data(posnegRichards_calls)
data(posneg_data)
modtofit <- as.character( posnegRichards_calls [[1]] [1,] )
modtofit <- parse(text = modtofit)
\donttest{
modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
eval(modtofit)}

# view call for model 20
posnegRichards_calls [[1]] [20,]

}
\keyword{datasets}
