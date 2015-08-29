\name{posnegRichards.calls}
\docType{data}
\alias{posnegRichards.calls}
\title{List of calls for fitting 33 SSposnegRichards models in \code{\link{nlsList}}}
\description{
  The \code{posnegRichards.calls} list has two components of 17 and 16 rows and 1 column, 
  respectively, called 'Examples of calls for FlexParamCurve models that estimate parameter m'
  (models with 4 estimable first curve parameters) and "Examples of calls for FlexParamCurve models
  that fix parameter m" (models with 3 estimable second curve parameters, i.e. M is fixed
  to value in \eqn{pnmodelparams}. Individual calls can be accessed by indexing first the component
  number and then the model number - see examples below.
  Note that model 17 is formulated differently (see \code{\link{SSposnegRichards}})
}
\usage{posnegRichards.calls}
\format{
  This object of class \code{\link{list}} containing the components:
  \describe{
    \item{Examples of calls for FlexParamCurve models that estimate parameter m}{
      a list of 16 possible reductions (nos. 1-16) of the FlexParamCurve double-Richards model 
      that estimate parameter m. Also includes a custom model (17; see \code{\link{SSposnegRichards}}).
    }
    \item{Examples of calls for FlexParamCurve models that fix parameter m}{
      a list of 16 possible reductions (nos. 21-36) of the FlexParamCurve double-Richards model 
      that do not estimate parameter m but instead fix it to a mean across the dataset 
      or user-specified value.
    }    
  }
}
\details{   
  A list object to provide users with examples of how to fit 33 different \code{\link{nlsList}}
  models using the selfStart function \code{\link{SSposnegRichards}}.
}
\examples{
# see all possible calls
posnegRichards.calls
# extract the call for fitting a nls model with 8-parameter double-Richards curve (model 1)
#for an example just fit a subset of the data, 3 group levels (individuals)
   data <- subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
modtofit <- as.character( 
               posnegRichards.calls [[2]] [row.names(posnegRichards.calls [[2]]) == "22",] )
#change the data source
modtofit <- sub("posneg.data","data",modtofit)
modtofit <- parse(text = modtofit)
#create list for fixed parameters
modpar(posneg.data$age, posneg.data$mass, pn.options = "myoptions")
#create a new nlsList object called richards22.lis
eval(modtofit)
#view object
richardsR22.lis

# view call for model 1
posnegRichards.calls [[1]] [row.names(posnegRichards.calls [[1]]) == "1",]

# view call for model 21
posnegRichards.calls [[2]] [row.names(posnegRichards.calls [[2]]) == "21",]
}
\keyword{datasets}
