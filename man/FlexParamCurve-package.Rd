\name{FlexParamCurve-package}
\alias{FlexParamCurve-package}
\docType{package}
\title{Tools to Fit Flexible Parametric Curves}
\description{selfStart functions and model selection tools to fit parametric curves in 
  \code{\link{nls}}, \code{\link{nlsList}} and \code{\link{nlme}} frameworks.}
\details{
    User level functions include:
    
    \eqn{pn_mod_compare}
    \preformatted{
      all-model selection for positive-negative Richards nlsList models
        }
    \eqn{pn_modselect_step}
    \preformatted{
      backward stepwise model selection for positive-negative Richards nlsList models
        }
    \eqn{SSposnegRichards}
    \preformatted{
      selfStart function for estimating parameters of 36 possible reductions of the 8-parameter
      positive-negative Richards model
        }
    \eqn{posnegRichards_eqn}
    \preformatted{
      function for evaluating 36 possible reductions of the 8-parameter 
      positive-negative Richards model
        }
    \eqn{modpar}
    \preformatted{
      estimates mean parameters (and parameter bounds) for 8-parameter positive-negative Richards
      models or 4-parameter Richards models and saves in objects pnmodelparams
      and pnmodelparamsbounds. (required prior to use of the above functions)
        }
    \eqn{change.pnparameters}
    \preformatted{      
      simple function to update pnmodelparams and pnmodelparamsbounds
      with user specified values
        }
    \eqn{extraF}
    \preformatted{
      performs extra sum-of-squares F test for two nested nlsList models
        }
    \eqn{extaF.nls}
    \preformatted{
      performs extra sum-of-squares F test for two nested nls models
        }
 \tabular{ll}{Package: \tab FlexParamCurve\cr
Type: \tab Package\cr
Title: \tab Tools to Fit Flexible Parametric Curves\cr
Version: \tab 1.0\cr
Date: \tab 2011-09-15\cr
Author: \tab Stephen Oswald\cr
Maintainer: \tab Stephen Oswald <ozsao23@hotmail.com>\cr
License: \tab GPL-2\cr
Depends: \tab nlme\cr
Enhances: \tab nlme\cr
LazyLoad: \tab yes\cr}
}
\examples{
# run all-model selection for posneg_data object without need to run any previous functions
\donttest{
data(posneg_data)
modseltable <- pn_mod_compare(posneg_data$age, posneg_data$mass,
    posneg_data$id, existing = FALSE)}

# run backwards stepwise model selection for logist_data object
\donttest{
data(logist_data)
modseltable <- pn_modselect_step(logist_data$age, logist_data$mass,
    logist_data$id, existing = FALSE)}

# estimate fixed parameters use data object posneg_data
\donttest{
modpar(posneg_data$age,posneg_data$mass)}

# change fixed values of M and constrain hatching mass to 45.5 in a growth curve
\donttest{
change.pnparameters(M=1,RM=0.5,first_y=45.5)}
    
# fit nlsList object using 6 parameter model with values M and RM
# fixed to value in pnmodelparams and then fit nlme model
\donttest{
richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, M = 1, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = 1,
      modno = 22), data = posneg_data)
 richardsR22.nlme <- nlme(richardsR22.lis)}
 
# fit reduced nlsList model and then compare performance with extraF
\donttest{
richardsR30.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, M = 1, RAsym = 1, Rk = 1, Ri = Ri, RM = 1,
      modno = 30), data = posneg_data)
 extraF(richardsR30.lis,richardsR22.lis)}
}
\seealso{\code{\link{nlme}}
\code{\link{SSlogis}}
\code{\link[richards]{richards}}
}
\author{Stephen Oswald <ozsao23@hotmail.com>}

\keyword{
Curve fit
Parametric curves
Growth
nlme
nlsList
double logistic
logistic equation
positive negative curve
}


