\name{FlexParamCurve-package}
\alias{FlexParamCurve-package}
\docType{package}
\title{Tools to Fit Flexible Parametric Curves}
\description{selfStart functions and model selection tools to fit parametric curves in 
  \code{\link{nls}}, \code{\link{nlsList}} and \code{\link{nlme}} frameworks.}
\details{

    General approach for using package (also see examples below)
    1) Run \eqn{modpar} to produce initial parameter estimates and estimates of parameter bounds for your dataset. 
These are used to accomodate fixed parameters and are saved in user-specified \code{\link{list}} object
All parameters and options in this list can be edited manually or using \eqn{change.pnparameters}. The
list could be created manually given that the elements were labelled sufficiently.

    2) Determine most appropriate model (number of necessary parameters) for your data
using  \eqn{pn.mod.compare} or \eqn{pn.modselect.step} (these rank competing model and then compare nested models using
\eqn{extraF}). This may take some time as many \eqn{nlsList} objects are fitted.
Note that if you perform this step, then you do not need to perform step 1.
If you are sure of your model (e.g. it is a simple logistic) Step 2 may be unnecessary.

    3) Fit \eqn{nls} or \eqn{nlsList} or \eqn{nlme} models using \eqn{SSposnegRichards} specifying
the appropriate model number and the list of parameters and options. Parameter bounds can be refined to improve fits by altering 
this list, either manually or using \code{\link{change.pnparameters}}.

    4) Plot your curves using eqn{posnegRichards.eqn} specifying the appropriate model number and list of parameters/options.


    User level functions include:
    
    \eqn{pn.mod.compare}
    \preformatted{
      all-model selection for positive-negative Richards nlsList models
        }
    \eqn{pn.modselect.step}
    \preformatted{
      backward stepwise model selection for positive-negative Richards nlsList models
        }
    \eqn{SSposnegRichards}
    \preformatted{
      selfStart function for estimating parameters of 36 possible reductions of the 8-parameter
      positive-negative Richards model (double-Richards)
        }
    \eqn{posnegRichards.eqn}
    \preformatted{
      function for evaluating 36 possible reductions of the 8-parameter 
      positive-negative Richards model (double-Richards)
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
Version: \tab 1.4-3\cr
Date: \tab 2012-09-23\cr
Author: \tab Stephen Oswald\cr
Maintainer: \tab Stephen Oswald <ozsao23@hotmail.com>\cr
License: \tab GPL-2\cr
Depends: \tab nlme\cr
Enhances: \tab nlme\cr
LazyLoad: \tab yes\cr}
}
\references{## Oswald, S.A. et al. (Early view) FlexParamCurve: R package for flexible
     fitting of nonlinear parametric curves. Methods in Ecology and Evolution.
     doi: 10.1111/j.2041-210X.2012.00231.x (see also tutorial and introductory videos at:
     http://www.methodsinecologyandevolution.org/view/0/podcasts.html
     (posted September 2012 - if no longer at this link, check the archived videos at:
     http://www.methodsinecologyandevolution.org/view/0/VideoPodcastArchive.html#allcontent
}     
\examples{
# run all-model selection for posneg.data object (Step 2) without need to run any previous functions
\donttest{
data(posneg.data)
modseltable <- pn.mod.compare(posneg.data$age, posneg.data$mass,
    posneg.data$id, existing = FALSE, pn.options = "myoptions")}

# run backwards stepwise model selection (Step 2) for logist.data object without need to run any previous functions
\donttest{
data(logist.data)
modseltable <- pn.modselect.step(logist.data$age, logist.data$mass,
    logist.data$id, existing = FALSE, pn.options = "myoptions")}

# estimate fixed parameters use data object posneg.data (Step 1)
\donttest{
modpar(posneg.data$age,posneg.data$mass, pn.options = "myoptions")}

# change fixed values of M and constrain hatching mass to 45.5 in a growth curve (Step 1)
\donttest{
change.pnparameters(M=1,RM=0.5,first.y=45.5, pn.options = "myoptions")}
    
# fit nlsList object using 6 parameter model with values M and RM (Step 3)
# fixed to value in pnmodelparams and then fit nlme model
\donttest{
richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, RAsym = RAsym, Rk = Rk, Ri = Ri,
      modno = 22, pn.options = "myoptions"), data = posneg.data)
 richardsR22.nlme <- nlme(richardsR22.lis, random = pdDiag(Asym + Infl ~ 1))}
 
# fit reduced nlsList model and then compare performance with extraF (manual version of Step 2)
\donttest{
richardsR20.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, modno = 20, pn.options = "myoptions"), data = posneg.data)
 extraF(richardsR20.lis,richardsR22.lis)}
 
 # fit and plot a logistic curve (M=1) to data, note that all parameters set to equal 1 here are ignored
 # note code here forces \eqn{modpar} to only estimate 4 curve parameters (simple Richards curve)
\donttest{
data(logist.data)
modpar(logist.data$age,logist.data$mass,force4par=TRUE, pn.options = "myoptions")#create list for fixed parameters
change.pnparameters(M=1) # set M to 1 for subsequent fit
richardsR20.nls <- nls(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, modno = 20, pn.options = "myoptions"), data = logist.data)
plot(logist.data$age , logist.data$mass, xlab = "age", ylab = "mass", pch = ".", cex = 0.7) 
par <- coef( richardsR20.nls )
#(Step 4)
curve(posnegRichards.eqn(x, Asym = par[1], K = par[2], Infl = par[3], modno = 20), add= TRUE, pn.options = "myoptions")}
}
\seealso{\code{\link{nlme}}
\code{\link{SSlogis}}
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


