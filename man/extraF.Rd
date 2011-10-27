\name{extraF}
\alias{extraF}
\title{Compare Two \eqn{nlsList} Models Using Extra Sum-of-Squares F-Tests}

\usage{
extraF(submodel = 1,
genmodel = 1)
}
\arguments{
  \item{submodel}{\eqn{nlsList} model with fewer curve parameters (reduced model)}
  \item{genmodel}{\eqn{nlsList} model with more curve parameters (general model)}
}
\description{Function to compare two nested models using extra sum-of-squares F-Tests.}
\details{
Models must be entered in the correct order with the reduced model appearing
first in the call and the more general model appearing later. These must be nested models,
i.e. the general model must contain all of the curve parameters in the reduced model and more.
Entering models with the same number of parameters will produce NAs in the output, but
the function will produce seemingly adequate output with non-nested models. The user must
check that models are nested prior to use.

This function is primarily designed to be called by the model selection functions
\code{\link{pn_modselect_step}} and \code{\link{pn_mod_compare}} but can be used independently.

Extra sum-of-squares is obtained from: \preformatted{F = (SS1 - SS2)(df1 - df2) / (SS2 / df2)}
where SS = sum-of-squares and df = degrees of freedom, for the more reduced model (1) and the
more general model (2), respectively.

If the F value is significant then the more general model provides a significant improvement
over the reduced model, but if the models are not significantly different then the reduced
parameter model is to be preferred.}
\value{A \code{\link{data.frame}} listing the names of the models compared, F,
numerator degrees of freedom,
demonimator degrees of freedom, P value and the residual sum of squares for both the general
and reduced models}
\references{Ritz, C. and Streibigg, J. C. (2008) \eqn{Nonlinear regression with R.}
Springer-Verlag, New York.}
\author{Stephen Oswald <ozsao23@hotmail.com>}



\seealso{\code{\link{extraF.nls}}
\code{\link{nlsList}}
\code{\link{pn_modselect_step}}
\code{\link{pn_mod_compare}}}
\examples{
   #compare two nested nlsList models (4 vs 8 parameter models)
\donttest{
   data(posneg_data)
   modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
   # (only first 4 group levels in data used for example's sake)
   subdata<-subset(posneg_data, as.numeric(row.names (posneg_data) ) < 53)
   richardsR2.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
     Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri,modno = 2)
                        , data = subdata)
   richardsR12.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
     Infl = Infl, M = M, modno = 12)
                        , data = subdata)
   extraF(richardsR12.lis, richardsR1.lis)}
}
