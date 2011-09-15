\name{SSposnegRichards}
\alias{SSposnegRichards}
\title{Self-Starting Positive-Negative Richards Model}

\usage{
SSposnegRichards(x,
Asym,
K,
Infl,
M,
RAsym,
Rk,
Ri,
RM,
modno)
}
\arguments{
  \item{x}{a numeric vector of the primary predictor
variable at which to evaluate the model}
  \item{Asym}{a numeric value for the asymptote of the
positive (increasing) curve}
  \item{K}{a numeric value for the rate parameter of
the positive (increasing) curve}
  \item{Infl}{a numeric value for the point of inflection
of the positive (increasing) curve}
  \item{M}{a numeric value for the shape parameter of
the positive (increasing) curve}
  \item{RAsym}{a numeric value for  the asymptote of the
negative (decreasing) curve}
  \item{Rk}{a numeric value for the rate parameter of
the negative (decreasing) curve}
  \item{Ri}{a numeric value for the point of inflection
of the negative (decreasing) curve}
  \item{RM}{a numeric value for the shape parameter of
the negative (decreasing) curve}
  \item{modno}{a numeric value (currently integer only)
between 1 and 36 specifying the identification
number of the equation to be fitted}
}
\description{This selfStart function evaluates a range of flexible logistic
  functions and their gradients.  It also has an initial attribute that creates
  initial estimates of the parameters
  for the model specified.}
\details{This selfStart function evaluates a range of flexible logistic
functions and their gradients.  It also has an initial attribute that creates
initial estimates of the parameters
for the model specified. All 36 possible equations are biomodal,
in that they have a positive (increasing) trajectory followed by a negative
(decreasing) trajectory. These equations have also been described as
double-Richards curves, or positive-negative Richards curves.

The 36 possible equations are all based on the subtraction of one Richards
curve from another, producing:
\eqn{y = A / ([1+ m exp(-k (t-i))]1/m) - A' / ([1+ m' exp(-k' (t-i' ))]1/m' )},
where A=Asym, k=K, i=Infl, m=M,
A'=RAsym, k'=Rk, i'=Ri, m'=RM; as described in the Arguments section above.

All 36 possible equations are simply reformulations of this equation, in each
case fixing a parameter or multiple parameters to the mean parameter across
all individuals in the dataset (such as produced by a \code{\link{nls}}
model). Thus, a model in which one parameter is fixed has a 7-parameter equation,
and one in which four are fixed has a 4-parameter equation, thus reducing 
complexity and computation when a parameter does not vary across groups (e.g individuals)
[the most appropriate equation can be determined using model selection in
\code{\link{pn_modselect_step}} or \code{\link{pn_mod_compare}}].
Any models that require parameter fixing (i.e. all except #1)
extract appropriate values from the object \eqn{pnmodelparams} for the fixed
parameters. This object is created by running
\code{\link{modpar}} and can be adjusted manually or by using
\code{\link{change.pnparameters}} to user required specification.

Each of the 36 equations is identified by an integer value for modno (1 to 36).
All equations (except 17 - see below) contain parameters Asym, K, and Infl.
The list below summarizes which of the other 5 parameters are contained in
which of the models (Y indicates that the parameter is estimated, blank indicates
it is fixed).
\preformatted{
modno   M   RAsym   Rk   Ri   RM   NOTES
1       Y     Y     Y    Y    Y    8 parameter model
2       Y     Y     Y    Y
3       Y                Y    Y
4       Y     Y               Y
5       Y                     Y
6       Y     Y          Y    Y
7       Y     Y
8       Y     Y     Y         Y
9       Y           Y         Y
10      Y                Y
11      Y     Y
12      Y
13      Y     Y          Y
14      Y           Y    Y
15      Y     Y     Y
16      Y           Y
17      see below
18      see below
19      see below
20                                 only 3 parameters
21            Y     Y    Y    Y    7 parameter model, 4 recession params
22            Y     Y    Y
23                       Y    Y
24            Y               Y
25                            Y
26            Y          Y    Y
27            Y
28            Y     Y         Y
29                  Y         Y
30                       Y
31            Y
32
33            Y          Y
34                  Y    Y
35            Y     Y
36                  Y           }

modno 17 represents a different parameterization for a custom model:
(Asym/ 1 + exp(Infl - x)/ M) - (RAsym / 1 + exp(Ri - x)/ RM), in
which M and RM actually represent scale parameters not shape parameters.

modnos 18 and 19 are identical to 1 and 12, respectively, but are used
in the initial estimation of \eqn{pnmodelparams} and are not intended for
use in themselves.

The call for \code{\link{SSposnegRichards}} is slightly different from
conventional selfStart models in that it requires a value for modno
and any fixed parameters to be set equal to 1. Here is an example:
\preformatted{richardsR2.lis <- nlsList(mass ~ SSposnegRichards(age,
Asym = Asym, K = K, Infl = Infl, M = M, RAsy m= RAsym, Rk = Rk, Ri = Ri,
RM = 1, modno = 2), data = posneg_data)}

Examples for all models can be found in the list object
\code{\link{posnegRichards_calls}}.

If specified using \code{\link{modpar}} optional constraints
may be placed to specify response values at the minimum value and/or
maximum values of the predictor. Such constraint allows realistic fits for
datasets which are missing data
at either end of the curve (e.g. hatching weight for some growth curves).

Estimates are produced by splitting the two curves into separate positive
and negative curves and estimating parameters for each curve separately
in a similar manner to \code{\link{SSlogis}}. Each curve is fit first by
\code{\link{optim}} using the parameter bounds in \eqn{pnmodelparamsbounds}
(see \code{\link{modpar}}) and a subsequent refinement is attempted using
\code{\link{nls}} with more restrictive parameter bounds. Finally, both curves
are annealed
and parameters are again estimated using restrictive bounds and starting values
already determined during separate estimates. Equations for which the positive
curve was inestimable are not estimated further, but if negative curve estimation
or overall curve estimation fail, partial estimates are used: either default
negative parameters (RAsym = 0.05*Asym, Rk = K, Ri = Infl, RM = M) anealed to positive
curves or separate estimates annealed; both with compensation for interation
between asymptotes.

}
\value{a numeric vector of the same length as x containing parameter
estimates for equation specified (fixed parameters
are not return but are substituted in calls to \code{\link{nls}}
\code{\link{nlsList}} and \code{\link{nlme}} with the fixed parameters
stored in \eqn{pnmodelparams}; see \code{\link{modpar}}}
\references{#1# Nelder, J.A. (1962) Note: an alternative form of a generalized
logistic equation. Biometrics, 18, 614-616.
#2#
Huin, N. & Prince, P.A. (2000) Chick growth in albatrosses: curve fitting
with a twist. Journal of Avian Biology, 31, 418-425.
#3#
Pinheiro, J. & Bates, D. (2000) \eqn{Mixed-Effects Models in S and S-Plus.}
Springer Verlag, Berlin.}
\author{Stephen Oswald <ozsao23@hotmail.com>}
\note{Any models that require parameter fixing (i.e. all except #1)
extract appropriate values from the object \eqn{pnmodelparams} for the fixed
parameters. This object is created by running
\code{\link{modpar}} and can be adjusted manually or by using
\code{\link{change.pnparameters}} to user required specification.

Output may have many errors and warnings especially during a \code{\link{nlsList}}
fit, in which the function is called repeatedly: once for each group level in the
dataset. Warnings indicate conditions for which default parameters or incomplete estimates
are used - see Details section - and errors occur from insufficient data or singularities.
As a result of possible interaction and correlation between the parameters in some models,
singularities may be common, but do not be alarmed by repeated error messages, as
examination of a fitted \code{\link{nlsList}} model may releave a large number of
well estimated group levels, thus the elimation of unsuitable outlying groups only. Also,
because very few of the 36 equations are likely to be suitable for the majority
of datasets,
consideration of the model being fitted is crucial when examining the output. Functions
\code{\link{pn_modselect_step}} and \code{\link{pn_mod_compare}} provide the ability
for model selection of these equations through stepwise backward deletion or all
model comparison, respectively. These potentially offer powerful ways to determine the
best equation for your dataset.

To increase the ability of optimization routines to deal with
a wide variety of values, particularly negative values for M or RM,
only real component of complex numbers are modelled and integer versions
of M and RM are used during estimation if floating values cause conversion
issues - see \code{\link[richards]{richards}} for similar discussion.
}


\seealso{\code{\link{SSlogis}}
\code{\link{SSgompertz}}
\code{\link[richards]{SSrichards}}
\code{\link{posnegRichards_eqn}}}
\examples{
 require(graphics)
# retrieve mean estimates of 8 parameters using getInitial
# and posneg_data object
 data(posneg_data)
 modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
 getInitial(mass ~ SSposnegRichards(age, Asym, K, Infl, M, RAsym,
                                 Rk, Ri, RM, modno = 1), data = posneg_data)
                                 
# retrieve mean estimates and produce plot to illustrate fit
 pars<-as.numeric( getInitial(mass ~ SSposnegRichards(age, Asym, K, Infl,
       M, RAsym, Rk, Ri, RM, modno = 22), data = posneg_data) )
 plot(posneg_data$age, posneg_data$mass)
 curve( posnegRichards_eqn(x, Asym = pars[1], K = pars[2], Infl = pars[3], 
        M = 1, RAsym = pars[4], Rk = pars[5], Ri = pars[6], RM = 1, 
        modno = 2), lty = 3, xlim = c(0, 200) , add = TRUE)

# retrieve mean estimates and produce plot to illustrate fit for custom model 17
\donttest{
\dontrun{    pars<-as.numeric( getInitial(mass ~ SSposnegRichards(age, Asym, K, Infl,
       M, RAsym, Rk, Ri, RM, modno = 17), data = datansd) )
 plot(datansd$jday21March, datansd$moosensd)
 curve( posnegRichards_eqn(x, Asym = pars[1], K = 1, Infl = pars[2], 
        M = pars[3], RAsym = pars[4], Rk = 1, Ri = pars[5], RM = pars[6], 
        modno = 17), lty = 3, xlim = c(0, 200) , add = TRUE)}}
    

# fit nls object using 8 parameter model
# note: ensure data object is a groupedData object
 richardsR1.nls <- nls(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1)
                        , data = posneg_data)

# fit nlsList object using 8 parameter model
# note: ensure data object is a groupedData object
 subdata<-subset(posneg_data, as.numeric(row.names (posneg_data) ) < 53)
 richardsR1.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1)
                        , data = subdata)

# fit nlsList object using 7 parameter model, with RM fixed
# to value in pnmodelparams
 richardsR2.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = 1 , modno = 2)
                        ,data = subdata)

# fit nlsList object using 7 parameter model with value M
# fixed to value in pnmodelparams and then fit nlme model
\donttest{
 richardsR21.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, M = 1, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM,
      modno = 21), data = subdata)
 richardsR21.nlme <- nlme(richardsR21.lis)}

# estimate nlme model directly (not recommended: use
# nlme.nlsList instead to fit from nlsList object - see above example)
\donttest{
 richards22.nlme <- nlme(mass ~ SSposnegRichards(age, Asym = Asym,K = K,
     Infl = Infl, M = 1, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = 1, modno = 22),
            data = subdata,
            fixed = Asym ~ 1,
            random = Asym ~ 1)}
            
# fit nls object using simple logistic model, with 
# M, RAsym, Rk, Ri, and RM fixed to values in pnmodelparams
\donttest{
 data(logist_data)
 richardsR20.lis <- nls(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = 1, RAsym = 1, Rk = 1, Ri = 1, RM = 1 , modno = 20)
                        ,data = logist_data)}                        
}
