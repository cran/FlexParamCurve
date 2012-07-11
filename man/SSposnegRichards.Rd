\name{SSposnegRichards}
\alias{SSposnegRichards}
\title{Self-Starting Positive-Negative Richards Model (double-Richards)}

\usage{
SSposnegRichards(x,
Asym = NA,
K = NA,
Infl = NA,
M = NA,
RAsym = NA,
Rk = NA,
Ri = NA,
RM = NA,
modno,
pn.options)
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
  \item{pn.options}{character string specifying name of 
\code{\link{list}} object populated with starting 
parameter estimates, fitting options and bounds}

}
\description{This selfStart function evaluates a range of flexible logistic
    functions.  It also has an initial attribute that creates
    initial estimates of the parameters
    for the model specified.}
\details{
     This selfStart function evaluates a range of flexible logistic
     functions.  It also has an initial attribute that creates
     initial estimates of the parameters
     for the model specified. Equations can fit both monotonic and non-monotonic 
     curves (two different trajectories). These equations have also been described as
     double-Richards curves, or positive-negative Richards curves. **From version 1.2 
     onwards this function can fit curves that exhibit negative followed by positive
     trajectories or double positive or double negative trajectories.***
    
     The 32 possible equations (plus custom model #17) are all based on the subtraction of one Richards
     curve from another, producing:
     \eqn{y = A / ([1+ m exp(-k (t-i))]1/m) + A' / ([1+ m' exp(-k' (t-i' ))]1/m' )},
     where A=Asym, k=K, i=Infl, m=M,
     A'=RAsym, k'=Rk, i'=Ri, m'=RM; as described in the Arguments section above.
    
     All 32 possible equations are simply reformulations of this equation, in each
     case fixing a parameter or multiple parameters to (by default) the mean parameter across
     all individuals in the dataset (such as produced by a \code{\link{nls}}
     model). Thus, a model in which one parameter is fixed has a 7-parameter equation,
     and one in which four are fixed has a 4-parameter equation, thus reducing 
     complexity and computation and compensatory parameter changes when a parameter does not 
     vary across group levels (e.g individuals)
     [the most appropriate equation can be determined using model selection in
     \code{\link{pn.modselect.step}} or \code{\link{pn.mod.compare}}].
     Any models that require parameter fixing (i.e. all except #1)
     extract appropriate values from the list object specified by \code{pn.options} for the fixed
     parameters. This object is most easily created by running
     \code{\link{modpar}} and can be adjusted manually or by using
     \code{\link{change.pnparameters}} to user-required specification.
    
     Each of the 32 equations is identified by an integer value for modno (1 to 36).Models
     21-36 are the same as 1-16 except that in the former the first curve parameter \code{m}
     is fixed and not estimated. All equations (except 17 - see below) contain parameters Asym, K, and Infl.
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
     7       Y           Y    Y    Y
     8       Y     Y     Y         Y
     9       Y           Y         Y
     10      Y                Y
     11      Y     Y
     12      Y			 4 parameter, standard Richards model
     13      Y     Y          Y
     14      Y           Y    Y
     15      Y     Y     Y
     16      Y           Y
     17      see below
     18      see below
     19      see below
     20      see below                           
     21            Y     Y    Y    Y    7 parameter model, 4 recession params
     22            Y     Y    Y		6 parameters (used for double-logistic, double-Gompertz or double-Von Bertalannfy)
     23                       Y    Y
     24            Y               Y
     25                            Y
     26            Y          Y    Y
     27                  Y    Y    Y
     28            Y     Y         Y
     29                  Y         Y
     30                       Y
     31            Y
     32					only 3 parameters (used for logistic, Gompertz or Von Bertalannfy - see below)
     33            Y          Y
     34                  Y    Y
     35            Y     Y
     36                  Y           }
    
     modno 17 represents a different parameterization for a custom model:
     (Asym/ 1 + exp(Infl - x)/ M) - (RAsym / 1 + exp(Ri - x)/ RM), in
     which M and RM actually represent scale parameters not shape parameters.
    
     modnos 18, 19 and 20 are reserved for internal use by \code{\link{modpar}}.
    
     To access common 3 parameter sigmoidal models use modno = 20, fixing
     parameters (using \code{\link{change.pnparameters}}) to M = 1 for logistic,
     M = 0.1 for Gompertz, and M= -0.3 for von Bertalanffy. The same settings can be
     used with modno = 2 to fit the double-logistic, double-Gompertz or double-Von Bertalannfy.
     Note that to fit only 3 or 4 parameter curves, the option force4par = TRUE should be specified
     when running \code{\link{modpar}}.
    
     The call for \code{\link{SSposnegRichards}} only differs from
     conventional selfStart models in that it requires a value for modno and a list of fitting options
     and values from \code{\link{modpar}} to which to fix parameters in the reduced models. 
     Depending on the model chosen, different combinations of the 8 possible
     parameters are required: if one is missing the routine will stop with an
     appropriate error, if an extra one is added, it will be ignored (provided
     that it is labelled, e.g. M = 1; this is good practice to prevent accidental
     misassignments). 
     
     Here are two examples (7 parameter and 3 parameter):
     \preformatted{richardsR2.lis <- nlsList(mass ~ SSposnegRichards(age,
         Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri,
         , modno = 2), data = posneg.data)
         #correct call includes all necessary parameters}
         
     \preformatted{richardsR20.lis <- nlsList(mass ~ SSposnegRichards(age,
         Asym = Asym, K = K, Infl = Infl, modno = 2), data = logist.data)
         #incorrect call missing required parameters, function terminates and generates an error message}
    
     Examples for all models can be found in the list object
     \code{\link{posnegRichards.calls}}.
    
     If specified using \code{\link{modpar}} optional constraints
     may be placed to specify response values at the minimum value and/or
     maximum values of the predictor. Such constraint allows realistic fits for
     datasets which are missing data
     at either end of the curve (e.g. hatching weight for some growth curves).
    
     Estimates are produced by splitting the two curves into separate positive
     and negative curves and estimating parameters for each curve separately
     in a similar manner to \code{\link{SSlogis}}. Each curve is fit first by
     \code{\link{optim}} using the parameter bounds in \code{pnmodelparamsbounds}
     (see \code{\link{modpar}}) and a subsequent refinement is attempted using
     \code{\link{nls}} with more restrictive parameter bounds. Finally, both curves
     are annealed
     and parameters are again estimated using restrictive bounds and starting values
     already determined during separate estimates. Equations for which the positive
     curve was inestimable are not estimated further, but if negative curve estimation
     or overall curve estimation fail, partial estimates are used: either default
     negative parameters (RAsym = 0.05*Asym, Rk = K, Ri = Infl, RM = M) annealed to positive
     curves or separate estimates annealed; both with compensation for interation
     between asymptotes.
    
     From version 1.2 onwards, this function can now fit two component models, where
     the first curve is used up to the x-value (e.g. age) of intersection and the second curve is used
     afterwards. Confusingly, these are also called "double Richards", "double Gompertz"
     or "double logistic": see Murphy et al. (2009) or Ross et al. (1995) for examples.
     To specify such models set twocomponent.x = TRUE (this will estimate the x of 
     intersection) when running \code{\link{modpar}}. Alternatively, if known, the 
     x of intersection can be set directly by setting
     twocomponent.x = # (where # is the x of intersection). Whne \{code{\link{modpar}}
     is run this option will be saved in \code{pnmodelparas} and can be changed at will,
     either manually or using \code{\link{change.pnparameters}}.
    
     From version 1.2 onwards, this function can now fit bilogistic (and more generally
     biRichards) curves, where the final curve is effectively either two positive curves
     or two negative curves. See Meyer (1994) for examples. This functionality is default 
     and does not need to be specified. 
}
\value{a numeric vector of the same length as x containing parameter
     estimates for equation specified (fixed parameters
     are not return but are substituted in calls to \code{\link{nls}}
     \code{\link{nlsList}} and \code{\link{nlme}} with the fixed parameters
     stored in \code{pnmodelparams}; see \code{\link{modpar}}}
\references{## Oswald, S.A. et al. (in review) FlexParamCurve: R package for flexible
     fitting of nonlinear parametric curves. Application paper submitted to Methods
     in Ecology and Evolution.
	#1# Nelder, J.A. (1962) Note: an alternative form of a generalized
     logistic equation. Biometrics, 18, 614-616.
         #2#
     Huin, N. & Prince, P.A. (2000) Chick growth in albatrosses: curve fitting
     with a twist. Journal of Avian Biology, 31, 418-425.
         #3#
     Meyer, P. (1994) Bi-logistic growth. Technological Forecasting and Social 
     Change. 47: 89-102
         #4#
     Murphy, S. et al. (2009) Importance of biological parameters in assessing 
     the status of \code{Delphinus delphis}. Marine Ecology Progress Series 388: 273-291.
         #5#
     Pinheiro, J. & Bates, D. (2000) \code{Mixed-Effects Models in S and S-Plus.}
     Springer Verlag, Berlin.
         #6#
     Ross, J.L. et al. (1994) Age, growth, mortality, and reproductive biology
     of red drums in North Carolina waters. Transactions of the American Fisheries
     Society 124: 37-54.}
\author{Stephen Oswald <ozsao23@hotmail.com>}
\note{Any models that require parameter fixing (i.e. all except #1)
     extract appropriate values from the object \code{pnmodelparams} for the fixed
     parameters. This object is created by running
     \code{\link{modpar}} and can be adjusted manually or by using
     \code{\link{change.pnparameters}} to user required specification.
    
     Output may show errors and warnings especially during a \code{\link{nlsList}}
     fit, in which the function is called repeatedly: once for each group level in the
     dataset. Warnings indicate conditions for which default parameters or incomplete estimates
     are used - see Details section - and errors occur from insufficient data or singularities.
     As a result of possible interaction and correlation between the parameters in some models,
     singularities may be common, but do not be alarmed by repeated error messages, as
     examination of a fitted \code{\link{nlsList}} model may releave a large number of
     well estimated group levels, thus the elimation of unsuitable outlying groups only. Also,
     because very few of the 32 equations are likely to be suitable for the majority
     of datasets, consideration of the model being fitted is crucial when examining the output. Functions
     \code{\link{pn.modselect.step}} and \code{\link{pn.mod.compare}} provide the ability
     for model selection of these equations through stepwise backward deletion or all
     model comparison, respectively. These offer powerful ways to determine the
     best equation for your dataset.
    
     To increase the ability of optimization routines to deal with
     a wide variety of values, particularly negative values for M or RM,
     only real component of complex numbers are modelled and integer versions
     of M and RM are used during estimation if floating values cause conversion
     issues.
    
     Speed of the function depends on the complexity of the data being fit.
}


\seealso{\code{\link{SSlogis}}
\code{\link{SSgompertz}}
\code{\link{posnegRichards.eqn}}}
\examples{
 set.seed(3) #for compatability issues
 require(graphics)
    # retrieve mean estimates of 8 parameters using getInitial
    # and posneg.data object
    data(posneg.data)
    modpar(posneg.data$age, posneg.data$mass,verbose=TRUE, pn.options = "myoptions")
    getInitial(mass ~ SSposnegRichards(age, Asym, K, Infl, M, 
        RAsym, Rk, Ri, RM, modno = 1, pn.options = "myoptions"), data = posneg.data)

    # retrieve mean estimates and produce plot to illustrate fit for 
    # curve with RM fixed
    pars <- coef(nls(mass ~ SSposnegRichards(age, 
        Asym = Asym, K = K, Infl = Infl, RAsym = RAsym, Rk = Rk, M = M,
        	Ri = Ri, modno = 2, pn.options = "myoptions"), data = posneg.data))
    plot(posneg.data$age, posneg.data$mass, pch=".")
    curve(posnegRichards.eqn(x, Asym = pars[1], K = pars[2], 
        Infl = pars[3], M = pars[4], RAsym = pars[5], Rk = pars[6], 
        Ri = pars[7],  modno = 2, pn.options = "myoptions"), xlim = c(0, 
        200), add = TRUE)
    
    # retrieve mean estimates and produce plot to illustrate fit for custom model 17
    \donttest{
    \dontrun{    pars<-as.numeric( getInitial(mass ~ SSposnegRichards(age, Asym, K, Infl,
           M, RAsym, Rk, Ri, RM, modno = 17, pn.options = "myoptions"), data = datansd) )
     plot(datansd$jday21March, datansd$moosensd)
     curve( posnegRichards.eqn(x, Asym = pars[1], K = 1, Infl = pars[2], 
            M = pars[3], RAsym = pars[4], Rk = 1, Ri = pars[5], RM = pars[6], 
            modno = 17, pn.options = "myoptions"), lty = 3, xlim = c(0, 200) , add = TRUE)}}
        
    # fit nls object using 8 parameter model
    # note: ensure data object is a groupedData object
    richardsR1.nls <- nls(mass ~ SSposnegRichards(age, Asym = Asym, 
        K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, 
        RM = RM, modno = 1, pn.options = "myoptions"), data = posneg.data)
        
    # fit nlsList object using 8 parameter model
    # note: ensure data object is a groupedData object
    # also note: not many datasets require all 8 parameters
    \donttest{
    \dontrun{
    richardsR1.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, 
        K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, 
        RM = RM, modno = 1, pn.options = "myoptions"), data = posneg.data)
    summary(richardsR1.lis)}}

    # fit nlsList object using 6 parameter model with value M and RM
    # fixed to value in pnmodelparams and then fit nlme model
    # note data is subset to provide estimates for a few individuals
    # as an example
    subdata <- subset(posneg.data,posneg.data$id == as.character(36)
   		| posneg.data$id == as.character(9) 
   		| posneg.data$id == as.character(32) 
   		| posneg.data$id == as.character(43))
    richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, 
        K = K, Infl = Infl, RAsym = RAsym, Rk = Rk, Ri = Ri, 
        modno = 22, pn.options = "myoptions"), data = subdata)
    summary(richardsR22.lis )
    richardsR22.nlme <- nlme(richardsR22.lis, random = pdDiag(Asym + Infl ~ 1) )
    summary(richardsR22.nlme)
         
    # fit nls object using simple logistic model, with 
    # M, RAsym, Rk, Ri, and RM fixed to values in pnmodelparams
    data(logist.data)
    modpar(logist.data$age, logist.data$mass ,force4par = TRUE, pn.options = "myoptions")
    change.pnparameters(M = 1, pn.options = "myoptions") #set to logistic (M =1) prior to fit
    richardsR32.nls <- nls(mass ~ SSposnegRichards(age, Asym = Asym, 
        K = K, Infl = Infl, modno = 32, pn.options = "myoptions"), data = logist.data)
    coef(richardsR32.nls)
                
    # fit a two component model - enter your own data in place of "mydata"
    \donttest{
    \dontrun{
    modpar(mydata$x,mydata$y,twocomponent.x=TRUE, pn.options = "myoptions") # if x of intersection unknown
    modpar(mydata$x,mydata$y,twocomponent.x=75, pn.options = "myoptions") # if x of intersection = 75
    richardsR1.nls <- nls(y~ SSposnegRichards(x, Asym = Asym, K = K,
      Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1, pn.options = "myoptions")
                       , data = mydata)
    coef(richardsR1.nls)}
    }}
