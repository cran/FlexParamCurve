\name{modpar}
\alias{modpar}
\title{Estimate Values to be Used for Fixed FlexParamCurve Parameters}
\description{This function creates the object \code{pnmodelparams}
which holds estimates of values for all 8 FlexParamCurve
parameters used for fitting and solving positive-negative Richards curves with
\code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
respectively. }
\usage{
modpar(x,
y,
first_y = NA,
x_at_first_y = NA,
last_y = NA,
x_at_last_y = NA,
twocomponent_age = NA,
verbose = FALSE,
force8par = FALSE,
force4par = FALSE)
}
\arguments{
  \item{x}{a numeric vector of primary predictor variable}
  \item{y}{a numeric vector of response variable}
  \item{first_y}{the value of y at minimum x  when it is required to be constrained}
  \item{x_at_first_y}{the final value of x - 0 value is used if not specified when last_y is not NA}
  \item{last_y}{the value of y at maximum x when it is required to be constrained}
  \item{x_at_last_y}{the final value of x - must be specified if last_y is not NA}
  \item{twocomponent_age}{a numerical specifying the age of intersection if a double model of 
two separate components is to be fitted. Alternatively a logical of value
= TRUE if the same type of model is to be fitted but the age of 
intersection is unknown}
  \item{verbose}{logical indicating whether information on successful optimization and
parameters should be returned during when using SSposnegRichards }
  \item{force8par}{logical specifying whether parameters of the negative Richards
curve should be set to defaults if they cannot be estimated}
  \item{force4par}{logical specifying whether parameters of the negative Richards
should be ignored - effectively using simple Richards curve}
}
\details{
This function creates the object \code{pnmodelparams}
which holds estimates of values for all 8 FlexParamCurve
parameters used for fitting and solving positive-negative Richards curves with
\code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
respectively. Running this function concurrently creates the object
\code{pnmodelparamsbounds} which holds the maximum and minimum values
for parameters to be used by \code{\link{optim}} and \code{\link{nls}}
during parameter estimation. For definitions of parameters see either
\code{\link{SSposnegRichards}} or \code{\link{posnegRichards_eqn}}. These
objects are updated globally and
do not need to be assigned by the user: the output of a return value simply
offers a way to save values elsewhere for posterity.

Estimates are produced by fitting positive-negative or double Richards curves in
\code{\link{nls}} using
\code{\link{SSposnegRichards}} for the full 8 parameter model (R1).
If this fails, the function \code{\link{getInitial}} is called to
attempt to produce initial estimates using the same 8 parameter model.
If this also fails, estimates are attempted in the same way using the
4 parameter (positive only) model (R12). In this case, only the positive
parameters are returned (NAs are substituted for negative parameters)
unless argument force8par=TRUE, in which case negative parameters are
defaulted to: RAsym = 0.05*Asym, Rk = K, Ri = Infl, RM = M.

This function can now fit bilogistic (and more generally
biRichards) curves, where the final curve is effectively either two positive curves
or two negative curves, as well as negative-positive curves. This functionality is default 
and does not need to be specified.

Parameter bounds estimated here for use in \code{\link{optim}} and \code{\link{nls}}
fits within \code{\link{SSposnegRichards}} are
applicable to a wide range of curves, although user may
change these manually in \code{\link{list}} object \code{pnmodelparamsbounds}
Bounds are estimated by \code{\link{modpar}} by adding or subtracting  multiples
of fixed parameter values to estimated mean parameter values:
-Asym*0.5 and +Asym*2.5,
-K*0.5 and +K*0.5,
-Infl*2.5 and +Infl*10
-M*2 and +M*2
-RAsym*0.5 and +RAsym*2.5,
-Rk*0.5 and +Rk*0.5,
-Ri*2.5 and +Ri*5
-RM*2 and +RM*2.

Use force8par = TRUE if initial call to \code{\link{modpar}} produces estimates for
only 4 parameters and yet an 8 parameter model is desired for \code{\link{SSposnegRichards}}
or \code{\link{posnegRichards_eqn}}.

Use force4par = TRUE if you desire to produces estimates only for the four parameters of
a single Richards curves. This should also be used if you wish to fit simple logistic
Gompertz or von Bertalanffy curves: see \code{\link{SSposnegRichards}} for more details.

When specified, first_y and last_y are saved in \code{pnmodelparams} to instruct
\code{\link{SSposnegRichards}} to add this as the first or last value of the response, respectively,
during estimation.

To fit two-component double-curves, in which one curve equation is used up to (and including)
the age of intersection and a separate equation is used for ages greater than the age of intersection
the argument twocomponent_age should be set to the value for the age of intersection. If this argument
is anything other than NA then a two-component model will be fitted when \code{\link{SSposnegRichards}}
is called. This option will be saved in \code{pnmodelparams} and can be changed at will.

If the argument verbose = TRUE then details concerning the optimization processes within 
\code{\link{SSposnegRichards}} are printed on screen whenever \code{\link{SSposnegRichards}} is called.
These include whether optimization of the positive, negative parts of the curve or simultaneous optimizations
are successful, if these have been further refined by nls, whether default parameters were used or the 
parameterization was aborted and what parameter values were finally exported by \code{\link{SSposnegRichards}}. 
This option will be saved in \code{pnmodelparams} and can be changed at will.
}
\value{a \code{\link{list}} of estimated fixed values for all
above arguments}

\author{Stephen Oswald <ozsao23@hotmail.com>}




\examples{
# estimate fixed parameters use data object posneg_data
   	 data(posneg_data)
   	 modpar(posneg_data$age, posneg_data$mass)

# estimate fixed parameters use data object posneg_data (only first 
# 4 group levels for example's sake) and specify a fixed hatching 
# mass for curve optimization using \code{\link{SSposnegRichards}}
\donttest{
	 modpar(posneg_data$age, posneg_data$mass)
   	 subdata <- subset(posneg_data,posneg_data$id == as.character(36)
	    		| posneg_data$id == as.character(9) 
	    		| posneg_data$id == as.character(32) 
	    		| posneg_data$id == as.character(43))
	 richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, 
	         K = K, Infl = Infl, RAsym = RAsym, Rk = Rk, Ri = Ri, 
        	 modno = 22), data = subdata)}
   		 
# fit a two component model - enter your own data in place of "mydata"
\donttest{
\dontrun{
 	modpar<-(mydata$x,mydata$y,twocomponent_age=TRUE) # if age of intersection unknown
 	modpar<-(mydata$x,mydata$y,twocomponent_age=75) # if age of intersection = 75
 	richardsR1.nls <- nls(y~ SSposnegRichards(x, Asym = Asym, K = K,
  		 Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1)
                        , data = mydata)
}
}

# force an 8 parameter estimate on logistic data
\donttest{data(logist_data)
modpar<-(logist_data$age,logist_data$mass,force8par=TRUE)
}

# force an 4 parameter model on logistic data
\donttest{data(logist_data)
modpar<-(logist_data$age,logist_data$mass,force4par=TRUE)
}

# troubleshoot the fit of a model 
\donttest{data(posneg_data)
modpar<-(posneg_data$age,posneg_data$mass,verbose=TRUE)
}

}
