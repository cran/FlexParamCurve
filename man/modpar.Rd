\name{modpar}
\alias{modpar}
\title{Estimate Values to be Used for Fixed FlexParamCurve Parameters}
\description{This function creates the object \code{pnmodelparams}
which holds estimates of values for all 8 FlexParamCurve
parameters used for fitting and solving positive-negative Richards curves with
\code{\link{SSposnegRichards}} and \code{\link{posnegRichards.eqn}},
respectively. }
\usage{
modpar(x,
y,
first.y = NA,
x.at.first.y = NA,
last.y = NA,
x.at.last.y = NA,
twocomponent.x = NA,
verbose = FALSE,
force8par = FALSE,
force4par = FALSE,
pn.options = NA)
}
\arguments{
  \item{x}{a numeric vector of primary predictor variable}
  \item{y}{a numeric vector of response variable}
  \item{first.y}{the value of y at minimum x  when it is required to be constrained}
  \item{x.at.first.y}{the final value of x - 0 value is used if not specified when last.y is not NA}
  \item{last.y}{the value of y at maximum x when it is required to be constrained}
  \item{x.at.last.y}{the final value of x - must be specified if last.y is not NA}
  \item{twocomponent.x}{a numerical specifying the x-value (e.g. age) of intersection if a double model of 
two separate components is to be fitted. Alternatively a logical of value
= TRUE if the same type of model is to be fitted but the x of 
intersection is unknown}
  \item{verbose}{logical indicating whether information on successful optimization and
parameters should be returned during when using SSposnegRichards }
  \item{force8par}{logical specifying whether parameters of the negative Richards
curve should be set to defaults if they cannot be estimated}
  \item{force4par}{logical specifying whether parameters of the negative Richards
should be ignored - effectively using simple Richards curve}
 \item{pn.options}{character string specifying name of \code{\link{list}} object populated with starting 
parameter estimates, fitting options and bounds or the destination for \code{modpar} to write a new list}

}
\details{
This function creates a formatted \code{\link{list}} object as named by the argument pn.options. This list
holds estimates of values for all 8 FlexParamCurve parameters, fitting options and parameter bounds used for
fitting and solving double-Richards curves with \code{\link{SSposnegRichards}} and \code{\link{posnegRichards.eqn}},
respectively. Parameter bounds are the maximum and minimum parameters values that can be used by \code{\link{optim}}
and \code{\link{nls}} during parameter estimation. For definitions of parameters see either \code{\link{SSposnegRichards}} 
or \code{\link{posnegRichards.eqn}}. The list is written automatically by the function but 
it is also output as a return value for assignation in the usual way [myoptions<- change.pnparameters(...)].

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

This function can now fit biphasic (and more generally
double-Richards) curves, where the final curve is effectively either two positive curves
or two negative curves, as well as negative-positive curves. This functionality is default 
and does not need to be specified.

Parameter bounds estimated here for use in \code{\link{optim}} and \code{\link{nls}}
fits within \code{\link{SSposnegRichards}} are
applicable to a wide range of curves, although user may
change these manually in \code{\link{list}} object specified by \code{pn.options}.
Bounds are estimated by \code{\link{modpar}} by adding or subtracting multiples
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
or \code{\link{posnegRichards.eqn}}.

Use force4par = TRUE if you desire to produces estimates only for the four parameters of
a single Richards curves. This should also be used if you wish to fit simple logistic
Gompertz or von Bertalanffy curves: see \code{\link{SSposnegRichards}} for more details.

When specified, first.y and last.y are saved in list object specified by \code{pn.options} to instruct
\code{\link{SSposnegRichards}} to add this as the first or last value of the response, respectively,
during estimation.

To fit two-component double-curves, in which one curve equation is used up to (and including)
the x of intersection and a separate equation is used for x-values greater than the x of intersection
the argument twocomponent.x should be set to the value for the x of intersection. If this argument
is anything other than NA then a two-component model will be fitted when \code{\link{SSposnegRichards}}
is called. This option will be saved in list object specified by \code{pn.options} and can be changed at will.

If the argument verbose = TRUE then details concerning the optimization processes within 
\code{\link{SSposnegRichards}} are printed on screen whenever \code{\link{SSposnegRichards}} is called.
These include whether optimization of the first or second parts of the curve or simultaneous optimizations
are successful, if these have been further refined by nls, whether default parameters were used or the 
parameterization was aborted and what parameter values were finally exported by \code{\link{SSposnegRichards}}. 
This option will be saved in the list object specified  by \code{pn.options} and can be changed at will.
}
\value{a \code{\link{list}} of estimated fixed values for all
above arguments}

\author{Stephen Oswald <ozsao23@hotmail.com>}




\examples{
# estimate fixed parameters use data object posneg.data
   	 data(posneg.data)
   	 modpar(posneg.data$age, posneg.data$mass, pn.options = "myoptions")

# estimate fixed parameters use data object posneg.data (only first 
# 4 group levels for example's sake) and specify a fixed hatching 
# mass for curve optimization using \code{\link{SSposnegRichards}}
\donttest{
	 modpar(posneg.data$age, posneg.data$mass, pn.options = "myoptions")
   	 subdata <- subset(posneg.data,posneg.data$id == as.character(36)
	    		| posneg.data$id == as.character(9) 
	    		| posneg.data$id == as.character(32) 
	    		| posneg.data$id == as.character(43))
	 richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, 
	         K = K, Infl = Infl, RAsym = RAsym, Rk = Rk, Ri = Ri, 
        	 modno = 22, pn.options = "myoptions"), data = subdata)}
   		 
# fit a two component model - enter your own data in place of "mydata"
\donttest{
\dontrun{
 	modpar(mydata$x,mydata$y,twocomponent.x=TRUE, pn.options = "myoptions") # if x of intersection unknown
 	modpar(mydata$x,mydata$y,twocomponent.x=75, pn.options = "myoptions") # if x of intersection = 75
 	richardsR1.nls <- nls(y~ SSposnegRichards(x, Asym = Asym, K = K,
  		 Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1, pn.options = "myoptions")
                        , data = mydata)
}
}

# force an 8 parameter estimate on logistic data
\donttest{data(logist.data)
modpar(logist.data$age,logist.data$mass,force8par=TRUE, pn.options = "myoptions")
}

# force an 4 parameter model on logistic data
\donttest{data(logist.data)
modpar(logist.data$age,logist.data$mass,force4par=TRUE, pn.options = "myoptions")
}

# troubleshoot the fit of a model 
\donttest{data(posneg.data)
modpar(posneg.data$age,posneg.data$mass,verbose=TRUE, pn.options = "myoptions")
}

}
