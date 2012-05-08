\name{posnegRichards.eqn}
\alias{posnegRichards.eqn}
\title{Equations of the FlexParamCurve Family}
\description{Function to solve any of the equations in the FlexParamCurve family,
depending on user-specified parameters and model choice}
\usage{
posnegRichards.eqn(x,
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
  \item{x}{a numeric vector of the primary predictor variable}
  \item{Asym}{a numeric value for the asymptote of the positive (increasing) curve}
  \item{K}{a numeric value for the rate parameter of the positive (increasing) curve}
  \item{Infl}{a numeric value for the point of inflection of the positive (increasing) curve}
  \item{M}{a numeric value for the shape parameter of the positive (increasing) curve}
  \item{RAsym}{a numeric value for the asymptote of the negative (decreasing) curve}
  \item{Rk}{a numeric value for the rate parameter of the negative (decreasing) curve}
  \item{Ri}{a numeric value for the point of inflection of the negative (decreasing) curve}
  \item{RM}{a numeric value for the shape parameter of the negative (decreasing) curve}
  \item{modno}{a numeric value (currently integer only) between 1 and 36 specifying the identification
number of the equation to be fitted}
  \item{pn.options}{a character vector specifying a list of parameters and options for plotting}
}
\details{This function fits 1 of 32 possible FlexParamCurve equations (plus custon model #17). Equations 
can fit both monotonic and non-monotonic curves (two different trajectories).

These equations have also been described as double-Richards curves, or positive-negative Richards curves.
From version 1.2 onwards this function can fit curves that exhibit negative followed by positive
trajectories or double-positive or double-negative trajectories. This function can now also fit two 
component (biphasic) models, where the first curve is used up to the x-value (e.g. age) of intersection and the 
second curve is used afterwards, thus the curves are not joined as in standard models (see \code{\link{SSposnegRichards}}
for details.

The 32 possible equations are all based on the subtraction of one Richards curve from another, producing:
\eqn{y = A / ([1+ m exp(-k (t-i))]1/m) - A' / ([1+ m' exp(-k' (t-i' ))]1/m' )}, where A=Asym, k=K, i=Infl, m=M,
A'=RAsym, k'=Rk, i'=Ri, m'=RM; as described in the Arguments section above.

All 32 possible equations are simply reformulations of this equation, in each case fixing a parameter or
multiple parameters to (by default) the mean parameter across all individuals in the dataset (such as produced by a \code{\link{nls}}
model). All models are detailed in the \code{\link{SSposnegRichards}} help file. Any models that require parameter fixing
(i.e. all except model #1) extract appropriate values from the specified list passed by name to pn.options for the fixed parameters. 
This object is most easily created by running \code{\link{modpar}} and can be adjusted manually or by using 
\code{\link{change.pnparameters}} to user required specification.}
\value{the solution of the equation specified (by modno), given the user-entered parameters}

\author{Stephen Oswald <ozsao23@hotmail.com>}
\note{Any models that require parameter fixing (i.e. all except model #1) extract appropriate values from the specified
list passed to pn.options for the fixed parameters. This object is most easily created by running \code{\link{modpar}} 
and can be adjusted manually or by using \code{\link{change.pnparameters}} to user required specification.}

\seealso{\code{\link{SSposnegRichards}}
\code{\link{modpar}}    }
\examples{
require(graphics)
# calculate y (dependent variable) for a given x for an 8-parameter double-Richards model
   	 data(posneg.data)
   	 modpar(posneg.data$age, posneg.data$mass, pn.options = "myoptions") #create pnmodelparams for fixed parameters
   	 x = 10
   	 y <- posnegRichards.eqn(x, 1000, 0.5, 25, 1, 100, 0.5, 125, 
   	 1, modno = 1, pn.options = "myoptions")
   	 print( c(x = x, y = y) )

# calculate y (dependent variable) for a given x for a 4-parameter Richards model, note that second curve parameters are unneeded
# and replaced with value from pnmodelparams
	 x = 10
   	 y <- posnegRichards.eqn(x, 1000, 0.5, 25, 1, 
   	 1, modno = 12, pn.options = "myoptions")
   	 print( c(x = x, y = y) )
   	 
# plot a logistic curve (M=1),  note that second curve parameters are unneeded 	     
   	 plot(1:200, posnegRichards.eqn(1:200, Asym = 1000, K = 0.5, Infl = 25, M = 1, 
   	  modno = 12, pn.options = "myoptions"), xlim = c(1, 200), xlab = "x", 
   		 ylab = "y", pch = 1, cex = 0.7)

# plot a double-logistic curve (M=1, RM=1),  note that second curve parameters are unneeded 	     
   	 plot(1:200, posnegRichards.eqn(1:200, Asym = 1000, K = 0.5, Infl = 25, M = 1, RAsym = -100, Rk = 0.5, Ri = 125, RM = 1,
   	  modno = 1, pn.options = "myoptions"), xlim = c(1, 200), xlab = "x", 
   		 ylab = "y", pch = 1, cex = 0.7)
}