\name{posnegRichards_eqn}
\alias{posnegRichards_eqn}
\title{Equations of the FlexParamCurve Family}
\description{Function to solve any of the equations in the FlexParamCurve family,
depending on user-specified parameters and model choice}
\usage{
posnegRichards_eqn(x,
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
}
\details{This function fits 1 of 36 possible FlexParamCurve equations. All equations are biomodal,
in that they have a positive (increasing) trajectory followed by a negative (decreasing) trajectory.

These equations have also been described as double-Richards curves, or positive-negative Richards curves.

The 36 possible equations are all based on the subtraction of one Richards curve from another, producing:
\eqn{y = A / ([1+ m exp(-k (t-i))]1/m) - A' / ([1+ m' exp(-k' (t-i' ))]1/m' )}, where A=Asym, k=K, i=Infl, m=M,
A'=RAsym, k'=Rk, i'=Ri, m'=RM; as described in the Arguments section above.

All 36 possible equations are simply reformulations of this equation, in each case fixing a parameter or
multiple parameters to the mean parameter across all individuals in the dataset (such as produced by a \code{\link{nls}}
model). All models are detailed in the \code{\link{SSposnegRichards}} help file. Any models that require parameter fixing
(i.e. all except R1)
extract appropriate values from the object \eqn{pnmodelparams} for the fixed parameters. This object is created by running
\code{\link{modpar}} and can be adjusted manually or by using \code{\link{change.pnparameters}} to user required specification.}
\value{the solution of the equation specified (by modno), given the user-entered parameters}

\author{Stephen Oswald <ozsao23@hotmail.com>}
\note{Any models that require parameter fixing (i.e. all except #1)
extract appropriate values from the object \eqn{pnmodelparams} for the fixed parameters. This object is created by running
\code{\link{modpar}} and can be adjusted manually or by using \code{\link{change.pnparameters}} to user required specification.}


\seealso{\code{\link{SSposnegRichards}}
\code{\link{modpar}}}
\examples{
    require(graphics)
    # calculate y from an 8-parameter model
    data(posneg_data)
    modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
    y <- posnegRichards_eqn(10, 1000, 0.5, 25, 1, 100, 0.5, 125, 1, modno = 1)
    
    # calculate y from a 4-parameter positive only model, note that all negative parameters specified here are ignored
    # and replaced with value from pnmodelparams
    y <- posnegRichards_eqn(10 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12)
    
    # plot a logistic curve (M=1), note that all negative parameters specified here are ignored
    plot(1:200 ,posnegRichards_eqn(1:200 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12),xlim=c(1, 200),
    xlab = "x", ylab = "y",pch = 1, cex = 0.7)
}
