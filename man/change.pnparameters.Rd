\name{change.pnparameters}
\alias{change.pnparameters}
\title{Change Fixed Parameter Values}
\description{Function to alter values of parameters to be used by \code{\link{SSposnegRichards}}
or \code{\link{posnegRichards.eqn}}
as the fixed values in equations where parameters are fixed}
\usage{
change.pnparameters(Asym = NA,
K = NA,
Infl = NA,
M = NA,
RAsym = NA,
Rk = NA,
Ri = NA,
RM = NA,
Amin = NA,
Amax = NA,
Kmin = NA,
Kmax = NA,
Imin = NA,
Imax = NA,
Mmin = NA,
Mmax = NA,
RAmin = NA,
RAmax = NA,
Rkmin = NA,
Rkmax = NA,
Rimin = NA,
Rimax = NA,
RMmin = NA,
RMmax = NA,
first.y = NA,
x.at.first.y = NA,
last.y = NA,
x.at.last.y = NA,
twocomponent.x = NA,
verbose = NA,
force4par = NA,
pn.options,
Envir = .GlobalEnv)
}
\arguments{
  \item{Asym}{a numeric value for the asymptote of the positive (increasing) curve}
  \item{K}{a numeric value for the rate parameter of the positive (increasing) curve}
  \item{Infl}{a numeric value for the point of inflection of the positive (increasing) curve}
  \item{M}{a numeric value for the shape parameter of the positive (increasing) curve}
  \item{RAsym}{a numeric value for  the asymptote of the negative (decreasing) curve}
  \item{Rk}{a numeric value for the rate parameter of the negative (decreasing) curve}
  \item{Ri}{a numeric value for the point of inflection of the negative (decreasing) curve}
  \item{RM}{a numeric value for the shape parameter of the negative (decreasing) curve}
    \item{Amin}{a numeric value for the minimum bound of Asym}
    \item{Amax}{a numeric value for the maximum bound of Asym}
    \item{Kmin}{a numeric value for the minimum bound of K}
    \item{Kmax}{a numeric value for the maximum bound of K}
    \item{Imin}{a numeric value for the minimum bound of Infl}
    \item{Imax}{a numeric value for the maximum bound of Infl}
    \item{Mmin}{a numeric value for the minimum bound of M}
    \item{Mmax}{a numeric value for the maximum bound of M}
    \item{RAmin}{a numeric value for the minimum bound of RAsym}
    \item{RAmax}{a numeric value for the maximum bound of RAsym}
    \item{Rkmin}{a numeric value for the minimum bound of Rk}
    \item{Rkmax}{a numeric value for the maximum bound of Rk}
    \item{Rimin}{a numeric value for the minimum bound of Ri}
    \item{Rimax}{a numeric value for the maximum bound of Ri}
    \item{RMmin}{a numeric value for the minimum bound of RM}
    \item{RMmax}{a numeric value for the maximum bound of RM}
  \item{first.y}{the value of y at minimum x when it is required to be constrained}
  \item{x.at.first.y}{the final value of x - 0 value is used if not specified when last.y is not NA}
  \item{last.y}{the value of y at maximum x when it is required to be constrained}
  \item{x.at.last.y}{the final value of x - this is option is currently disabled}
  \item{twocomponent.x}{a numerical specifying the x-value (e.g. age) of intersection if a double model of 
two separate components is to be fitted. Alternatively a logical of value
= TRUE if the same type of model is to be fitted but the x of 
intersection is unknown}
  \item{verbose}{logical indicating whether information on successful optimization and
parameters should be returned during when using SSposnegRichards }
  \item{force4par}{logical specifying whether parameters of the negative Richards
should be ignored - effectively using simple Richards curve}
  \item{pn.options}{required character string specifying the name of a list object currently 
  populated with starting parameter estimates, fitting options and bounds to be modified}
   \item{Envir}{a valid R environment to find pn.options in, by default this is the global environment} 
}
\details{This function provides a simple way for the user to update
the a user-named list that holds fixed values and options
for fitting and solving positive-negative Richards curves with
\code{\link{SSposnegRichards}} and \code{\link{posnegRichards.eqn}},
respectively. Running this function also concurrently updates the parameterbounds
in the same list which are vthe maximum and minimum values
for parameters to be used by \code{\link{optim}} and \code{\link{nls}}
during parameter estimation
in \code{\link{SSposnegRichards}}. The list is written automatically by the function but 
it is also output as a return value for assignation in the usual way [myoptions<- change.pnparameters(...)].

The list specified by pn.options must exist before this function is called. Use \code{\link{modpar}}
to estimate values for all parameters and easily generate a suitable list. See  \code{\link{modpar}} for details of bounding.}
\value{a \code{\link{list}} of values for all above arguments,
with new values substituted where specified in the call}

\author{Stephen Oswald <steve.oswald@psu.edu>}
\note{Requires \code{\link{modpar}} to be have been run prior to execution   

  Version 1.5 saves many variables, and internal variables in the package environment:
  FlexParamCurve:::FPCEnv. By default, the pn.options file is copied to the environment
  specified by the functions (global environment by default). Model selection routines
  also copy from FPCenv to the global environment all nlsList models fitted during 
  model selection to provide backcompatibility with code for earlier versions. The user
  can redirect the directory to copy to by altering the Envir argument when calling the 
  function.
 }


\seealso{\code{\link{modpar}} \code{\link{SSposnegRichards}} \code{\link{posnegRichards.eqn}}}
\examples{
# change all fixed values except K and Rk
    modpar(posneg.data$age, posneg.data$mass, pn.options = "myoptions")
    change.pnparameters(Asym = 10000, Infl = 80, M = 5, RAsym = 10000, 
        Ri = 240, RM = 5, pn.options = "myoptions")

# change fixed values of M and constrain hatching mass to 45.5 in a growth curve
    change.pnparameters(M = 1, RM = 0.5, first.y = 45.5, pn.options = "myoptions")
}
