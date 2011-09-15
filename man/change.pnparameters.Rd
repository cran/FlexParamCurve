\name{change.pnparameters}
\alias{change.pnparameters}
\title{Change Fixed Parameter Values}
\description{Function to alter values of parameters to be used by \code{\link{SSposnegRichards}}
or \code{\link{posnegRichards_eqn}}
as the fixed values in equations where parameters are fixed}
\usage{
change.pnparameters(Asym = NA,
K = NA,
Infl = NA,
M = NA,
RAsym = NA,
RK = NA,
Ri = NA,
RM = NA,
first_y = NA,
last_y = NA,
x_at_last_y = NA)
}
\arguments{
  \item{Asym}{a numeric value for the asymptote of the positive (increasing) curve}
  \item{K}{a numeric value for the rate parameter of the positive (increasing) curve}
  \item{Infl}{a numeric value for the point of inflection of the positive (increasing) curve}
  \item{M}{a numeric value for the shape parameter of the positive (increasing) curve}
  \item{RAsym}{a numeric value for  the asymptote of the negative (decreasing) curve}
  \item{RK}{a numeric value for the rate parameter of the negative (decreasing) curve}
  \item{Ri}{a numeric value for the point of inflection of the negative (decreasing) curve}
  \item{RM}{a numeric value for the shape parameter of the negative (decreasing) curve}
  \item{first_y}{the value of y at minimum x when it is required to be constrained}
  \item{last_y}{the value of y at maximum x when it is required to be constrained}
  \item{x_at_last_y}{the final value of x - this is option is currently disabled}
}
\details{This function provides a simple way for the user to update
the object \eqn{pnmodelparams} that holds fixed values and options
for fitting and solving positive-negative Richards curves with
\code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
respectively. Running this function concurrently updates the object
\eqn{pnmodelparamsbounds} which holds the maximum and minimum values
for parameters to be used by \code{\link{optim}} and \code{\link{nls}}
during parameter estimation
in \code{\link{SSposnegRichards}}. These objects are updated globally and
do not need to be assigned by the user: the output of a return value simply
offers a way to save values elsewhere for posterity.

Both objects must exist before this function is called. Use \code{\link{modpar}}
to estimate values for all parameters and create \eqn{pnmodelparams} and
\eqn{pnmodelparamsbounds}. See  \code{\link{modpar}} for details of bounding.}
\value{a \code{\link{list}} of values for all above arguments,
with new values substituted where specified in the call}

\author{Stephen Oswald <ozsao23@hotmail.com>}
\note{Requires \code{\link{modpar}} to be have been run prior to execution}


\seealso{\code{\link{modpar}} \code{\link{SSposnegRichards}} \code{\link{posnegRichards_eqn}}}
\examples{
# change all fixed values except K and Rk
  data(posneg_data)
  modpar(posneg_data$age,posneg_data$mass)
  change.pnparameters(Asym=10000,Infl=80,M=5,RAsym=10000,Ri=240,RM=5)

# change fixed values of M and constrain hatching mass to 45.5 in a growth curve
  change.pnparameters(M=1,RM=0.5,first_y=45.5)
}
