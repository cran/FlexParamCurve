\name{modpar}
\alias{modpar}
\title{Estimate Values to be Used for Fixed FlexParamCurve Parameters}

\usage{
modpar(x,
y,
first_y = NA,
last_y = NA,
x_at_last_y = NA,
force8par = FALSE)
}
\arguments{
  \item{x}{a numeric vector of primary predictor variable}
  \item{y}{a numeric vector of response variable}
  \item{first_y}{the value of y at minimum x  when it is required to be constrained}
  \item{last_y}{the value of y at maximum x when it is required to be constrained}
  \item{x_at_last_y}{the final value of x - this is option is currently disabled}
  \item{force8par}{logical specifying whether parameters of the negative Richards
curve should be set to defaults if they cannot be estimated}
}
\description{This function creates the object \eqn{pnmodelparams}
which holds estimates of values for all 8 FlexParamCurve
parameters used for fitting and solving positive-negative Richards curves with
\code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
respectively.}
\details{Running this function concurrently creates the object
\eqn{pnmodelparamsbounds} which holds the maximum and minimum values
for parameters to be used by \code{\link{optim}} and \code{\link{nls}}
during parameter estimation. For definitions of parameters see either
\code{\link{SSposnegRichards}} or \code{\link{posnegRichards_eqn}}. These
objects are updated globally and
do not need to be assigned by the user: the output of a return value simply
offers a way to save values elsewhere for posterity.

Estimates are produced by fitting positive-negative Richards curves in
\code{\link{nls}} using
\code{\link{SSposnegRichards}} for the full 8 parameter model (R1).
If this fails, the function \code{\link{getInitial}} is called to
attempt to produce initial estimates using the same 8 parameter model.
If this also fails, estimates are attempted in the same way using the
4 parameter (positive only) model (R12). In this case, only the positive
parameters are returned (NAs are substituted for negative parameters)
unless argument force8mod=TRUE, in which case negative parameters are
defaulted to: RAsym = 0.05*Asym, Rk = K, Ri = Infl, RM = M.

Parameter bounds estimated here for use in \code{\link{optim}} and \code{\link{nls}}
fits within \code{\link{SSposnegRichards}} are
applicable to a wide range of curves, although user may
change these manually in \code{\link{list}} object \eqn{pnmodelparamsbounds}
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

Use force8mod=TRUE if initial call to \code{\link{modpar}} produces estimates for
only 4 parameters and yet an 8 parameter model is desired for \code{\link{SSposnegRichards}}
or \code{\link{posnegRichards_eqn}}.

When specified, first_y and last_y are saved in \eqn{pnmodelparams} to instruct
\code{\link{SSposnegRichards}} to add this as the first or last value of the response, respectively,
during estimation.}
\value{a \code{\link{list}} of estimated fixed values for all
above arguments}

\author{Stephen Oswald <ozsao23@hotmail.com>}




\examples{
    # estimate fixed parameters use data object posneg_data
        data(posneg_data)
        modpar(posneg_data$age,posneg_data$mass)

    # estimate fixed parameters use data object posneg_data (only first 
    # 4 group levels for example's sake) and specify a fixed hatching 
    # mass for curve optimization using \code{\link{SSposnegRichards}}
        modpar(posneg_data$age,posneg_data$mass)
        subdata<-subset(posneg_data, as.numeric(row.names (posneg_data) ) < 53)
        richardsR1.lis<-nlsList(mass~SSposnegRichards(age,Asym=Asym,K=K,Infl=Infl,M=M,RAsym=RAsym,Rk=Rk,Ri=Ri,RM=RM,modno=1)
                        ,data=subdata)
}
