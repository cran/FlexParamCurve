\name{pn.mod.compare}
\alias{pn.mod.compare}
\title{Compare All Possible Positive-Negative Richards \eqn{nlslist} Models}

\usage{
pn.mod.compare(x,
y,
grp,
forcemod = 0,
existing = FALSE,
penaliz = "1/sqrt(n)",
pn.options)
}
\arguments{
  \item{x}{a numeric vector of the primary predictor}
  \item{y}{a numeric vector of the response variable}
  \item{grp}{a factor of same length as x and y that distinguishes groups within
the dataset}
  \item{forcemod}{optional numeric value to constrain model selection (see Details)}
  \item{existing}{optional logical value specifying whether some of the relevant models
have already been fitted}
  \item{penaliz}{optional character value to determine how models are ranked (see Details)}
  \item{pn.options}{required character string specifying name of 
\code{\link{list}} object populated with starting 
parameter estimates, fitting options and bounds}

}
\description{This function performs model selection for \code{\link{nlsList}} models fitted using
 \code{\link{SSposnegRichards}}.}
\details{First, whether parameter M should be fixed
(see \code{\link{SSposnegRichards}}) is determined by fitting models 12 and 20 and comparing
their perfomance using \code{\link{extraF}}.
If model 12 provides superior performance (variable values of  M) then 16 models that estimate M
are run
(models 1 through 16), otherwise the models with fixed M are fitted (models 21 through 36).
Fitting these \code{\link{nlsList}} models can be time-consuming (2-4 hours using the dataset
\code{\link{posneg.data}} that encompasses 100 individuals) and if several of the relevant
models are already fitted the option existing=TRUE can be used to avoid refitting models that
already exist globally (note that a model object in which no grouping levels were successfully
parameterized will be refitted, as will objects that are not of class nlsList).

Specifying forcemod=3 will force model selection to only consider fixed M models and setting
forcemod=4 will force model selection to consider models with varying values of M only.
If fitting both models
12 and 20 fails, fixed M models will be used by default.

Models are ranked by modified pooled residual square error. By default residual standard error
is divided by the square root of sample size. This exponentially penalizes models for which very few
grouping levels (individuals) are successfully parameterized (the few individuals that are
parameterized in these models are fit unsuprisingly well) using a function based on the relationship
between standard error and sample size. However, different users may have different preferences
and these can be specified in the argument penaliz (which residual
standard error is multiplied by). This argument must be a character value
that contains the character n (sample size) and must be a valid right hand side (RHS) of a formula:
e.g. 1*(n), (n)^2. It cannot contain more than one n but could be a custom function, e.g. FUN(n).}
\value{A list object with two components: $'Model rank table' contains the
statistics from \code{\link{extraF}} ranked by the  modified residual standard error,
and $'P values from pairwise extraF comparison' is a matrix of P values from
\code{\link{extraF}} for legitimate comparisons (nested and successfully fitted models).
The naming convention for models is a concatenation of 'richardsR', the modno and '.lis'
which is shortened in the matrix output, where the number of parameters has been
pasted in parentheses to allow users to easily distinguish the more general model from
the more reduced model
(see \code{\link{extraF}} and \code{\link{SSposnegRichards}}).}

\author{Stephen Oswald <ozsao23@hotmail.com>}
\note{If appropriate bounds (or starting parameters) are not available in the list specified by the variable supplied 
to \code{pn.options}, \code{\link{modpar}} will be called automatically prior to model selection. 
During selection, text is output to the screen to inform the user of the progress of model selection
(which model is being fitted, which were fit successfully)}


\seealso{\code{\link{extraF}}
\code{\link{SSposnegRichards}}
\code{\link{nlsList}}}
\examples{
#these examples will take a long while to run as they have to complete the 32 model comparison
#run model selection for posneg.data object (only first 3 group levels for example's sake)
\donttest{
   data(posneg.data)
   subdata <- subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
   modseltable <- pn.mod.compare(subdata$age, subdata$mass,
      subdata$id, existing = FALSE, pn.options = "myoptions")}
    
#fit nlsList model initially and then run model selection
#for posneg.data object when at least one model is already fit
# note forcemod is set to 3 so that models 21-36 are evaluated
\donttest{
    richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
      Infl = Infl, RAsym = RAsym, Rk = Rk, Ri = Ri , modno = 22)
                        ,data = posneg.data, pn.options = "myoptions")
   modseltable <- pn.mod.compare(subdata$age, subdata$mass,
      subdata$id, forcemod = 3, existing = TRUE, pn.options = "myoptions")}
 
#run model selection ranked by residual standard error*sample size
\donttest{
    modseltable <- pn.mod.compare(subdata$age, subdata$mass,
      subdata$id, penaliz='1*(n)', existing = TRUE, pn.options = "myoptions")}
}