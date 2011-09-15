change.pnparameters <-
structure(function # Change Fixed Parameter Values
### Function to alter values of parameters to be used by \code{\link{SSposnegRichards}}
### or \code{\link{posnegRichards_eqn}}
### as the fixed values in equations where parameters are fixed
                               (Asym = NA,
                                ### a numeric value for the asymptote of the positive (increasing) curve
                                K = NA,
                                ### a numeric value for the rate parameter of the positive (increasing) curve
                                Infl = NA,
                                ### a numeric value for the point of inflection of the positive (increasing) curve
                                M = NA,
                                ### a numeric value for the shape parameter of the positive (increasing) curve
                                RAsym = NA,
                                ### a numeric value for  the asymptote of the negative (decreasing) curve
                                RK = NA,
                                ### a numeric value for the rate parameter of the negative (decreasing) curve
                                Ri = NA,
                                ### a numeric value for the point of inflection of the negative (decreasing) curve
                                RM = NA,
                                ### a numeric value for the shape parameter of the negative (decreasing) curve
                                first_y = NA,
                                ### the value of y at minimum x when it is required to be constrained
                                last_y = NA,
                                ### the value of y at maximum x when it is required to be constrained
                                x_at_last_y = NA
### the final value of x - this is option is currently disabled
                                ) {
    newparams <- c(Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym,
        Rk = RK, Ri = Ri, RM = RM, first_y = first_y, last_y = last_y,
        x_at_last_y = x_at_last_y)
    adjmodelparams <- get("pnmodelparams", envir = .GlobalEnv)
    for (i in 1:11) {
        if (is.na(newparams[i]) == FALSE)
            adjmodelparams[i] <- newparams[i]
    }
    names(adjmodelparams) <- names(newparams)
    assign("pnmodelparams", adjmodelparams, envir = globalenv())
    Amax = pnmodelparams$Asym + (abs(pnmodelparams$Asym) * 2.5)
    Amin = pnmodelparams$Asym - (abs(pnmodelparams$Asym) * 0.5)
    Kmax = pnmodelparams$K + (abs(pnmodelparams$K) * 0.5)
    Kmin = pnmodelparams$K - (abs(pnmodelparams$K) * 0.5)
    Imax = pnmodelparams$Infl + (abs(pnmodelparams$Infl) * 2.5)
    Imin = pnmodelparams$Infl - (abs(pnmodelparams$Infl) * 1.5)
    while (abs(Imax * Kmax) > 700) Imax = Imax * 0.9
    while (abs(Imin * Kmax) > 700) Imin = Imin * 0.9
    Mmax = pnmodelparams$M + abs(pnmodelparams$M * 2)
    Mmin = pnmodelparams$M - abs(pnmodelparams$M * 2)
    if (is.na(pnmodelparams$RAsym)) {
        pnmodelparams$RAsym <- pnmodelparams$Asym
        pnmodelparams$Rk <- pnmodelparams$K
        pnmodelparams$Ri <- pnmodelparams$Infl
        pnmodelparams$RM <- pnmodelparams$M
    }
    RAmax = pnmodelparams$RAsym + (abs(pnmodelparams$RAsym) *
        2.5)
    RAmin = pnmodelparams$RAsym - (abs(pnmodelparams$RAsym) *
        0.5)
    Rkmax = pnmodelparams$Rk + (abs(pnmodelparams$Rk) * 0.5)
    Rkmin = pnmodelparams$Rk - (abs(pnmodelparams$Rk) * 0.5)
    Rimax = pnmodelparams$Ri + (abs(pnmodelparams$Ri) * 1.25)
    Rimin = pnmodelparams$Ri - (abs(pnmodelparams$Ri) * 0.5)
    while (abs(Rimax * Rkmax) > 700) Rimax = Rimax * 0.9
    while (abs(Rimin * Rkmax) > 700) Rimin = Rimin * 0.9
    RMmax = pnmodelparams$RM + abs(pnmodelparams$RM * 2)
    RMmin = pnmodelparams$RM - abs(pnmodelparams$RM * 2)
    value2 <- c(Amin, Amax, Kmin, Kmax, Imin, Imax, Mmin, Mmax,
        RAmin, RAmax, Rkmin, Rkmax, Rimin, Rimax, RMmin, RMmax)
    skel1 <- rep(list(1), 16)
    value3 <- relist(value2, skel1)
    names(value3) <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
        "Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
        "Rimin", "Rimax", "RMmin", "RMmax")
    assign("pnmodelparamsbounds", value3, envir = globalenv())
    return(adjmodelparams)
    ##value<< a \code{\link{list}} of values for all above arguments,
    ## with new values substituted where specified in the call
    ##details<< This function provides a simple way for the user to update
    ## the object \eqn{pnmodelparams} that holds fixed values and options
    ## for fitting and solving positive-negative Richards curves with
    ## \code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
    ## respectively. Running this function concurrently updates the object
    ## \eqn{pnmodelparamsbounds} which holds the maximum and minimum values
    ## for parameters to be used by \code{\link{optim}} and \code{\link{nls}}
    ##during parameter estimation
    ## in \code{\link{SSposnegRichards}}. These objects are updated globally and
    ## do not need to be assigned by the user: the output of a return value simply
    ## offers a way to save values elsewhere for posterity.
    ##
    ## Both objects must exist before this function is called. Use \code{\link{modpar}}
    ## to estimate values for all parameters and create \eqn{pnmodelparams} and
    ##  \eqn{pnmodelparamsbounds}. See  \code{\link{modpar}} for details of bounding.
    ##seealso<< \code{\link{modpar}} \code{\link{SSposnegRichards}} \code{\link{posnegRichards_eqn}}
    ##note<< Requires \code{\link{modpar}} to be have been run prior to execution
}
, ex = function(){
# change all fixed values except K and Rk
data(posneg_data)
modpar(posneg_data$age,posneg_data$mass)
change.pnparameters(Asym=10000,Infl=80,M=5,RAsym=10000,Ri=240,RM=5)

# change fixed values of M and constrain hatching mass to 45.5 in a growth curve
change.pnparameters(M=1,RM=0.5,first_y=45.5)
}
)
