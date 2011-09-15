modpar <-
structure(function # Estimate Values to be Used for Fixed FlexParamCurve Parameters
                    (x,
                     ### a numeric vector of primary predictor variable
                     y,
                     ### a numeric vector of response variable
                     first_y = NA,
                     ###  the value of y at minimum x  when it is required to be constrained
                     last_y = NA,
                     ###  the value of y at maximum x when it is required to be constrained
                     x_at_last_y = NA,
                     ###  the final value of x - this is option is currently disabled
                     force8par = FALSE
                     ### logical specifying whether parameters of the negative Richards
                     ### curve should be set to defaults if they cannot be estimated
                     ) {
    ##description<< This function creates the object \eqn{pnmodelparams}
    ## which holds estimates of values for all 8 FlexParamCurve
    ## parameters used for fitting and solving positive-negative Richards curves with
    ## \code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
    ## respectively. 
    options(warn = -1)
    skel <- rep(list(1), 11)
    initval <- c(rep(NA, 11))
    initval <- relist(initval, skel)
    names(initval) <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
        "Ri", "RM", "first_y", "last_y", "x_at_last_y")
    assign("pnmodelparams", initval, envir = globalenv())
    skel1 <- rep(list(1), 16)
    initval1 <- c(rep(NA, 16))
    initval1 <- relist(initval1, skel1)
    names(initval1) <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
        "Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
        "Rimin", "Rimax", "RMmin", "RMmax")
    assign("pnmodelparamsbounds", initval1, envir = globalenv())
    xy <- data.frame(x, y)
    assign("pnnlsxy", xy, envir = globalenv())
    xy <- na.omit(xy)
    value <- NA
    print("try 8 parameter curve fit in nls")
    try(value <- coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
        K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri,
        RM = RM, modno = 18), data = xy)), silent = TRUE)
    if (is.na(value[1]) == TRUE) {
        print("8 parameter nls fit failed, use getInitial to retrieve 8 parameters")
        try(value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
            K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk,
            Ri = Ri, RM = RM, modno = 18), data = xy), silent = TRUE)
    } else {
        print("8 parameter nls fit successful")
    }
    if (is.na(value[1]) == TRUE) {
        print("getInitial failed, use 4 parameter nls fit - no estimates available for recession")
        print("if force8par==TRUE recession parameters estimated as RAsym=Asym*0.05, Rk=K, Ri=Infl, RM=M")
        try({
            value <- coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
                K = K, Infl = Infl, M = M, RAsym = 1, Rk = 1,
                Ri = 1, RM = 1, modno = 19), data = xy))
            if (force8par == TRUE) {
                value <- c(value, value[1] * 0.05, value[2],
                  value[3], value[4])
                names(value) <- c("Asym", "K", "Infl", "M", "RAsym",
                  "Rk", "Ri", "RM")
            }
        }, silent = FALSE)
        if (class(value) == "try-error")
            stop("estimates not available for data provided. Please check data or provide estimates manually, see ?modpar")
    } else {
        print("8 parameter getInitial successful")
    }
    if (length(value) == 4) {
        value <- c(value, rep(NA, 4))
        names(value) <- c(names(value[1:4]), "RAsym", "Rk", "Ri",
            "RM")
    }
    optvals <- c(first_y, last_y, x_at_last_y)
    names(optvals) <- c("first_y", "last_y", "x_at_last_y")
    value <- c(unlist(value), optvals)
    skel <- rep(list(1), 11)
    value1 <- relist(value, skel)
    names(value1) <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
        "Ri", "RM", "first_y", "last_y", "x_at_last_y")
    assign("pnmodelparams", value1, envir = globalenv())
    ##details<< This function creates the object \eqn{pnmodelparams}
    ## which holds estimates of values for all 8 FlexParamCurve
    ## parameters used for fitting and solving positive-negative Richards curves with
    ## \code{\link{SSposnegRichards}} and \code{\link{posnegRichards_eqn}},
    ## respectively. Running this function concurrently creates the object
    ## \eqn{pnmodelparamsbounds} which holds the maximum and minimum values
    ## for parameters to be used by \code{\link{optim}} and \code{\link{nls}}
    ## during parameter estimation. For definitions of parameters see either
    ## \code{\link{SSposnegRichards}} or \code{\link{posnegRichards_eqn}}. These
    ## objects are updated globally and
    ## do not need to be assigned by the user: the output of a return value simply
    ## offers a way to save values elsewhere for posterity.
    ##
    ## Estimates are produced by fitting positive-negative Richards curves in
    ## \code{\link{nls}} using
    ## \code{\link{SSposnegRichards}} for the full 8 parameter model (R1).
    ## If this fails, the function \code{\link{getInitial}} is called to
    ## attempt to produce initial estimates using the same 8 parameter model.
    ## If this also fails, estimates are attempted in the same way using the
    ## 4 parameter (positive only) model (R12). In this case, only the positive
    ## parameters are returned (NAs are substituted for negative parameters)
    ## unless argument force8mod=TRUE, in which case negative parameters are
    ## defaulted to: RAsym = 0.05*Asym, Rk = K, Ri = Infl, RM = M.
    ##
    ## Parameter bounds estimated here for use in \code{\link{optim}} and \code{\link{nls}}
    ## fits within \code{\link{SSposnegRichards}} are
    ## applicable to a wide range of curves, although user may
    ## change these manually in \code{\link{list}} object \eqn{pnmodelparamsbounds}
    ## Bounds are estimated by \code{\link{modpar}} by adding or subtracting  multiples
    ## of fixed parameter values to estimated mean parameter values:
    ## -Asym*0.5 and +Asym*2.5,
    ## -K*0.5 and +K*0.5,
    ## -Infl*2.5 and +Infl*10
    ## -M*2 and +M*2
    ## -RAsym*0.5 and +RAsym*2.5,
    ## -Rk*0.5 and +Rk*0.5,
    ## -Ri*2.5 and +Ri*5
    ## -RM*2 and +RM*2.
    ##
    ## Use force8mod=TRUE if initial call to \code{\link{modpar}} produces estimates for
    ## only 4 parameters and yet an 8 parameter model is desired for \code{\link{SSposnegRichards}}
    ## or \code{\link{posnegRichards_eqn}}.
    ##
    ## When specified, first_y and last_y are saved in \eqn{pnmodelparams} to instruct
    ## \code{\link{SSposnegRichards}} to add this as the first or last value of the response, respectively,
    ## during estimation.
    Amax = pnmodelparams$Asym + (abs(pnmodelparams$Asym) * 2.5)
    Amin = pnmodelparams$Asym - (abs(pnmodelparams$Asym) * 0.5)
    Kmax = pnmodelparams$K + (abs(pnmodelparams$K) * 0.5)
    Kmin = pnmodelparams$K - (abs(pnmodelparams$K) * 0.5)
    Imax = pnmodelparams$Infl + (abs(pnmodelparams$Infl) * 10)
    Imin = pnmodelparams$Infl + (abs(pnmodelparams$Infl) * -2.5)
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
    Rimax = pnmodelparams$Ri + (abs(pnmodelparams$Ri) * 5)
    Rimin = pnmodelparams$Ri + (abs(pnmodelparams$Ri) * -2.5)
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
    options(warn = 0)
    return(value1)
    ##value<< a \code{\link{list}} of estimated fixed values for all
    ## above arguments
}
, ex = function(){
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
)
