change.pnparameters <-
structure(function 
                               (Asym = NA,
                               K = NA,
                                Infl = NA,
                                M = NA,
                                RAsym = NA,
                               Rk = NA,
                                Ri = NA,
                                RM = NA,
                                first_y = NA,
                                x_at_first_y = NA,
                                last_y = NA,
                                x_at_last_y = NA,
				twocomponent_age = NA,
				verbose = NA,
				force4par = NA
				) {
    newparams <- list(Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym,
        Rk = Rk, Ri = Ri, RM = RM, first_y = first_y, x_at_first_y = x_at_first_y, 
        last_y = last_y,
        x_at_last_y = x_at_last_y, twocomponent_age = twocomponent_age, 
        verbose = verbose, force4par = force4par)
    adjmodelparams <- get("pnmodelparams", envir = .GlobalEnv)
    for (i in 1:15) {
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
}
, ex = function(){
data(posneg_data)
modpar(posneg_data$age,posneg_data$mass)
change.pnparameters(Asym=10000,Infl=80,M=5,RAsym=10000,Ri=240,RM=5)

change.pnparameters(M=1,RM=0.5,first_y=45.5)
}
)
