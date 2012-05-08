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
                                first.y = NA,
                                x.at.first.y = NA,
                                last.y = NA,
                                x.at.last.y = NA,
				twocomponent.x = NA,
				verbose = NA,
				force4par = NA,
				pn.options
				) {
    newparams <- list(Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym,
        Rk = Rk, Ri = Ri, RM = RM, first.y = first.y, x.at.first.y = x.at.first.y, 
        last.y = last.y,
        x.at.last.y = x.at.last.y, twocomponent.x = twocomponent.x, 
        verbose = verbose, force4par = force4par)
    pnoptnm<- as.character( pn.options[1] ) 
    adjmodelparams <- get(pnoptnm, envir = .GlobalEnv)
     for (i in 1:15) {
        if (is.na(newparams[i]) == FALSE)
            adjmodelparams[names(adjmodelparams) %in% names(newparams[i])] <- newparams[i]
    }
    adjmodelparams <- adjmodelparams[names(adjmodelparams) %in% names(newparams)]
    pnmodelparams <- adjmodelparams 
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
    tmplate <- get(pnoptnm, envir = .GlobalEnv)
    adjmodelparams[names(adjmodelparams) %in% names(value3)] <- value3[names(value3) %in% names(value3)]
    "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    tmplate[names(tmplate) %in% names(adjmodelparams)] <- adjmodelparams[names(adjmodelparams) %in% names(adjmodelparams)]
    assign(pnoptnm, tmplate, .GlobalEnv) 
    return(tmplate)
}
, ex = function(){
data(posneg.data)
modpar(posneg.data$age,posneg.data$mass)
change.pnparameters(Asym=10000,Infl=80,M=5,RAsym=10000,Ri=240,RM=5)

change.pnparameters(M=1,RM=0.5,first.y=45.5)
}
)
