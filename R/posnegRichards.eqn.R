posnegRichards.eqn <-
structure(function(x, Asym = NA,
    K = NA, Infl = NA, M = NA, RAsym = NA, Rk = NA, Ri = NA, RM = NA,
    modno, pn.options) {
    params<-list(Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym,
        Rk = Rk, Ri = Ri, RM = RM, first.y = NA, x.at.first.y = NA, 
        last.y = NA, x.at.last.y = NA, twocomponent.x = NA, 
        verbose = FALSE, force4par = FALSE)
    pnmodelparams <- rep(list(NA),15)
    pnmodelparams[14]<-FALSE
    pnmodelparams[15]<-FALSE
    pnopt<- get(as.character( pn.options[1] ), .GlobalEnv)
    names(pnmodelparams)<-names(params)
    pnmodelparams[ names(pnmodelparams) %in% names(pnopt) ] <- pnopt[ names(pnopt) %in% names(pnmodelparams) ]
    if (is.na(Asym[1])) stop(paste("Parameter Asym required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    if (is.na(Infl[1])) stop(paste("Parameter Infl required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    if (is.na(K[1]) & modno != 17) stop(paste("Parameter K required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    fractM <- M
    fractRM <- RM
    if (modno > 19) {
        fractM <- 1/pnmodelparams$M
        M <- pnmodelparams$M
    } else {
        if (is.na(M[1])) stop(paste("Parameter M required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
    }
    if (modno != 17) {
    if (modno == 2 | modno == 22 | modno == 10 | modno == 30 |
        modno == 11 | modno == 31 | modno == 12 | modno == 32 |
        modno == 19 | modno == 13 | modno == 33 | modno == 14 |
        modno == 34 | modno == 15 | modno == 35 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        fractRM <- 1/pnmodelparams$RM
        RM <- pnmodelparams$RM
    } else {
    if (is.na(RM[1])) stop(paste("Parameter RM required for modno = ", modno,
    	" but is absent from user-provided call", sep = ""))
    }
    if (modno == 3 | modno == 23 | modno == 5 | modno == 25 |
        modno == 7 | modno == 27 | modno == 9 | modno == 29 |
        modno == 10 | modno == 30 | modno == 12 | modno == 32 |
        modno == 19 | modno == 14 | modno == 34 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        RAsym <- pnmodelparams$RAsym
    } else {
        if (is.na(RAsym[1])) stop(paste("Parameter RAsym required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
    }
    if (modno == 3 | modno == 23 | modno == 4 | modno == 24 |
        modno == 5 | modno == 25 | modno == 6 | modno == 26 |
        modno == 10 | modno == 30 | modno == 11 | modno == 31 |
        modno == 12 | modno == 32 | modno == 19 | modno == 13 |
        modno == 33 | modno == 20 | modno == 32) {
        Rk <- pnmodelparams$Rk
    } else {
        if (is.na(Rk[1])) stop(paste("Parameter Rk required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
    }
    if (modno == 4 | modno == 24 | modno == 5 | modno == 25 |
        modno == 8 | modno == 28 | modno == 9 | modno == 29 |
        modno == 11 | modno == 31 | modno == 12 | modno == 32 |
        modno == 19 | modno == 15 | modno == 35 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        Ri <- pnmodelparams$Ri
    } else {
        if (is.na(Ri[1])) stop(paste("Parameter Ri required for modno = ", modno,
        	" but is absent from user-provided call", sep = ""))
    }
    try({
    if (Re(as.complex(1 + M[1] * exp(-K[1] * (max(x) - Infl[1])))) < 0) {
     	    fractM <- round(1/M)
    } else {
            fractM <- 1/M
    }
    if (Re(as.complex(1 + RM[1] * exp(-Rk[1] * (max(x) - Ri[1])))) < 0) {
     	    fractRM <- round(1/RM)
    } else {
       	    fractRM <- 1/RM
    }
    },silent = TRUE)
    }
    force4par<-pnmodelparams$force4par
    if(modno == 12 | modno == 32) {
    force4par = TRUE
    print("note: the model selected does not have any second curve parameters (i.e. is monotonic)")
    					}
    if (!is.na(pnmodelparams$twocomponent.x)) {
        if (force4par == TRUE)
            stop("Cannot force a two component Richards model to have a single component.\nSet force4par to FALSE")
        c((Asym/Re(as.complex(1 + M * exp(-K * ((x[x <= pnmodelparams$twocomponent.x]) -
            Infl)))^(fractM))), (RAsym/Re(as.complex(1 + RM *
            exp(-Rk * ((x[x > pnmodelparams$twocomponent.x]) -
                Ri)))^(fractRM))))
    } else {
        if (modno == 17) {
            if (force4par == TRUE) {
                Asym/Re(as.complex(1 + exp(Infl - x)/M))
            } else {
                (Asym/Re(as.complex(1 + exp(Infl - x)/M))) +
                  (RAsym/Re(as.complex(1 + exp(Ri - x)/RM)))
            }
        } else {
            if (force4par == TRUE) {
                (Asym/Re(as.complex(1 + M * exp(-K * ((x) - Infl)))^(fractM)))
            } else {
                (Asym/Re(as.complex(1 + M * exp(-K * ((x) - Infl)))^(fractM))) +
                  (RAsym/Re(as.complex(1 + RM * exp(-Rk * ((x) -
                    Ri)))^(fractRM)))
            }
        }
    }
}, ex = function(){
    require(graphics)
    data(posneg.data)
    modpar(posneg.data$age, posneg.data$mass)
    y <- posnegRichards.eqn(10, 1000, 0.5, 25, 1, 100, 0.5, 125, 1, modno = 1)

    y <- posnegRichards.eqn(10 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12)

    plot(1:200 ,posnegRichards.eqn(1:200 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12),xlim=c(1, 200),
    xlab = "x", ylab = "y",pch = 1, cex = 0.7)
    }
)
