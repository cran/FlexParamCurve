modpar <-
structure(function
		(x,
                     y,
                     first_y = NA,
                     x_at_first_y = NA,
                     last_y = NA,
                     x_at_last_y = NA,
                     twocomponent_age = NA,
                     verbose = FALSE,
                     force8par = FALSE,
                     force4par = FALSE
                     ) {
    options(warn = -1)
    if(!is.na(twocomponent_age) & force4par == TRUE) 
    	stop("Cannot force a two component Richards model to have a single component.\nSet force4par to FALSE")
    detl <- TRUE
    if(verbose == TRUE) detl <- FALSE
    skel <- rep(list(1), 15)
    initval <- c(rep(NA, 15))
    initval <- relist(initval, skel)
    names(initval) <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
        "Ri", "RM", "first_y", "x_at_first_y", "last_y", "x_at_last_y",
        "twocomponent_age","verbose","force4par")
    initval$first_y <- first_y
    initval$x_at_first_y <- x_at_first_y
    initval$last_y <- last_y
    initval$x_at_last_y <- x_at_last_y
    initval$twocomponent_age <- twocomponent_age
    initval$verbose <- verbose
    initval$force4par <- force4par
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
    evlfit<-function(val1,val2){
 	richards <- function(x, Asym, K, Infl, M) Asym/Re(as.complex(1 +
    	    M * exp(-K * (x - Infl)))^(1/M))
	 SSposnegRichardsF <- function(x, Asym, K, Infl, M, RAsym,
  	      Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x -
  	      Infl)))^(1/M))) + (RAsym/Re(as.complex(1 + RM * exp(-Rk *
  	      (x - Ri)))^(1/RM)))
	if(is.na(val1[5])){
		evl1<- sum((y-richards(x,as.numeric(val1[1]),as.numeric(val1[2]),
               	 as.numeric(val1[3]),as.numeric(val1[4])))^2)
		evl2<- sum((y-richards(x,as.numeric(val2[1]),as.numeric(val2[2]),
               	 as.numeric(val2[3]),as.numeric(val2[4])))^2)

	}else{
		evl1<- sum((y-SSposnegRichardsF(x,as.numeric(val1[1]),as.numeric(val1[2]),
               	 as.numeric(val1[3]),as.numeric(val1[4]),
		as.numeric(val1[5]),as.numeric(val1[6]),
		as.numeric(val1[6]),as.numeric(val1[8])))^2)
		evl2<- sum((y-SSposnegRichardsF(x,as.numeric(val2[1]),as.numeric(val2[2]),
               	 as.numeric(val2[3]),as.numeric(val2[4]),
		as.numeric(val2[5]),as.numeric(val2[6]),
		as.numeric(val2[6]),as.numeric(val2[8])))^2)
	}
	if(evl1<=evl2) {
	valfin<-val1
	}else{
	valfin<-val2
	}
	return(valfin)
	}
    value <- NA
    succ <- FALSE
    if(force4par == TRUE & is.na(twocomponent_age)) {
      try(value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
    	             K = K, Infl = Infl, M = M, modno = 19), data = xy), silent = detl)
            savvalue<-value
            value<-NA
    print("try 4 parameter curve fit in nls")
    try({value <- coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
                    K = K, Infl = Infl, M = M, modno = 19), data = xy))
         }, silent = detl)
         if(is.na(value[1]) == FALSE) print("4 parameter nls fit successful")
         if(is.na(value[1]) == TRUE) {
          print("4 parameter nls fit failed, use getInitial to retrieve 4 parameters")
         try(value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
	             K = K, Infl = Infl, M = M, modno = 19), data = xy), silent = detl)
            if(is.na(value[1]) == TRUE) {
            stop("estimates not available for data provided. Please check data or provide estimates manually, see ?modpar")
         					} else {
         					print("4 parameter getInitial successful")
         					}
         }
    if(!is.na(value[1]) & !is.na(savvalue[1])) value<-evlfit(savvalue,value)
    }else{
    try(value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
                K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk,
            Ri = Ri, RM = RM, modno = 18), data = xy), silent = detl)
            tst<-get("pnmodelparamsbounds", envir = .GlobalEnv)
     if (is.na(value[1]) == TRUE){
        try(value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
         	             K = K, Infl = Infl, M = M, modno = 19), data = xy), silent = detl)
            bndsvals<-get("pnmodelparamsbounds", envir = .GlobalEnv)
            bndsvals[9:16]<-bndsvals[1:8]
            assign("pnmodelparamsbounds",bndsvals, envir = .GlobalEnv)
            			} else {
            assign("pnmodelparamsbounds",tst, envir = .GlobalEnv)			
            			}
         savvalue<-value
         value<-NA
    if(force4par == TRUE & !is.na(twocomponent_age)) print("Cannot force a two component model to have 4 parameters")
    print("try 8 parameter curve fit in nls")
    try(value <- coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
        K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri,
        RM = RM, modno = 18), data = xy)), silent = detl)
    if (is.na(value[1]) == TRUE) {
        print("8 parameter nls fit failed, use getInitial to retrieve 8 parameters")
        try(value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
            K = K, Infl = Infl, M = M, RAsym = RAsym, Rk = Rk,
            Ri = Ri, RM = RM, modno = 18), data = xy), silent = detl)
             if (is.na(value[1]) == FALSE) {
             		succ <- TRUE
             		print("8 parameter getInitial successful")

             		}
    } else {
        print("8 parameter nls fit successful")
        succ <- TRUE
    }
    if(!is.na(value[1]) & !is.na(savvalue[1])) value<-evlfit(savvalue,value)
    if (is.na(value[1]) == TRUE & is.na(twocomponent_age)) {
        print("getInitial failed, use 4 parameter nls fit - no estimates available for recession")
        print("if force8par==TRUE recession parameters estimated as RAsym=Asym*0.05, Rk=K, Ri=Infl, RM=M")
        try({
            value <- coef(nls(y ~ SSposnegRichards(x, Asym = Asym,
                K = K, Infl = Infl, M = M, modno = 19), data = xy))
            if (force8par == TRUE) {
                value <- c(value, value[1] * 0.05, value[2],
                  value[3], value[4])
                names(value) <- c("Asym", "K", "Infl", "M", "RAsym",
                  "Rk", "Ri", "RM")
            }
        }, silent = detl)
        if(is.na(value[1]) == TRUE) {
          print("4 parameter nls fit failed, use getInitial to retrieve 4 parameters")
          try({
 	    value <- getInitial(y ~ SSposnegRichards(x, Asym = Asym,
	             K = K, Infl = Infl, M = M, modno = 19), data = xy)
	    if (force8par == TRUE) {
	        value <- c(value, value[1] * 0.05, value[2],
	          value[3], value[4])
	        names(value) <- c("Asym", "K", "Infl", "M", "RAsym",
	        "Rk", "Ri", "RM")
	     }
          }, silent = detl)
          if(is.na(value[1]) == TRUE) print("4 parameter getInitial failed")
        } else {
        print("4 parameter nls successful")
        }
        if(is.na(value[1]) == TRUE)
            stop("estimates not available for data provided. Please check data or provide estimates manually, see ?modpar")
    } else {
         if (succ != TRUE)
                  stop("estimates not available for data provided. Please check data or provide estimates manually, see ?modpar")
    }
    }
    if (length(value) == 4) {
        value <- c(value, rep(NA, 4))
        names(value) <- c(names(value[1:4]), "RAsym", "Rk", "Ri",
            "RM")
    }
    optvals <- c(first_y, x_at_first_y, last_y, x_at_last_y, twocomponent_age, verbose, force4par)
    names(optvals) <- c("first_y", "x_at_first_y", "last_y", "x_at_last_y", "twocomponent_age", "verbose", "force4par")
    value <- c(unlist(value), optvals)
    skel <- rep(list(1), 15)
    value1 <- relist(value, skel)
    names(value1) <- c("Asym", "K", "Infl", "M", "RAsym", "Rk",
        "Ri", "RM", "first_y", "x_at_first_y", "last_y", "x_at_last_y",
        "twocomponent_age", "verbose", "force4par")
    lodpar <- get("pnmodelparams", envir = .GlobalEnv)
     value1$twocomponent_age <- lodpar$twocomponent_age
    value1$verbose <- initval$verbose
    value1$force4par <- initval$force4par
    assign("pnmodelparams", value1, envir = globalenv())
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
    }
, ex = function(){
        data(posneg_data)
        modpar(posneg_data$age,posneg_data$mass)

        modpar(posneg_data$age,posneg_data$mass)
        subdata<-subset(posneg_data, as.numeric(row.names (posneg_data) ) < 53)
        richardsR1.lis<-nlsList(mass~SSposnegRichards(age,Asym=Asym,K=K,
        	Infl=Infl,M=M,RAsym=RAsym,Rk=Rk,Ri=Ri,RM=RM,modno=1)
                        ,data=subdata)

})
