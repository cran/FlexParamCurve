posnegRichardsInit <-
function(mCall, LHS, data, modno) {
    xy <- sortedXyData(mCall[["x"]], LHS, data)
    modno <- mCall[["modno"]]
    modelparams <- NA
    try(modelparams <- get("pnmodelparams", envir = .GlobalEnv),
        silent = TRUE)
    if (is.na(modelparams[1]) == TRUE) {
        if (modno != 18 & modno != 19) {
            stop("ERROR: global parameters empty: run function modpar before using selfStart functions: see ?modpar")
        }
    }
    if (is.na(modelparams$first_y) == FALSE) {
        xy <- rbind(xy, c(0, modelparams$first_y))
        xy <- xy[order(xy$x), ]
    }
    if (is.na(modelparams$last_y) == FALSE) {
        xy <- rbind(xy, c(modelparams$last_y, modelparams$last_y))
    }
    sub <- subset(xy, xy$y == max(xy$y))
    tAsym <- mean(sub$x)
    SSposnegRichardsF <- function(x, Asym, K, Infl, M, RAsym,
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x -
        Infl)))^(1/M))) - (RAsym/Re(as.complex(1 + RM * exp(-Rk *
        (x - Ri)))^(1/RM)))
    SSposnegRichardsFM <- function(x, Asym, K, Infl, M, RAsym,
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x -
        Infl)))^round(1/M))) - (RAsym/Re(as.complex(1 + RM *
        exp(-Rk * (x - Ri)))^(1/RM)))  #integer version for -M
    SSposnegRichardsFRM <- function(x, Asym, K, Infl, M, RAsym,
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x -
        Infl)))^(1/M))) - (RAsym/Re(as.complex(1 + RM * exp(-Rk *
        (x - Ri)))^round(1/RM)))
    SSposnegRichardsFMRM <- function(x, Asym, K, Infl, M, RAsym,
        Rk, Ri, RM) (Asym/Re(as.complex(1 + M * exp(-K * (x -
        Infl)))^round(1/M))) - (RAsym/Re(as.complex(1 + RM *
        exp(-Rk * (x - Ri)))^round(1/RM)))
    SSposnegRichardsF17 <- function(x, Asym, Infl, M, RAsym,
        Ri, RM) (Asym/Re(as.complex(1 + exp(Infl - x)/M))) -
        (RAsym/Re(as.complex(1 + exp(Ri - x)/RM)))
    richards <- function(x, Asym, K, Infl, M) Asym/Re(as.complex(1 +
        M * exp(-K * (x - Infl)))^(1/M))
    richards3 <- function(x, Asym, K, Infl) Asym/Re(as.complex(1 +
        pnmodelparams$M * exp(-K * (x - Infl)))^(1/pnmodelparams$M))
    richards173 <- function(x, Asym, Infl, M) Asym/Re(as.complex(1 +
        exp(Infl - x)/M))
    try(modelparamsbounds <- get("pnmodelparamsbounds", envir = .GlobalEnv),
        silent = TRUE)
    if (is.na(modelparamsbounds[1]) == TRUE) {
        if (modno != 18 & modno != 19) {
            stop("ERROR: fix parameter bounds empty: run function modpar before using selfStart functions: see ?modpar")
        }
    }
    if (modno != 18 & modno != 19) {
        Amax = modelparamsbounds$Amax
        Amin = modelparamsbounds$Amin
        Kmax = modelparamsbounds$Kmax
        Kmin = modelparamsbounds$Kmin
        Imax = modelparamsbounds$Imax
        Imin = modelparamsbounds$Imin
        Mmax = modelparamsbounds$Mmax
        Mmin = modelparamsbounds$Mmin
        RAmax = modelparamsbounds$RAmax
        RAmin = modelparamsbounds$RAmin
        Rkmax = modelparamsbounds$Rkmax
        Rkmin = modelparamsbounds$Rkmin
        Rimax = modelparamsbounds$Rimax
        Rimin = modelparamsbounds$Rimin
        RMmax = modelparamsbounds$RMmax
        RMmin = modelparamsbounds$RMmin
    }
    xyE <- subset(xy, xy$x <= tAsym)
    if (nrow(xyE) < 5) {
        stop("ERROR: too few distinct input values to fit the positive Richards model, aborting")
    } else {
        rng <- range(xyE$y)
        drng <- diff(rng)
        xyE$prop <- (xyE$y - rng[1L] + 0.05 * drng)/(1.1 * drng)
        for (xk in 1:3) {
            if (xk == 1)
                {
                  transy <- (log(xyE$prop/(1 - xyE$prop)))
                  irL <- lm(x ~ transy, data = xyE)
                  ssqL <- sum((irL$residuals)^2)
                }
            if (xk == 2) {
                transy <- (-log(8/(27 * xyE$prop)))
                irV <- lm(x ~ transy, data = xyE)
                ssqV <- sum((irV$residuals)^2)
            }
            if (xk == 3) {
                transy <- (-log(-log(xyE$prop)))
                irG <- lm(x ~ transy, data = xyE)
                ssqG <- sum((irG$residuals)^2)
            }
        }
        if (ssqL < ssqV & ssqL < ssqG) {
            ir <- as.vector(irL$coefficients)
            M <- 1
        } else {
            if (ssqV < ssqL & ssqV < ssqG) {
                ir <- as.vector(irV$coefficients)
                M = -0.3
            } else {
                ir <- as.vector(irG$coefficients)
                M = 0.1
            }
        }
        if (modno > 19) {
            if (M > 0.3) {
                M <- 1
            } else {
                if (M > 0) {
                  M <- 0.1
                } else {
                  M <- (-0.3)
                }
            }
        }
        Infl = ir[1L]
        K = 1/(ir[2L])
        Asym = rng[2L]
        if (modno == 17)
            {
                M <- ((Infl - mean(xyE$x)) * M) - (Infl - mean(xyE$x))
            }
        val <- c(Asym, K, Infl, M)
        if (modno > 19 | modno == 17) {
            if (modno > 19)
                val <- c(Asym, K, Infl)
            if (modno == 17)
                val <- c(Asym, Infl, M)
            func21 <- function(val) {
                func23 <- function(x, K, Infl) {
                  x = x + (0+0i)
                  (1 + M * exp(-K * (x - Infl)))^(1/M)
                }
                func173 <- function(x, Infl, M) {
                  x = x + (0+0i)
                  (1 + exp(Infl - x)/M)
                }
                Asym = val[1]
                if (modno == 17) {
                  K = 1
                  Infl = val[2]
                  M = val[3]
                  P1 <- Re(func173(xyE$x, Infl, M))
                } else {
                  K = val[2]
                  Infl = val[3]
                  P1 <- Re(func23(xyE$x, K, Infl))
                }
                P1[is.na(P1)] <- 1e-290 * pnmodelparams$Asym
                P1[P1 == Inf] <- 1e-290 * pnmodelparams$Asym
                y1 <- Asym/P1
                evl <- sum((xyE$y - y1)^2)
                if (evl == Inf) {
                  evl <- 1e+290
                }
                try(if (min(Im(as.complex(1 + M * exp(-K * ((0:200) -
                  Infl)))^(1/M))) < 0) {
                  evl <- 1e+200
                }, silent = T)
                return(evl)
            }
            oppar = 0
            is.na(oppar) <- TRUE
            func25 <- function(Infl) {
                repl <- Kmax
                if (M < 1e-20) {
                  op <- data.frame((0:200) - Infl)
                  if (nrow(subset(op, op == 0)) > 0)
                    repl <- (-(2.170611e-16 - 1e-06 + (log(abs(M)))/(Infl +
                      1e-10)))
                }
                return(repl)
            }
            if (Kmin < 1e-05)
                Kmin = 1e-05
            if (modno == 17) {
                oppar <- optim(c(Asym, Infl, M), func21, method = "L-BFGS-B",
                  lower = c(Amin, Imin, Mmin), upper = c(Amax,
                    Imax, Mmax), control = list(maxit = 2000,
                    parscale = c(10000, 100, 10)))
            } else {
                oppar <- optim(c(Asym, K, Infl), func21, method = "L-BFGS-B",
                  lower = c(Amin, Kmin, Imin), upper = c(Amax,
                    func25(Infl), Imax), control = list(maxit = 2000,
                    parscale = c(10000, 0.01, 100)))
            }
            if (is.na(oppar[1]) == FALSE) {
                if (oppar$convergence < 52) {
                  val <- (c(oppar$par))
                } else {
                  stop("ERROR: positive optimization failed,aborting")
                }
            } else {
                stop("ERROR: positive optimization failed,aborting")
            }
            if (modno == 17) {
                Asym = val[1]
                Infl = val[2]
                M = val[3]
                value <- c(Asym, Infl, M)
            } else {
                Asym = val[1]
                K = val[2]
                Infl = val[3]
                value <- c(Asym, K, Infl)
            }
            Kmax = K + (abs(K) * 0.25)
            Kmin = K - (abs(K) * 0.25)
            if (Kmax < 1e-05)
                Kmax = 0.05
            if (Kmin < 1e-05)
                Kmin = 1e-05
            Imax = Imax - (abs(Imax - Infl) * 0.75)
            Imin = Imin + (abs(Imin - Infl) * 0.75)
            if (modno == 17) {
            } else {
                while (abs(Imax * Kmax) > 700) Imax = Imax *
                  0.9
                while (abs(Imin * Kmax) > 700) Imin = Imin *
                  0.9
            }
            Mmax = Mmax - (abs(Mmax - M) * 0.75)
            Mmin = Mmin + (abs(Mmin - M) * 0.75)
            pars <- 1
            is.na(pars) <- TRUE
            if (modno == 17) {
                try(pars <- as.vector(coef(nls(y ~ richards173(x,
                  Asym, Infl, M), data = xyE, start = list(Asym = Asym,
                  Infl = Infl, M = M), nls.control(maxiter = 1000,
                  tol = 1), algorithm = "port", lower = list(Asym = Amin,
                  Infl = Imin, M = Mmin), upper = list(Asym = Amax,
                  I = Imax, M = Mmax)))), silent = TRUE)
            } else {
                try(pars <- as.vector(coef(nls(y ~ richards3(x,
                  Asym, K, Infl), data = xyE, start = list(Asym = Asym,
                  K = K, Infl = Infl), nls.control(maxiter = 1000,
                  tol = 1), algorithm = "port", lower = list(Asym = Amin,
                  K = Kmin, Infl = Imin), upper = list(Asym = Amax,
                  K = Kmax, I = Imax)))), silent = TRUE)
            }
            if (is.na(pars[1]) == FALSE) {
                if (modno == 17) {
                  value <- c(pars[1L], K, pars[2L], pars[3L],
                    (Asym * 0.05), K, pars[2L], tAsym)
                } else {
                  value <- c(pars[1L], pars[2L], pars[3L], M,
                    (Asym * 0.05), pars[2L], tAsym, M)
                }
            } else {
                print("Warning: nls failed for positive Richards model, using optim parameters")
                value <- c(Asym, K, Infl, M, (Asym * 0.05), K,
                  tAsym, M)
            }
        } else {
            func1 <- function(val) {
                Asym = val[1]
                K = val[2]
                Infl = val[3]
                M = val[4]
                func3 <- function(x, K, Infl, M) {
                  x = x + (0+0i)
                  (1 + M * exp(-K * (x - Infl)))^(1/M)
                }
                y1 <- Asym/Re(func3(xyE$x, K, Infl, M))
                evl <- sum((xyE$y - y1)^2)
                try(if (min(Im(as.complex(1 + M * exp(-K * ((0:200) -
                  Infl)))^(1/M))) < 0) {
                  evl <- 1e+200
                }, silent = T)
                return(evl)
            }
            oppar = 0
            is.na(oppar) <- TRUE
            func5 <- function(Infl, M) {
                repl <- Kmax
                if (M < 1e-20) {
                  op <- data.frame((0:200) - Infl)
                  if (nrow(subset(op, op == 0)) > 0)
                    repl <- (-(2.170611e-16 - 1e-06 + (log(abs(M)))/(Infl +
                      1e-10)))
                }
                return(repl)
            }
            if (modno == 18 | modno == 19) {
                Amax = Asym + (abs(Asym) * 2.5)
                Amin = Asym - (abs(Asym) * 0.5)
                Kmax = K + (abs(K) * 0.5)
                Kmin = K - (abs(K) * 0.5)
                Imax = Infl + (abs(Infl) * 10)
                Imin = Infl + (abs(Infl) * -2.5)
                while (abs(Imax * Kmax) > 700) Imax = Imax *
                  0.9
                while (abs(Imin * Kmax) > 700) Imin = Imin *
                  0.9
                Mmax = M + abs(M * 2)
                Mmin = M - abs(M * 2)
                RAmax = Amax
                RAmin = Amin
                Rkmax = Kmax
                Rkmin = Kmin
                Rimax = Imax
                Rimin = Imin
                RMmax = Mmax
                RMmin = Mmin
                skel <- rep(list(1), 11)
                exportparams <- c(Asym, K, Infl, M, (Asym * 0.05),
                  K, tAsym, M, NA, NA, NA)
                exportparams <- relist(exportparams, skel)
                names(exportparams) <- c("Asym", "K", "Infl",
                  "M", "RAsym", "Rk", "Ri", "RM", "first_y",
                  "last_y", "x_at_last_y")
                skel <- rep(list(1), 16)
                exportparamsbounds <- c(Amin, Amax, Kmin, Kmax,
                  Imin, Imax, Mmin, Mmax, RAmin, RAmax, Rkmin,
                  Rkmax, Rimin, Rimax, RMmin, RMmax)
                exportparamsbounds <- relist(exportparamsbounds,
                  skel)
                names(exportparamsbounds) <- c("Amin", "Amax",
                  "Kmin", "Kmax", "Imin", "Imax", "Mmin", "Mmax",
                  "RAmin", "RAmax", "Rkmin", "Rkmax", "Rimin",
                  "Rimax", "RMmin", "RMmax")
                assign("pnmodelparams", exportparams, envir = globalenv())
                assign("pnmodelparamsbounds", exportparamsbounds,
                  envir = globalenv())
            }
            if (Kmin < 1e-05)
                Kmin = 1e-05
            oppar <- optim(c(Asym, K, Infl, M), func1, method = "L-BFGS-B",
                lower = c(Amin, Kmin, Imin, Mmin), upper = c(Amax,
                  func5(Infl, M), Imax, Mmax), control = list(maxit = 2000,
                  parscale = c(10000, 0.01, 100, 10)))
            if (is.na(oppar[1]) == FALSE) {
                if (oppar$convergence < 52) {
                  val <- (c(oppar$par))
                } else {
                  stop("ERROR: positive optimization failed,aborting")
                }
            } else {
                stop("ERROR: positive optimization failed,aborting")
            }
            Asym = val[1]
            K = val[2]
            Infl = val[3]
            M = val[4]
            pars = 0
            is.na(pars) <- TRUE
            Amax = Asym + (abs(Asym) * 2.5)
            Amin = Asym - (abs(Asym) * 0.5)
            Kmax = K + (abs(K) * 0.25)
            Kmin = K - (abs(K) * 0.25)
            if (Kmax < 1e-05)
                Kmax = 0.05
            if (Kmin < 1e-05)
                Kmin = 1e-05
            Imax = Infl + (abs(Infl) * 5)
            Imin = Infl + (abs(Infl) * -1.5)
            while (abs(Imax * Kmax) > 700) Imax = Imax * 0.9
            while (abs(Imin * Kmax) > 700) Imin = Imin * 0.9
            Mmax = M + abs(M * 0.75)
            Mmin = M - abs(M * 0.75)
            try(pars <- as.vector(coef(nls(y ~ richards(x, Asym,
                K, Infl, M), data = xyE, start = list(Asym = Asym,
                K = K, Infl = Infl, M = M), nls.control(maxiter = 1000,
                tol = 1), algorithm = "port", lower = list(Asym = Amin,
                K = Kmin, Infl = Imin, M = Mmin), upper = list(Asym = Amax,
                K = Kmax, I = Imax, M = Mmax)))), silent = TRUE)
            if (is.na(pars[1]) == FALSE) {
                value <- c(pars[1L], pars[2L], pars[3L], pars[4L],
                  (Asym * 0.05), pars[2L], tAsym, pars[4L])
            } else {
                print("Warning: nls failed for positive Richards model, using optim parameters")
                value <- c(Asym, K, Infl, M, (Asym * 0.05), K,
                  tAsym, M)
            }
        }
    }
    xyL <- subset(xy, xy$x >= tAsym)
    if (nrow(xyL) < 3 | modno == 12 | modno == 20 | modno ==
        32 | modno == 19) {
        if (nrow(xyL) < 3) {
            print("Warning: too few distinct input values to fit the negative Richards model, defaulting to RAsym= 5% of peak mass, Rk=k, Ri=age at peak, RM=M")
        }
        value <- c(value[1L], value[2L], value[3L], value[4L],
            (Asym * 0.05), value[2L], tAsym, value[4L])
        inputval <- 0
        is.na(inputval) <- TRUE
        if (modno == 1 | modno == 21 | modno == 18) {
            val1 <- data.frame(value[5L], value[6L], value[7L],
                value[8L])
            names(val1) <- c("RAsym", "Rk", "Ri", "RM")
            inputmin <- c(val1$RAmin, val1$Rkmin, val1$Rimin,
                val1$RMmin)
            inputmax <- c(val1$RAmax, val1$Rkmax, val1$Rimax,
                val1$RMmax)
        }
        if (modno == 2 | modno == 22) {
            val1 <- data.frame(value[5L], value[6L], value[7L])
            names(val1) <- c("RAsym", "Rk", "Ri")
            inputmin <- c(val1$RAmin, val1$Rkmin, val1$Rimin)
            inputmax <- c(val1$RAmax, val1$Rkmax, val1$Rimax)
        }
        if (modno == 3 | modno == 23) {
            val1 <- data.frame(value[7L], value[8L])
            names(val1) <- c("Ri", "RM")
            inputmin <- c(val1$Rimin, val1$RMmin)
            inputmax <- c(val1$Rimax, val1$RMmax)
        }
        if (modno == 4 | modno == 24) {
            val1 <- data.frame(value[5L], value[8L])
            names(val1) <- c("RAsym", "RM")
            inputmin <- c(val1$RAmin, val1$RMmin)
            inputmax <- c(val1$RAmax, val1$RMmax)
        }
        if (modno == 5 | modno == 25) {
            val1 <- data.frame(value[8L])
            names(val1) <- c("RM")
            inputmin <- c(val1$RMmin)
            inputmax <- c(val1$RMmax)
        }
        if (modno == 6 | modno == 26 | modno == 17) {
            val1 <- data.frame(value[5L], value[7L], value[8L])
            names(val1) <- c("RAsym", "Ri", "RM")
            inputmin <- c(val1$RAmin, val1$Rimin, val1$RMmin)
            inputmax <- c(val1$RAmax, val1$Rimax, val1$RMmax)
        }
        if (modno == 7 | modno == 27) {
            val1 <- data.frame(value[6L], value[7L], value[8L])
            names(val1) <- c("Rk", "Ri", "RM")
            inputmin <- c(val1$Rkmin, val1$Rimin, val1$RMmin)
            inputmax <- c(val1$Rkmin, val1$Rimax, val1$RMmax)
        }
        if (modno == 8 | modno == 28) {
            val1 <- data.frame(value[5L], value[6L], value[8L])
            names(val1) <- c("RAsym", "Rk", "RM")
            inputmin <- c(val1$RAmin, val1$Rkmin, val1$RMmin)
            inputmax <- c(val1$RAmax, val1$Rkmax, val1$RMmax)
        }
        if (modno == 9 | modno == 29) {
            val1 <- data.frame(value[6L], value[8L])
            names(val1) <- c("Rk", "RM")
            inputmin <- c(val1$Rkmin, val1$RMmin)
            inputmax <- c(val1$Rkmax, val1$RMmax)
        }
        if (modno == 10 | modno == 30) {
            val1 <- data.frame(value[7L])
            names(val1) <- c("Ri")
            inputmin <- c(val1$Rimin)
            inputmax <- c(val1$Rimax)
        }
        if (modno == 11 | modno == 31) {
            val1 <- data.frame(value[5L])
            names(val1) <- c("RAsym")
            inputmin <- c(val1$RAmin)
            inputmax <- c(val1$RAmax)
        }
        if (modno == 12 | modno == 32 | modno == 19) {
        }
        if (modno == 13 | modno == 33) {
            val1 <- data.frame(value[5L], value[7L])
            names(val1) <- c("RAsym", "Ri")
            inputmin <- c(val1$RAmin, val1$Rimin)
            inputmax <- c(val1$RAmax, val1$Rimax)
        }
        if (modno == 14 | modno == 34) {
            val1 <- data.frame(value[6L], value[7L])
            names(val1) <- c("Rk", "Ri")
            inputmin <- c(val1$Rkmin, val1$Rimin)
            inputmax <- c(val1$Rkmax, val1$Rimax)
        }
        if (modno == 15 | modno == 35) {
            val1 <- data.frame(value[5L], value[6L])
            names(val1) <- c("RAsym", "Rk")
            inputmin <- c(val1$RAmin, val1$Rkmin)
            inputmax <- c(val1$RAmax, val1$Rkmax)
        }
        if (modno == 16 | modno == 36) {
            val1 <- data.frame(value[6L])
            names(val1) <- c("Rk")
            inputmin <- c(val1$Rkmin)
            inputmax <- c(val1$Rkmax)
        }
        if (modno == 20) {
        }
    } else {
        xyL$y <- (xyL$y - (xyL$y[1] + 1)) * -1
        rng <- range(xyL$y)
        drng <- diff(rng)
        xyL$prop <- (xyL$y - rng[1L] + 0.05 * drng)/(1.1 * drng)
        for (xk in 1:3) {
            if (xk == 1)
                {
                  transy <- (log(xyL$prop/(1 - xyL$prop)))
                  irL <- lm(x ~ transy, data = xyL)
                  ssqL <- sum((irL$residuals)^2)
                }
            if (xk == 2) {
                transy <- (-log(8/(27 * xyL$prop)))
                irV <- lm(x ~ transy, data = xyL)
                ssqV <- sum((irV$residuals)^2)
            }
            if (xk == 3) {
                transy <- (-log(-log(xyL$prop)))
                irG <- lm(x ~ transy, data = xyL)
                ssqG <- sum((irG$residuals)^2)
            }
        }
        if (ssqL < ssqV & ssqL < ssqG) {
            ir <- as.vector(irL$coefficients)
            RM <- 1
        } else {
            if (ssqV < ssqL & ssqV < ssqG) {
                ir <- as.vector(irV$coefficients)
                RM = -0.3
            } else {
                ir <- as.vector(irG$coefficients)
                RM = 0.1
            }
        }
        Ri <- ir[1L]
        Rk <- 1/(ir[2L])
        RAsym <- rng[2L]
        if (modno == 18 | modno == 19) {
            RAmax = RAsym + (abs(RAsym) * 2.5)
            RAmin = RAsym - (abs(RAsym) * 0.5)
            Rkmax = Rk + (abs(Rk) * 0.5)
            Rkmin = Rk - (abs(Rk) * 0.5)
            Rimax = Ri + (abs(Ri) * 10)
            Rimin = Ri + (abs(Ri) * -2.5)
            while (abs(Rimax * Rkmax) > 700) Rimax = Rimax *
                0.9
            while (abs(Rimin * Rkmax) > 700) Rimin = Rimin *
                0.9
            RMmax = RM + abs(RM * 2)
            RMmin = RM - abs(M * 2)
            skel <- rep(list(1), 11)
            exportparams <- c(Asym, K, Infl, M, (Asym * 0.05),
                K, tAsym, M, NA, NA, NA)
            exportparams <- relist(exportparams, skel)
            names(exportparams) <- c("Asym", "K", "Infl", "M",
                "RAsym", "Rk", "Ri", "RM", "first_y", "last_y",
                "x_at_last_y")
            skel <- rep(list(1), 16)
            exportparamsbounds <- c(Amin, Amax, Kmin, Kmax, Imin,
                Imax, Mmin, Mmax, RAmin, RAmax, Rkmin, Rkmax,
                Rimin, Rimax, RMmin, RMmax)
            exportparamsbounds <- relist(exportparamsbounds,
                skel)
            names(exportparamsbounds) <- c("Amin", "Amax", "Kmin",
                "Kmax", "Imin", "Imax", "Mmin", "Mmax", "RAmin",
                "RAmax", "Rkmin", "Rkmax", "Rimin", "Rimax",
                "RMmin", "RMmax")
            assign("pnmodelparams", exportparams, envir = globalenv())
            assign("pnmodelparamsbounds", exportparamsbounds,
                envir = globalenv())
        }
        if (modno == 1 | modno == 21 | modno == 18) {
            inputval <- c(RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM)
            inputmin <- c(RAmin, Rkmin, Rimin, RMmin)
            inputmax <- c(RAmax, Rkmax, Rimax, RMmax)
        }
        if (modno == 2 | modno == 22) {
            RM <- modelparams$RM
            inputval <- c(RAsym = RAsym, Rk = Rk, Ri = Ri)
            inputmin <- c(RAmin, Rkmin, Rimin)
            inputmax <- c(RAmax, Rkmax, Rimax)
        }
        if (modno == 3 | modno == 23) {
            Rk <- modelparams$Rk
            RAsym <- modelparams$RAsym
            inputval <- c(Ri = Ri, RM = RM)
            inputmin <- c(Rimin, RMmin)
            inputmax <- c(Rimax, RMmax)
        }
        if (modno == 4 | modno == 24) {
            Rk <- modelparams$Rk
            Ri <- modelparams$Ri
            inputval <- c(RAsym = RAsym, RM = RM)
            inputmin <- c(RAmin, RMmin)
            inputmax <- c(RAmax, RMmax)
        }
        if (modno == 5 | modno == 25) {
            RAsym <- modelparams$RAsym
            Rk <- modelparams$Rk
            Ri <- modelparams$Ri
            inputval <- c(RM = RM)
            inputmin <- c(RMmin)
            inputmax <- c(RMmax)
        }
        if (modno == 6 | modno == 26 | modno == 17) {
            Rk <- modelparams$Rk
            inputval <- c(RAsym = RAsym, Ri = Ri, RM = RM)
            inputmin <- c(RAmin, Rimin, RMmin)
            inputmax <- c(RAmax, Rimax, RMmax)
        }
        if (modno == 7 | modno == 27) {
            RAsym <- modelparams$RAsym
            inputval <- c(Rk = Rk, Ri = Ri, RM = RM)
            inputmin <- c(Rkmin, Rimin, RMmin)
            inputmax <- c(Rkmax, Rimax, RMmax)
        }
        if (modno == 8 | modno == 28) {
            Ri <- modelparams$Ri
            inputval <- c(RAsym = RAsym, Rk = Rk, RM = RM)
            inputmin <- c(RAmin, Rkmin, RMmin)
            inputmax <- c(RAmax, Rkmax, RMmax)
        }
        if (modno == 9 | modno == 29) {
            RAsym <- modelparams$RAsym
            Ri <- modelparams$Ri
            inputval <- c(Rk = Rk, RM = RM)
            inputmin <- c(Rkmin, RMmin)
            inputmax <- c(Rkmax, RMmax)
        }
        if (modno == 10 | modno == 30) {
            Rk <- modelparams$Rk
            RAsym <- modelparams$RAsym
            RM <- modelparams$RM
            inputval <- c(Ri = Ri)
            inputmin <- c(Rimin)
            inputmax <- c(Rimax)
        }
        if (modno == 11 | modno == 31) {
            Rk <- modelparams$Rk
            Ri <- modelparams$Ri
            RM <- modelparams$RM
            inputval <- c(RAsym = RAsym)
            inputmin <- c(RAmin)
            inputmax <- c(RAmax)
        }
        if (modno == 12 | modno == 32 | modno == 19) {
            print("error in model 12/32 code - please report")
        }
        if (modno == 13 | modno == 33) {
            Rk <- modelparams$Rk
            RM <- modelparams$RM
            inputval <- c(RAsym = RAsym, Ri = Ri)
            inputmin <- c(RAmin, Rimin)
            inputmax <- c(RAmax, Rimax)
        }
        if (modno == 14 | modno == 34) {
            RAsym <- modelparams$RAsym
            RM <- modelparams$RM
            inputval <- c(Rk = Rk, Ri = Ri)
            inputmin <- c(Rkmin, Rimin)
            inputmax <- c(Rkmax, Rimax)
        }
        if (modno == 15 | modno == 35) {
            Ri <- modelparams$Ri
            RM <- modelparams$RM
            inputval <- c(RAsym = RAsym, Rk = Rk)
            inputmin <- c(RAmin, Rkmin)
            inputmax <- c(RAmax, Rkmax)
        }
        if (modno == 16 | modno == 36) {
            RAsym <- modelparams$RAsym
            Ri <- modelparams$Ri
            RM <- modelparams$RM
            inputval <- c(Rk = Rk)
            inputmin <- c(Rkmin)
            inputmax <- c(Rkmax)
        }
        if (modno == 20) {
            print("error in model 20 code - please report")
        }
        func2 <- function(val1) {
            pnmodelparams <- get("pnmodelparams", envir = .GlobalEnv)
            value <- data.frame(RAsym = pnmodelparams$RAsym,
                Rk = pnmodelparams$Rk, Ri = pnmodelparams$Ri,
                RM = pnmodelparams$RM)
            val1 <- (data.frame(t(val1)))
            if (is.null(val1$RAsym) == FALSE)
                value$RAsym <- val1$RAsym
            if (is.null(val1$Rk) == FALSE)
                value$Rk <- val1$Rk
            if (is.null(val1$Ri) == FALSE)
                value$Ri <- val1$Ri
            if (is.null(val1$RM) == FALSE)
                value$RM <- val1$RM
            RAsym <- value$RAsym
            Rk <- value$Rk
            Ri <- value$Ri
            RM <- value$RM
            if (modno == 17) {
                func4 <- function(x, Ri, RM) {
                  x = x + (0+0i)
                  (1 + exp(Ri - x)/M)
                }
                P <- Re(func4(xyL$x, Ri, RM))
            } else {
                func41 <- function(x, Rk, Ri, RM) {
                  x = x + (0+0i)
                  (1 + RM * exp(-Rk * (x - Ri)))^(1/RM)
                }
                P <- Re(func41(xyL$x, Rk, Ri, RM))
            }
            P[is.na(P)] <- 1e-290 * pnmodelparams$RAsym
            P[P == Inf] <- 1e-290 * pnmodelparams$RAsym
            y1 <- RAsym/P
            evl <- sum((xyL$y - y1)^2)
            if (evl == Inf)
                evl <- 1e+290
            options(warn = -1)
            try(if (min(Im(as.complex(1 + RM * exp(-Rk * ((0:200) -
                Ri)))^(1/RM))) < 0) {
                evl <- 1e+290
            }, silent = TRUE)
            options(warn = 0)
            return(evl)
        }
        oppar = 0
        is.na(oppar) <- TRUE
        oppar <- optim(inputval, func2, method = "L-BFGS-B",
            lower = inputmin, upper = inputmax, control = list(maxit = 2000))
        if (is.null(oppar[1]) == FALSE) {
            if (oppar$convergence < 52) {
                val1 <- data.frame(t(c(oppar$par)))
            } else {
                print("Warning: negative optimization failed using default parameters (straight line)")
                val1 <- data.frame((Asym * 0.05), value[2L],
                  tAsym, value[4L])
                inputval <- data.frame(t(inputval))
                if (is.null(inputval$RAsym) == TRUE)
                  val1[1] = (-999)
                if (is.null(inputval$Rk) == TRUE)
                  val1[2] = (-999)
                if (is.null(inputval$Ri) == TRUE)
                  val1[3] = (-999)
                if (is.null(inputval$RM) == TRUE)
                  val1[4] = (-999)
                names(val1) <- c("RAsym", "Rk", "Ri", "RM")
                savname <- names(val1[, val1 > -999])
                val1 <- data.frame((val1))
                val1 <- val1[, val1 > -999]
                names(val1) <- savname
            }
        } else {
            print("Warning: negative optimization failed using default parameters (straight line)")
            val1 <- c((Asym * 0.05), value[2L], tAsym, value[4L])
            inputval <- data.frame(t(inputval))
            if (is.null(inputval$RAsym) == TRUE)
                val1[1] = (-999)
            if (is.null(inputval$Rk) == TRUE)
                val1[2] = (-999)
            if (is.null(inputval$Ri) == TRUE)
                val1[3] = (-999)
            if (is.null(inputval$RM) == TRUE)
                val1[4] = (-999)
            names(val1) <- c("RAsym", "Rk", "Ri", "RM")
            savname <- names(val1[, val1 > -999])
            val1 <- val1[, val1 > -999]
            val1 <- data.frame((val1))
            names(val1) <- savname
        }
    }
    if (modno == 12 | modno == 20 | modno == 32 | modno == 19) {
        if (modno == 12 | modno == 19) {
            value <- value[1:4]
        } else {
            value <- value[1:3]
        }
    } else {
        if (is.null(val1$RAsym) == FALSE) {
            RAsym <- val1$RAsym
        } else {
            RAsym <- modelparams$RAsym
        }
        if (is.null(val1$Rk) == FALSE)
            Rk <- {
                val1$Rk
            } else {
            Rk <- modelparams$Rk
        }
        if (is.null(val1$Ri) == FALSE)
            Ri <- {
                val1$Ri
            } else {
            Ri <- modelparams$Ri
        }
        if (is.null(val1$RM) == FALSE)
            RM <- {
                val1$RM
            } else {
            RM <- modelparams$RM
        }
        value <- data.frame(value[1L], value[2L], value[3L],
            value[4L], RAsym, Rk, Ri, RM)
        Asym = value[1]
        K = value[2]
        Infl = value[3]
        M = value[4]
        RAsym = value[5]
        Rk = value[6]
        Ri = value[7]
        RM = value[8]
        names(value) <- c("Asym", "K", "Infl", "M", "RAsym",
            "Rk", "Ri", "RM")
        if (is.null(val1$RAsym) == FALSE) {
            RAsym <- val1$RAsym
        } else {
            value[5] <- (-999)
        }
        if (is.null(val1$Rk) == FALSE)
            Rk <- {
                val1$Rk
            } else {
            value[6] <- (-999)
        }
        if (is.null(val1$Ri) == FALSE)
            Ri <- {
                val1$Ri
            } else {
            value[7] <- (-999)
        }
        if (is.null(val1$RM) == FALSE)
            RM <- {
                val1$RM
            } else {
            value[8] <- (-999)
        }
        if (modno > 19)
            value[4] <- (-999)
        if (modno == 17)
            value[2] <- (-999)
        savname <- names(value[, value > -999])
        value <- value[, value > -999]
        names(value) <- savname
        if (modno > 19) {
            posmin <- c(Amin, Kmin, Imin)
            posmax <- c(Amax, Kmax, Imax)
        } else {
            if (modno == 17) {
                posmin <- c(Amin, Imin, Mmin)
                posmax <- c(Amax, Imax, Mmax)
            } else {
                posmin <- c(Amin, Kmin, Imin, Mmin)
                posmax <- c(Amax, Kmax, Imax, Mmax)
            }
        }
        finalpars <- 0
        is.na(finalpars) <- TRUE
        richardsR <- function(Rparams) {
            val2 <- data.frame(Asym = modelparams$Asym, K = modelparams$K,
                Infl = modelparams$Infl, M = modelparams$M, RAsym = modelparams$RAsym,
                Rk = modelparams$Rk, Ri = modelparams$Ri, RM = modelparams$RM)
            val3 <- (data.frame(t(Rparams)))
            if (is.null(val3$Asym) == FALSE)
                val2$Asym <- val3$Asym
            if (is.null(val3$K) == FALSE)
                val2$K <- val3$K
            if (is.null(val3$Infl) == FALSE)
                val2$Infl <- val3$Infl
            if (is.null(val3$M) == FALSE)
                val2$M <- val3$M
            if (is.null(val3$RAsym) == FALSE)
                val2$RAsym <- val3$RAsym
            if (is.null(val3$Rk) == FALSE)
                val2$Rk <- val3$Rk
            if (is.null(val3$Ri) == FALSE)
                val2$Ri <- val3$Ri
            if (is.null(val3$RM) == FALSE)
                val2$RM <- val3$RM
            Asym <- val2$Asym
            K <- val2$K
            Infl <- val2$Infl
            M <- val2$M
            RAsym <- val2$RAsym
            Rk <- val2$Rk
            Ri <- val2$Ri
            RM <- val2$RM
            if (is.na(exp(-K * (min(xy$x) - Infl))) == TRUE |
                (exp(-K * (min(xy$x) - Infl))) == Inf) {
                K = modelparams$K
                Infl = modelparams$Infl
            }
            if (is.na(exp(-K * (max(xy$x) - Infl))) == TRUE |
                (exp(-K * (min(xy$x) - Infl))) == Inf) {
                K = modelparams$K
                Infl = modelparams$Infl
            }
            if (is.na(exp(-Rk * (min(xy$x) - Ri))) == TRUE |
                (exp(-Rk * (min(xy$x) - Ri))) == Inf) {
                Rk = modelparams$Rk
                Ri = modelparams$Ri
            }
            if (is.na(exp(-Rk * (max(xy$x) - Ri))) == TRUE |
                (exp(-Rk * (min(xy$x) - Ri))) == Inf) {
                Rk = modelparams$Rk
                Ri = modelparams$Ri
            }
            options(warn = -1)
            if (Re(as.complex(1 + M * exp(-K * (xy$x - Infl)))) <
                0) {
                if (Re(as.complex(1 + RM * exp(-Rk * (xy$x -
                  Ri)))) < 0) {
                  if (modno == 17) {
                    y1 <- SSposnegRichardsF17(xy$x, Asym, Infl,
                      M, RAsym, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)),
                      Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  } else {
                    y1 <- SSposnegRichardsFMRM(xy$x, Asym, K,
                      Infl, M, RAsym, Rk, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsFMRM((0:max(xy$x)),
                      Asym, K, Infl, M, RAsym, Rk, Ri, RM)) <
                      0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  }
                } else {
                  if (modno == 17) {
                    y1 <- SSposnegRichardsF17(xy$x, Asym, Infl,
                      M, RAsym, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)),
                      Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  } else {
                    y1 <- SSposnegRichardsFM(xy$x, Asym, K, Infl,
                      M, RAsym, Rk, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsFM((0:max(xy$x)),
                      Asym, K, Infl, M, RAsym, Rk, Ri, RM)) <
                      0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  }
                }
            } else {
                if (Re(as.complex(1 + RM * exp(-Rk * (xy$x -
                  Ri)))) < 0) {
                  if (modno == 17) {
                    y1 <- SSposnegRichardsF17(xy$x, Asym, Infl,
                      M, RAsym, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)),
                      Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  } else {
                    y1 <- SSposnegRichardsFRM(xy$x, Asym, K,
                      Infl, M, RAsym, Rk, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsFRM((0:max(xy$x)),
                      Asym, K, Infl, M, RAsym, Rk, Ri, RM)) <
                      0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  }
                } else {
                  if (modno == 17) {
                    y1 <- SSposnegRichardsF17(xy$x, Asym, Infl,
                      M, RAsym, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsF17((0:max(xy$x)),
                      Asym, Infl, M, RAsym, Ri, RM)) < 0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  } else {
                    y1 <- SSposnegRichardsF(xy$x, Asym, K, Infl,
                      M, RAsym, Rk, Ri, RM)
                    y1[is.na(y1)] <- 1e-290 * pnmodelparams$RAsym
                    y1[y1 == Inf] <- 1e-290 * pnmodelparams$RAsym
                    evl <- sum((xy$y - y1)^2)
                    if (evl == Inf)
                      evl <- 1e+290
                    try(if (min(Im(SSposnegRichardsF((0:max(xy$x)),
                      Asym, K, Infl, M, RAsym, Rk, Ri, RM)) <
                      0)) {
                      evl <- 1e+200
                    }, silent = TRUE)
                  }
                }
            }
            if (abs(evl) == Inf)
                evl <- 1e+290
            if (is.na(eval) == TRUE)
                evl <- 1e+290
            options(warn = 0)
            return(evl)
        }
        upbnds <- data.frame((value[1:length(value)]))
        dnbnds <- upbnds
        if ("Asym" %in% names(value)) {
            upbnds$Asym = Asym + (abs(pnmodelparamsbounds$Amax -
                pnmodelparamsbounds$Amin) * 0.1)
        }
        if ("Asym" %in% names(value)) {
            dnbnds$Asym = Asym - (abs(pnmodelparamsbounds$Amax -
                pnmodelparamsbounds$Amin) * 0.1)
        }
        if ("K" %in% names(value)) {
            upbnds$K = K + (abs(pnmodelparamsbounds$Kmax - pnmodelparamsbounds$Kmin) *
                0.25)
        }
        if ("K" %in% names(value)) {
            dnbnds$K = K - (abs(pnmodelparamsbounds$Kmax - pnmodelparamsbounds$Kmin) *
                0.25)
        }
        if ("Infl" %in% names(value)) {
            upbnds$Infl = Infl + (abs(pnmodelparamsbounds$Imax -
                pnmodelparamsbounds$Imin) * 0.35)
        }
        if ("Infl" %in% names(value)) {
            dnbnds$Infl = Infl - (abs(pnmodelparamsbounds$Imax -
                pnmodelparamsbounds$Imin) * 0.2)
        }
        if ("M" %in% names(value)) {
            upbnds$M = M + abs(pnmodelparamsbounds$Mmax - pnmodelparamsbounds$Mmin) *
                0.2
        }
        if ("M" %in% names(value)) {
            dnbnds$M = M - abs(pnmodelparamsbounds$Mmax - pnmodelparamsbounds$Mmin) *
                0.2
        }
        if ("RAsym" %in% names(value)) {
            upbnds$RAsym = RAsym + (abs(pnmodelparamsbounds$RAmax -
                pnmodelparamsbounds$RAmin) * 0.1)
        }
        if ("RAsym" %in% names(value)) {
            dnbnds$RAsym = RAsym - (abs(pnmodelparamsbounds$RAmax -
                pnmodelparamsbounds$RAmin) * 0.1)
        }
        if ("Rk" %in% names(value)) {
            upbnds$Rk = Rk + (abs(pnmodelparamsbounds$Rkmax -
                pnmodelparamsbounds$Rkmin) * 0.25)
        }
        if ("Rk" %in% names(value)) {
            dnbnds$Rk = Rk - (abs(pnmodelparamsbounds$Rkmax -
                pnmodelparamsbounds$Rkmin) * 0.25)
        }
        if ("Ri" %in% names(value)) {
            upbnds$Ri = Ri + (abs(pnmodelparamsbounds$Rimax -
                pnmodelparamsbounds$Rimin) * 0.35)
        }
        if ("Ri" %in% names(value)) {
            dnbnds$Ri = Ri - (abs(pnmodelparamsbounds$Rimax -
                pnmodelparamsbounds$Rimin) * 0.2)
        }
        if ("RM" %in% names(value)) {
            upbnds$RM = RM + abs(pnmodelparamsbounds$RMmax -
                pnmodelparamsbounds$RMmin) * 0.2
        }
        if ("RM" %in% names(value)) {
            dnbnds$RM = RM - abs(pnmodelparamsbounds$RMmax -
                pnmodelparamsbounds$RMmin) * 0.2
        }
        options(warn = -1)
        cnt <- 0
        finalpars <- 0
        is.na(finalpars) <- TRUE
        oppar1 <- 0
        is.na(oppar1) <- TRUE
        names(oppar1) <- "convergence"
        while (is.na(finalpars[1L]) == TRUE) {
            if (cnt > 1) {
                if ("Rk" %in% names(value)) {
                  if (cnt == 2) {
                    value$Rk <- pnmodelparamsbounds$Rkmin
                    dnbnds$Rk <- pnmodelparamsbounds$Rkmin -
                      abs(pnmodelparams$Rk) * 0.5
                    upbnds$Rk <- pnmodelparamsbounds$Rkmin +
                      abs(pnmodelparams$Rk) * 0.5
                  }
                } else {
                  cnt = cnt + 1
                }
                if ("Rk" %in% names(value)) {
                  if (cnt == 3) {
                    value$Rk <- pnmodelparamsbounds$Rkmax
                    dnbnds$Rk <- pnmodelparamsbounds$Rkmax -
                      abs(pnmodelparams$Rk) * 0.5
                    upbnds$Rk <- pnmodelparamsbounds$Rkmax +
                      abs(pnmodelparams$Rk) * 0.5
                  }
                } else {
                  cnt = cnt + 1
                }
                if (cnt == 4 & modno != 17) {
                  value$Rk <- pnmodelparams$Rk
                  dnbnds$Rk <- pnmodelparamsbounds$Rkmin
                  upbnds$Rk <- pnmodelparamsbounds$Rkmax
                  value$K <- pnmodelparamsbounds$Kmin
                  dnbnds$K <- pnmodelparamsbounds$Kmin - abs(pnmodelparams$K) *
                    0.5
                  upbnds$K <- pnmodelparamsbounds$Kmin + abs(pnmodelparams$K) *
                    0.5
                }
                if (cnt == 5 & modno != 17) {
                  value$K <- pnmodelparamsbounds$Kmax
                  dnbnds$K <- pnmodelparamsbounds$Kmax - abs(pnmodelparams$K) *
                    0.5
                  upbnds$K <- pnmodelparamsbounds$Kmax + abs(pnmodelparams$K) *
                    0.5
                }
            }
            oppar1 <- data.frame(52)
            names(oppar1) <- c("convergence")
            try(oppar1 <- (optim(value, richardsR, method = "L-BFGS-B",
                lower = c(dnbnds), upper = c(upbnds), control = list(maxit = 2000))),
                silent = FALSE)
            cnt <- cnt + 1
            if (oppar1$convergence < 52) {
                finalpars <- data.frame(t(c(oppar1$par)))
            } else {
                if (cnt > 4) {
                  finalpars <- 1
                } else {
                  is.na(finalpars) <- TRUE
                }
            }
        }
        if (cnt > 4)
            is.na(finalpars) <- TRUE
        options(warn = 0)
        if (is.na(finalpars[1] == TRUE)) {
            print("Warning: simultaneous optimization of pre- and post- peak curves failed, using separately fitted parameters")
            if (is.na(inputval[1]) == TRUE) {
                if (modno > 19) {
                  value <- c(Asym, K, Infl, val1)
                } else {
                  if (modno == 17) {
                    value <- c(Asym, Infl, M, val1)
                  } else {
                    value <- c(Asym, K, Infl, M, val1)
                  }
                }
            } else {
                if (modno > 19) {
                  value <- c(Asym, K, Infl, inputval)
                } else {
                  if (modno == 17) {
                    value <- c(Asym, Infl, M, inputval)
                  } else {
                    value <- c(Asym, K, Infl, M, inputval)
                  }
                }
            }
        } else {
            finpars <- finalpars
            if ("Rk" %in% names(finpars))
                {
                  if (finpars$Rk < 1e-04 & finpars$Rk > -1e-04) {
                    finpars$Rk <- finalpars$Rk
                    finpars[1L] <- finalpars[1L] + (finalpars$Rk/(finalpars[1L]/finalpars$Rk))
                  }
                }
            value <- c(finpars)
        }
    }
    if (modno > 19) {
        vnames <- c("Asym", "K", "Infl")
    } else {
        vnames <- c("Asym", "K", "Infl", "M")
    }
    if (modno == 1 | modno == 21 | modno == 18) {
        names(value) <- mCall[c(vnames, "RAsym", "Rk", "Ri",
            "RM")]
    }
    if (modno == 2 | modno == 22) {
        names(value) <- mCall[c(vnames, "RAsym", "Rk", "Ri")]
    }
    if (modno == 3 | modno == 23) {
        names(value) <- mCall[c(vnames, "Ri", "RM")]
    }
    if (modno == 4 | modno == 24) {
        names(value) <- mCall[c(vnames, "RAsym", "RM")]
    }
    if (modno == 5 | modno == 25) {
        names(value) <- mCall[c(vnames, "RM")]
    }
    if (modno == 6 | modno == 26) {
        names(value) <- mCall[c(vnames, "RAsym", "Ri", "RM")]
    }
    if (modno == 7 | modno == 27) {
        names(value) <- mCall[c(vnames, "Rk", "Ri", "RM")]
    }
    if (modno == 8 | modno == 28) {
        names(value) <- mCall[c(vnames, "RAsym", "Rk", "RM")]
    }
    if (modno == 9 | modno == 29) {
        names(value) <- mCall[c(vnames, "Rk", "RM")]
    }
    if (modno == 10 | modno == 30) {
        names(value) <- mCall[c(vnames, "Ri")]
    }
    if (modno == 11 | modno == 31) {
        names(value) <- mCall[c(vnames, "RAsym")]
    }
    if (modno == 12 | modno == 32 | modno == 19) {
        names(value) <- mCall[vnames]
    }
    if (modno == 13 | modno == 33) {
        names(value) <- mCall[c(vnames, "RAsym", "Ri")]
    }
    if (modno == 14 | modno == 34) {
        names(value) <- mCall[c(vnames, "Rk", "Ri")]
    }
    if (modno == 15 | modno == 35) {
        names(value) <- mCall[c(vnames, "RAsym", "Rk")]
    }
    if (modno == 16 | modno == 36) {
        names(value) <- mCall[c(vnames, "Rk")]
    }
    if (modno == 20) {
        names(value) <- mCall[c("Asym", "K", "Infl")]
    }
    if (modno == 17) {
        names(value) <- mCall[c("Asym", "Infl", "M", "RAsym",
            "Ri", "RM")]
    }
    value
}

