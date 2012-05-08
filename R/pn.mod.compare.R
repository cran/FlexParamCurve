pn.mod.compare <-
structure(function # Compare All Possible Positive-Negative Richards \eqn{nlslist} Models
                         (x,
                          ### a numeric vector of the primary predictor
                          y,
                          ### a numeric vector of the response variable
                          grp,
                          ### a factor of same length as x and y that distinguishes groups within
                          ### the dataset
                          forcemod = 0,
                          ### optional numeric value to constrain model selection (see Details)
                          existing = FALSE,
                          ### optional logical value specifying whether some of the relevant models
                          ### have already been fitted
                          penaliz = "1/sqrt(n)",
                          ### optional character value to determine how models are ranked (see Details)
                          pn.options
                          ### required character value for name of list object populated with starting 
                          ### parameter estimates, fitting 
                          ### options and bounds or destination for modpar to write a new list (see Details)
                          ) {
##description<< This function performs model selection for \code{\link{nlsList}} models fitted using
## \code{\link{SSposnegRichards}}.
##details<< First, whether parameter M should be fixed
## (see \code{\link{SSposnegRichards}}) is determined by fitting models 12 and 20 and comparing
## their perfomance using \code{\link{extraF}}.
## If model 12 provides superior performance (variable values of  M) then 16 models that estimate M
## are run
## (models 1 through 16), otherwise the models with fixed M are fitted (models 21 through 36).
## Fitting these \code{\link{nlsList}} models can be time-consuming (2-4 hours using the dataset
## \code{\link{posneg.data}} that encompasses 100 individuals) and if several of the relevant
## models are already fitted the option existing=TRUE can be used to avoid refitting models that
## already exist globally (note that a model object in which no grouping levels were successfully
## parameterized will be refitted, as will objects that are not of class nlsList).
##
## Specifying forcemod=3 will force model selection to only consider fixed M models and setting
## forcemod=4 will force model selection to consider models with varying values of M only.
## If fitting both models
## 12 and 20 fails, fixed M models will be used by default.
##
## Models are ranked by modified pooled residual square error. By default residual standard error
## is divided by the square root of sample size. This exponentially penalizes models for which very few
## grouping levels (individuals) are successfully parameterized (the few individuals that are
## parameterized in these models are fit unsuprisingly well) using a function based on the relationship
## between standard error and sample size. However, different users may have different preferences
## and these can be specified in the argument penaliz (which residual
## standard error is multiplied by). This argument must be a character value
## that contains the character n (sample size) and must be a valid right hand side (RHS) of a formula:
## e.g. 1*(n), (n)^2. It cannot contain more than one n but could be a custom function, e.g. FUN(n).
    library(nlme)
    pnoptnm <- as.character(pn.options)
    testit <- try(library(nlme), silent = TRUE)
    if (class(testit)[1] == "try-error")
        stop("please load library: nlme")
    checkpen <- try(unlist(strsplit(penaliz, "(n)")), silent = TRUE)
    if (length(checkpen) != 2 | class(checkpen)[1] == "try-error") {
        stop("penaliz parameter is ill defined: see ?pn.mod.compare")
    } else {
        checkpen <- try(eval(parse(text = sprintf("%s", paste(checkpen[1],
            "1", checkpen[2], sep = "")))))
        if (class(checkpen)[1] == "try-error")
            stop("penaliz parameter is ill defined: see ?pn.mod.compare")
    }
    datamerg <- data.frame(x, y, grp)
    userdata <- groupedData(y ~ x | grp, outer = ~grp, data = datamerg)
    testbounds <- 1
    testpar <- 1
    is.na(testbounds) <- TRUE
    is.na(testpar) <- TRUE
    testbounds <- try(get(pnoptnm, envir = .GlobalEnv)[16:32],
        silent = TRUE)
    testpar <- try(get(pnoptnm, envir = .GlobalEnv)[1:15],
        silent = TRUE)
    if (class(testbounds)[1] == "try-error" | class(testpar)[1] ==
        "try-error" | is.na(testbounds[1]) == TRUE | is.na(testpar[1]) ==
        TRUE)
        try(assign(pnoptnm, modpar(datamerg[,1], datamerg[,2], pn.options = ".modcmptemp", verbose=FALSE), .GlobalEnv)
        	, silent = FALSE)
    extraF <- try(get("extraF", envir = .GlobalEnv), silent = TRUE)
    if (class(extraF)[1] == "try-error") {
        stop("cannot find function: extraF - please reload FlexParamCurve")
    }
    assign(".modcmptemp",get(pnoptnm, envir = .GlobalEnv), .GlobalEnv)
    print("checking fit of positive section of the curve for variable M*************************************")
    richardsR12.lis <- try(get("richardsR12.lis", envir = .GlobalEnv),
        silent = TRUE)
    if (class(richardsR12.lis)[1] == "try-error" | existing ==
        FALSE)
        richardsR12.lis <- try(nlsList(y ~ SSposnegRichards(x,
            Asym = Asym, K = K, Infl = Infl, M = M, modno = 12, pn.options = .modcmptemp), data = userdata),
            silent = TRUE)
    print("checking fit of positive section of the curve for fixed M*************************************")
    pnmodelparams <- get(pnoptnm, envir = .GlobalEnv)[1:15]
    change.pnparameters <- try(get("change.pnparameters", envir = .GlobalEnv),
        silent = TRUE)
    chk <- try(unlist(summary(richardsR12.lis))["RSE"], silent = TRUE)
    	if (class(chk)[1] == "try-error" | class(richardsR12.lis)[1] != "nlsList" 
     			 |  is.null(nrow(coef(richardsR12.lis))) == TRUE ){
    	} else {    
    		dummy <- try(change.pnparameters(M = (fixef(richardsR12.lis)[4])),
    		silent = TRUE)
    	} 
    richardsR20.lis <- try(get("richardsR20.lis", envir = .GlobalEnv),
        silent = TRUE)
    if (class(richardsR20.lis)[1] == "try-error" | existing ==
        FALSE)
        richardsR20.lis <- try(nlsList(y ~ SSposnegRichards(x,
            Asym = Asym, K = K, Infl = Infl, modno = 20, pn.options = .modcmptemp), data = userdata),
            silent = TRUE)
    if ((class(richardsR20.lis)[1]) == "try-error") {
        print("3 parameter positive richards model failed*************************************")
        forcemod = 4
        richardsR20.lis <- 1
    }
    if ((class(richardsR12.lis)[1]) == "try-error")
        {
            print("4 parameter positive richards model failed*************************************")
            forcemod = 3
            richardsR12.lis <- 1
        }
    currentmodel <- 1
    testmod <- try(extraF(richardsR20.lis, richardsR12.lis), silent = TRUE)
    if (forcemod == 0) {
        if (class(testmod) == "try-error") {
            modelsig = 0.1
        } else {
            modelsig = testmod[4]
            if ((testmod[4]) > 0.05 & (testmod[5]) > (testmod[6])) {
                currentmodel <- richardsR20.lis
                mostreducednm <- substr("richardsR20.lis", 10,
                  11)
            } else {
                currentmodel <- richardsR12.lis
                mostreducednm <- substr("richardsR12.lis", 10,
                  11)
            }
        }
    }
    mostreducedmod <- currentmodel
    if (class(testmod) != "try-error") {
    mostreducednm <- substr("richardsR20.lis", 10,
                  11)
    mostreducedmod <- richardsR20.lis
    } else {
    mostreducednm <- "NONE"
    }
    if (forcemod == 3)
        {
            modelsig = 0.05
        }
    if (forcemod == 4)
        {
            modelsig = 0.1
        }
    if (modelsig > 0.05) {
        print("Variable M models most appropriate*************************************")
        modno <- c(1:16)
    } else {
        print("Fixed M models most appropriate*************************************")
        modno <- c(21:36)
    }
    runmod <- function(userdata, modelno, modsig, existing = FALSE) {
		savnm <- paste("richardsR", as.character(modelno), ".lis",
		    sep = "")
		if (modelno[1] <20) {
		    savM <- ",M = M"
        } else {
            savM <- ","
        }
        if (existing == FALSE) {
            if (modelno == 1 | modelno == 21)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,Rk=Rk,Ri=Ri,RM=RM,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 2 | modelno == 22)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,Rk=Rk,Ri=Ri,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 3 | modelno == 23)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",Ri=Ri,RM=RM,modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
            if (modelno == 4 | modelno == 24)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,RM=RM,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 5 | modelno == 25)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RM=RM,modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
            if (modelno == 6 | modelno == 26)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,Ri=Ri,RM=RM,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 7 | modelno == 27)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",Rk=Rk,Ri=Ri,RM=RM,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 8 | modelno == 28)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,Rk=Rk,RM=RM,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 9 | modelno == 29)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",Rk=Rk,RM=RM,modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
            if (modelno == 10 | modelno == 30)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",Ri=Ri,modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
            if (modelno == 11 | modelno == 31)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 12 | modelno == 32)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
            if (modelno == 13 | modelno == 33)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,Ri=Ri,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 14 | modelno == 34)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",Rk=Rk,Ri=Ri,modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
            if (modelno == 15 | modelno == 35)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",RAsym=RAsym,Rk=Rk,modno=",
                  modelno, ", pn.options = .modcmptemp),data=userdata)", sep = "")))),
                  silent = TRUE)
            if (modelno == 16 | modelno == 36)
                mod <- try(eval(parse(text = sprintf("%s", paste("nlsList(y~SSposnegRichards(x,Asym=Asym,K=K,Infl=Infl",
                  savM, ",Rk=Rk,modno=", modelno,
                  ", pn.options = .modcmptemp),data=userdata)", sep = "")))), silent = TRUE)
        } else {0
            mod <- try(get(savnm, envir = .GlobalEnv), silent = TRUE)
        }
        if (class(mod)[[1]] != "nlsList")
            mod <- NULL
        checkmod <- try(if (is.null(nrow(coef(mod))) == TRUE) {
            mod <- NULL
        }, silent = TRUE)
        if (class(checkmod)[1] == "try-error" | class(mod)[1] ==
            "NULL") {
            messagesav <- (paste("**********************  Model ",
                savnm, " has not been successfully fit, please trouble-shoot this model separately and then repeat function using existing=TRUE  *************************************************",
                sep = ""))
        } else {
            messagesav <- (paste("**********************  Model ",
                savnm, " fit successfully and saved in Global Environment   *************************************************",
                sep = ""))
            assign(savnm, mod, envir = .GlobalEnv)
            return(messagesav)
        }
    }
    skel <- rep(list(1), 16)
    initval <- c(rep(NA, 16))
    initval <- relist(initval, skel)
    for (i in 1:length(modno)) {
    	print("################  ################  ##################  #################  ###############  #########")
    	print(paste("Fitting model ",i," of ",length(modno),": richardsR",modno[i],".lis",sep=""))
        initval[i] <- runmod(userdata, modno[i], modelsig, existing = existing)
        print(initval[i])
    }
    print("#########################################################################################################")
    print("Model fitting completed, summary of fits follows:")
    print(initval)
    print("#########################################################################################################")
    print("Writing summary tables.....")
    modrank <- data.frame(`PN Richards Model` = rep(NA, length(modno)),
        `Ranking function value` = rep(-999, length(modno)),
        `No. Individuals Fitted` = rep(-999, length(modno)),
        RSE = rep(-999, length(modno)), `df model` = rep(-999,
            length(modno)), `df residual` = rep(-999, length(modno)),
        `No. curve params` = rep(-999, length(modno)))
    for (i in 1:length(modno)) {
        modnm <- paste("richardsR", as.character(modno[i]), ".lis",
            sep = "")
        RSEstr <- "RSE"
        dfstr <- "df"
        dfstr1 <- "df.residual"
        usefun <- unlist(strsplit(penaliz, "(n)"))
        model1 <- try(get(modnm, envir = .GlobalEnv), silent = TRUE)
        checkmod <- try(if (is.null(nrow(coef(model1))) == TRUE) {
            model1 <- NULL
        } else {
            model1 <- model1
        }, silent = TRUE)
        evfun <- parse(text = sprintf("%s", paste("summary(model1)[['",
            RSEstr, "']]*(", usefun[1], "1+sum( summary(model1)[['",
            dfstr, "']],na.rm=TRUE))", usefun[2], sep = "")))
        if (class(model1)[[1]] == "nlsList" & class(checkmod)[1] !=
            "try-error") {
            modrank[i, 2] <- eval(evfun)
            modrank[i, 3] <- nrow(subset(coef(model1)[1], is.na((coef(model1)[1])) ==
                FALSE))
            modrank[i, 4] <- eval(parse(text = sprintf("%s",
                paste("summary(model1)[['", RSEstr, "']]", sep = ""))))
            modrank[i, 5] <- eval(parse(text = sprintf("%s",
                paste("sum(summary(model1)[['", dfstr, "']][,1],na.rm=TRUE)",
                  sep = ""))))
            modrank[i, 6] <- eval(parse(text = sprintf("%s",
                paste("sum(summary(model1)[['", dfstr1, "']],na.rm=TRUE)",
                  sep = ""))))
            modrank[i, 7] <- length(coef(model1))
        } else {
        }
        modrank[i, 1] <- modnm
    }
    pvalmat <- data.frame(matrix(nrow = length(modno), ncol = length(modno)))
    for (i in 1:length(modno)) {
        for (j in 2:(length(modno) - 1)) {
            if (modrank[i, 2] != -999 & modrank[j, 2] != -999) {
                mod1 <- eval(parse(text = sprintf("%s", modrank[i,
                  1])))
                mod2 <- eval(parse(text = sprintf("%s", modrank[j,
                  1])))
                if (length(coef(mod1)) > length(coef(mod2))) {
                  nmmod1 <- names(coef(mod1))
                  nmmod2 <- names(coef(mod2))
                  if (length(nmmod2 %in% nmmod1) == length(nmmod2))
                    pvalmat[i, j] <- as.numeric(extraF(mod2,
                      mod1)[4])
                  if (length(nmmod2 %in% nmmod1) == length(nmmod2))
                    pvalmat[j, i] <- pvalmat[i, j]
                } else {
                  if (length(coef(mod1)) < length(coef(mod2))) {
                    nmmod1 <- names(coef(mod1))
                    nmmod2 <- names(coef(mod2))
                    if (length(nmmod1 %in% nmmod2) == length(nmmod1))
                      pvalmat[i, j] <- as.numeric(extraF(mod1,
                        mod2)[4])
                    if (length(nmmod1 %in% nmmod2) == length(nmmod1))
                      pvalmat[j, i] <- pvalmat[i, j]
                  } else {
                  }
                }
            } else {
            }
        }
        if (modrank[i, 2] != -999) {
        modtemp <- eval(parse(text = sprintf("%s", modrank[i,
            1])))
        modnm <- paste("R", as.character(modno[i]), ".lis", sprintf("(%s)",
            length(coef(modtemp))), sep = "")
        row.names(pvalmat)[i] <- paste("R", as.character(modno[i]),
            sprintf("(%s)", length(coef(modtemp))), sep = "")        
        } else {
        modtemp <- NULL
        modnm <- paste("R", as.character(modno[i]), ".lis", sprintf("(%s)",0), sep = "")
        row.names(pvalmat)[i] <- paste("R", as.character(modno[i]), ".lis", 
        	sprintf("(%s)",0), sep = "")
        }
        names(pvalmat)[i] <- modnm
    }
    pvalmat <- apply(pvalmat, 1, function(x) round(x, 3))
    modrank[modrank == -999] <- NA
    modrank[is.na(modrank[,2]),2] <- Inf
    modrank[is.na(modrank[,6]),6] <- Inf
    modrank <- modrank[order(modrank[, 2], modrank[, 6]), ]
    modrank[(modrank[,2])==Inf,2] <- NA
    modrank[(modrank[,6])==Inf,6] <- NA
    row.names(modrank) <- c(1:16)
    print("...done")
    outp <- list(modrank, pvalmat)
    names(outp) <- c("Model rank table", "P values from pairwise extraF comparisons")
    return(outp)
    ##value<<  A list object with two components: $'Model rank table' contains the
    ## statistics from \code{\link{extraF}} ranked by the  modified residual standard error,
    ## and $'P values from pairwise extraF comparison' is a matrix of P values from
    ## \code{\link{extraF}} for legitimate comparisons (nested and successfully fitted models).
    ## The naming convention for models is a concatenation of 'richardsR', the modno and '.lis'
    ## which is shortened in the matrix output, where the number of parameters has been
    ## pasted in parentheses to allow users to easily distinguish the more general model from
    ## the more reduced model
    ## (see \code{\link{extraF}} and \code{\link{SSposnegRichards}}).
    ##seealso<< \code{\link{extraF}}
    ## \code{\link{SSposnegRichards}}
    ## \code{\link{nlsList}}
    ##note<< If object \eqn{pnmodelparams} does not exist, \code{\link{modpar}}
    ## will be called automatically prior to model selection. During selection, text is output
    ## to the screen to inform the user of the progress of model selection
    ## (which model is being fitted, which were fitted successfully)
}
, ex = function(){
#run model selection for posneg.data object (only first 3 group levels for example's sake)
data(posneg.data)
subdata <- subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
modseltable <- pn.mod.compare(subdata$age, subdata$mass,
    subdata$id, existing = FALSE)
    
#fit nlsList model initially and then run model selection
#for posneg.data object when at least one model is already fit
# note forcemod is set to 3 so that models 21-36 are evaluated
subdata <- subset(posneg.data, as.numeric(row.names (posneg.data) ) < 40)
richardsR22.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, RAsym = RAsym, Rk = Rk, Ri = Ri, modno = 22)
                        ,data = posneg.data)
modseltable <- pn.mod.compare(subdata$age, subdata$mass,
    subdata$id, forcemod = 3, existing = TRUE)
 
#run model selection ranked by residual standard error*sample size
modseltable <- pn.mod.compare(subdata$age, subdata$mass,
    subdata$id, penaliz='1*(n)', existing = TRUE)
}
)