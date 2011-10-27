extraF <-
structure(function # Compare Two \eqn{nlsList} Models Using Extra Sum-of-Squares F-Tests
                   (submodel = 1,
                    ### \eqn{nlsList} model with fewer curve parameters (reduced model)
                    genmodel = 1
                    ### \eqn{nlsList} model with more curve parameters (general model)
                    ) {
    ##description<< Function to compare two nested models using extra sum-of-squares F-Tests.                
    ##details<< Models must be entered in the correct order with the reduced model appearing
    ## first in the call and the more general model appearing later. These must be nested models,
    ## i.e. the general model must contain all of the curve parameters in the reduced model and more.
    ## Entering models with the same number of parameters will produce NAs in the output, but
    ## the function will produce seemingly adequate output with non-nested models. The user must
    ## check that models are nested prior to use.
    ##
    ## This function is primarily designed to be called by the model selection functions
    ## \code{\link{pn_modselect_step}} and \code{\link{pn_mod_compare}} but can be used independently.
    ##
    ## Extra sum-of-squares is obtained from: \preformatted{F = (SS1 - SS2)(df1 - df2) / (SS2 / df2)}
    ## where SS = sum-of-squares and df = degrees of freedom, for the more reduced model (1) and the
    ## more general model (2), respectively.
    ##
    ## If the F value is significant then the more general model provides a significant improvement
    ## over the reduced model, but if the models are not significantly different then the reduced
    ## parameter model is to be preferred.
    assign("legitmodel", "legitmodelreset", envir = .GlobalEnv)
    chk <- try(unlist(summary(submodel))["RSE"], silent = TRUE)
    if (class(chk)[1] == "try-error" | class(submodel)[1] != "nlsList" )
        submodel <- 1
    chk1 <- try(unlist(summary(genmodel))["RSE"], silent = TRUE)
    if (class(chk1)[1] == "try-error" | class(genmodel)[1] != "nlsList" )
        genmodel <- 1
    if (class(submodel)[1] == "numeric" | class(genmodel)[1] ==
        "numeric") {
        if (class(submodel)[1] == "numeric" & class(genmodel)[1] ==
            "numeric") {
            assign("legitmodel", 1, envir = .GlobalEnv)
        output <- data.frame(NA, NA, NA, NA, NA, NA)
        } else {        
         if (class(submodel)[1] == "numeric") {
            is.na(submodel) <- TRUE
            assign("legitmodel", genmodel, envir = .GlobalEnv)
	 df1 <- as.numeric(unlist(summary(genmodel)["df.residual"]))
	 sdf<- unlist(coef(genmodel)[1])
	 lsdf <- length(sdf[!is.na(sdf)])
	 if(lsdf == 0) lsdf =1
	 df1 <- (df1 / lsdf ) * length(sdf)
	 df1 <- as.integer(round(df1))
         RSSa <- (as.numeric(unlist(summary(genmodel))["RSE"])^2) *
            df1         
         output <- data.frame(NA, df1, NA, NA, RSSa, NA)
         }
         if (class(genmodel)[1] == "numeric") {
            is.na(genmodel) <- TRUE
            assign("legitmodel", submodel, envir = .GlobalEnv)
	 df2 <- as.numeric(unlist(summary(submodel)["df.residual"]))
	 sdf<- unlist(coef(submodel)[1])
	 lsdf <- length(sdf[!is.na(sdf)])
	 if(lsdf == 0) lsdf =1
	 df2 <- (df2 / lsdf ) * length(sdf)
	 df2 <- as.integer(round(df2))			
         RSSb <- (as.numeric(unlist(summary(submodel))["RSE"])^2) *
            df2
         output <- data.frame(NA, df2, NA, NA, NA, RSSb)    
         }
        }
     } else {
        df1 <- as.numeric(unlist(summary(genmodel)["df.residual"]))
	sdf<- unlist(coef(genmodel)[1])
	lsdf <- length(sdf[!is.na(sdf)])
	if(lsdf == 0) lsdf =1
	df1 <- (df1 / lsdf ) * length(sdf)
	df1 <- as.integer(round(df1))
        RSSa <- (as.numeric(unlist(summary(genmodel))["RSE"])^2) *
            df1
        df2 <- as.numeric(unlist(summary(submodel)["df.residual"]))
	sdf<- unlist(coef(submodel)[1])
	lsdf <- length(sdf[!is.na(sdf)])
	if(lsdf == 0) lsdf =1
	df2 <- (df2 / lsdf ) * length(sdf)
	df2 <- as.integer(round(df2))
        RSSb <- (as.numeric(unlist(summary(submodel))["RSE"])^2) *
            df2
        F <- (RSSb - RSSa)/(df2 - df1)/(RSSa/df1)
        sigF <- 1 - pf(abs(F), df1, df2)
        output <- data.frame(F, df2 - df1, df2, sigF, RSSa, RSSb)
    }
    names(output) <- c("Fstat", "df_n", "df_d", "P", "RSS_gen",
        "RSS_sub")
    row.names(output) <- paste(paste(paste("test of", as.character(substitute(submodel)),
        sep = " "), "vs.", sep = " "), as.character(substitute(genmodel)),
        sep = " ")
    return(output)
##value<< A \code{\link{data.frame}} listing the names of the models compared, F,
## numerator degrees of freedom,
## demonimator degrees of freedom, P value and the residual sum of squares for both the general
## and reduced models
##references<< Ritz, C. and Streibigg, J. C. (2008) \eqn{Nonlinear regression with R.}
## Springer-Verlag, New York.
##seealso<< \code{\link{extraF.nls}}
##  \code{\link{nlsList}}
##  \code{\link{pn_modselect_step}}
##  \code{\link{pn_mod_compare}}
}
, ex = function(){
   #compare two nested nlsList models (4 vs 8 parameter models)
   data(posneg_data)
   modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
   # (only first 4 group levels in data used for example's sake)
   subdata<-subset(posneg_data, as.numeric(row.names (posneg_data) ) < 53)
   richardsR1.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1)
                        , data = subdata)
   richardsR12.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = 1, Rk = 1, Ri = 1, RM = 1, modno = 12)
                        , data = subdata)
   extraF(richardsR12.lis, richardsR1.lis)
}
)
