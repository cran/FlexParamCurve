posnegRichards_eqn <-
structure(function # Equations of the FlexParamCurve Family
                       ### Function to solve any of the equations in the FlexParamCurve family,
                       ### depending on user-specified parameters and model choice
                       (x,
                       ### a numeric vector of the primary predictor variable
                        Asym,
                        ### a numeric value for the asymptote of the positive (increasing) curve
                        K,
                        ### a numeric value for the rate parameter of the positive (increasing) curve
                        Infl,
                        ### a numeric value for the point of inflection of the positive (increasing) curve
                        M,
                        ### a numeric value for the shape parameter of the positive (increasing) curve
                        RAsym,
                        ### a numeric value for the asymptote of the negative (decreasing) curve
                        Rk,
                        ### a numeric value for the rate parameter of the negative (decreasing) curve
                        Ri,
                        ### a numeric value for the point of inflection of the negative (decreasing) curve
                        RM,
                        ### a numeric value for the shape parameter of the negative (decreasing) curve
                        modno
                        ### a numeric value (currently integer only) between 1 and 36 specifying the identification
### number of the equation to be fitted
                        ) {
    if (Re(as.complex(1 + M * exp(-K * (max(x) - Infl)))) < 0) {
        fractM <- round(1/M)
    } else {
        fractM <- 1/M
    }
    if (Re(as.complex(1 + RM * exp(-Rk * (max(x) - Ri)))) < 0) {
        fractRM <- round(1/RM)
    } else {
        fractRM <- 1/RM
    }
    if (modno > 19) {
        fractM <- 1/pnmodelparams$M
        M <- pnmodelparams$M
    }
    if (modno == 2 | modno == 22 | modno == 10 | modno == 30 |
        modno == 11 | modno == 31 | modno == 12 | modno == 32 |
        modno == 19 | modno == 13 | modno == 33 | modno == 14 |
        modno == 34 | modno == 15 | modno == 35 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        fractRM <- 1/pnmodelparams$RM
        RM <- pnmodelparams$RM
    }
    if (modno == 3 | modno == 23 | modno == 5 | modno == 25 |
        modno == 7 | modno == 27 | modno == 9 | modno == 29 |
        modno == 10 | modno == 30 | modno == 12 | modno == 32 |
        modno == 19 | modno == 14 | modno == 34 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        RAsym <- pnmodelparams$RAsym
    }
    if (modno == 3 | modno == 23 | modno == 4 | modno == 24 |
        modno == 5 | modno == 25 | modno == 6 | modno == 26 |
        modno == 10 | modno == 30 | modno == 11 | modno == 31 |
        modno == 12 | modno == 32 | modno == 19 | modno == 13 |
        modno == 33 | modno == 20 | modno == 32) {
        Rk <- pnmodelparams$Rk
    }
    if (modno == 4 | modno == 24 | modno == 5 | modno == 25 |
        modno == 8 | modno == 28 | modno == 9 | modno == 29 |
        modno == 11 | modno == 31 | modno == 12 | modno == 32 |
        modno == 19 | modno == 15 | modno == 35 | modno == 16 |
        modno == 36 | modno == 20 | modno == 32) {
        Ri <- pnmodelparams$Ri
    }
    if (modno == 17) {
        (Asym/Re(as.complex(1 + exp(Infl - x)/M))) - (RAsym/Re(as.complex(1 +
            exp(Ri - x)/RM)))
    } else {
        (Asym/Re(as.complex(1 + M * exp(-K * ((x) - Infl)))^(fractM))) -
            (RAsym/Re(as.complex(1 + RM * exp(-Rk * ((x) - Ri)))^(fractRM)))
    }
##value<< the solution of the equation specified (by modno), given the user-entered parameters
##details<< This function fits 1 of 36 possible FlexParamCurve equations. All equations are biomodal,
## in that they have a positive (increasing) trajectory followed by a negative (decreasing) trajectory.
## 
    ## These equations have also been described as double-Richards curves, or positive-negative Richards curves.
    ##
    ## The 36 possible equations are all based on the subtraction of one Richards curve from another, producing:
    ## \eqn{y = A / ([1+ m exp(-k (t-i))]1/m) - A' / ([1+ m' exp(-k' (t-i' ))]1/m' )}, where A=Asym, k=K, i=Infl, m=M,
    ## A'=RAsym, k'=Rk, i'=Ri, m'=RM; as described in the Arguments section above.
    ##
    ## All 36 possible equations are simply reformulations of this equation, in each case fixing a parameter or
    ## multiple parameters to the mean parameter across all individuals in the dataset (such as produced by a \code{\link{nls}}
    ## model). All models are detailed in the \code{\link{SSposnegRichards}} help file. Any models that require parameter fixing
    ## (i.e. all except R1)
    ## extract appropriate values from the object \eqn{pnmodelparams} for the fixed parameters. This object is created by running
    ## \code{\link{modpar}} and can be adjusted manually or by using \code{\link{change.pnparameters}} to user required specification.
    ##note<< Any models that require parameter fixing (i.e. all except #1)
    ## extract appropriate values from the object \eqn{pnmodelparams} for the fixed parameters. This object is created by running
    ## \code{\link{modpar}} and can be adjusted manually or by using \code{\link{change.pnparameters}} to user required specification.
    ##seealso<< \code{\link{SSposnegRichards}}
    ## \code{\link{modpar}}
    }
, ex = function(){
    require(graphics)
    # calculate y from an 8-parameter model
    data(posneg_data)
    modpar(posneg_data$age, posneg_data$mass) #create pnmodelparams for fixed parameters
    y <- posnegRichards_eqn(10, 1000, 0.5, 25, 1, 100, 0.5, 125, 1, modno = 1)
    
    # calculate y from a 4-parameter positive only model, note that all negative parameters specified here are ignored
    # and replaced with value from pnmodelparams
    y <- posnegRichards_eqn(10 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12)
    
    # plot a logistic curve (M=1), note that all negative parameters specified here are ignored
    plot(1:200 ,posnegRichards_eqn(1:200 ,1000 ,0.5 ,25 ,1 ,100 ,0.5 ,125 ,1 ,modno = 12),xlim=c(1, 200),
    xlab = "x", ylab = "y",pch = 1, cex = 0.7)
    }
)
