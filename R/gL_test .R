#' Generalized Location (gL) Test
#'
#' This function performs a generalized location (gL) test (Soave and Sun 2017) using the generalized least squares function, gls(), from the nlme package (Pinheiro and Bates 2000) to allow for correlated errors.
#' @param model a two-sided linear formula object describing the model, with the response (y) on the left and a ~ operator separating the covariates of interest on the right, separated by + operators.
#' @param data a data frame containing the variables named in model and correlation arguments.  This is required.
#' @param correlation an optional corStruct object describing the within-group correlation structure. The correlation structure must be called directly from the nlme pacakge using "nlme::" (see examples below). See the documentation of corClasses for a description of the available corStruct classes. If a grouping variable is to be used, it must be specified in the form argument to the corStruct constructor. Defaults to NULL, corresponding to uncorrelated errors.
#' @keywords gL
#' @export
#' @author David Soave
#' @import nlme
#' @details No missing data are allowed - function will return an "error". Outcome (phenotype) must be quantitative and covariate (genotype) may be discrete (categorical) or continuous.
#' @return gL_F the gL test statistic
#' @return numDF the gL test statistic numerator degrees of freedom
#' @return denDF the gL test statistic denominator degrees of freedom
#' @return gL_p the gL test p-value
#' @references Soave, D., Corvol, H., Panjwani, N., Gong, J., Li, W., Boelle, P.Y., Durie, P.R., Paterson, A.D., Rommens, J.M., Strug, L.J., and Sun, L. (2015). A Joint Location-Scale Test Improves Power to Detect Associated SNPs, Gene Sets, and Pathways. American journal of human genetics 97, 125-138.
#' @references Soave, D. and Sun, L. (2017). A Generalized Levene's Scale Test for Variance Heterogeneity in the Presence of Sample Correlation and Group Uncertainty. Biometrics (Accepted).
#' @seealso \code{\link{gS_test}}, \code{\link{gJLS_test}}
#' @examples
#' #################################################################################
#' ## Example simulating data from model [i] (Soave et al. 2015 AJHG)
#' #################################################################################
#'
#'
#' n<-2000  ## total sample size
#' pA<-0.3  ## MAF
#' pE1<-0.3  ## frequency of exposure E1
#'
#' ## Genotypes (XG)
#' genocount<-rmultinom(1,size=n,prob=c(pA*pA, 2*pA*(1-pA), (1-pA)*(1-pA)))
#' XG<-c(rep(0, genocount[1]), rep(1, genocount[2]), rep(2,genocount[3]))
#' XG<-sample(XG,size=length(XG),replace=FALSE)
#'
#' ## Exposures (E1)
#' E1<-rbinom(n,1,prob=pE1)
#'
#' ## Phenotype (y)
#' y<-0.01*XG+0.3*E1+0.5*XG*E1+rnorm(n,0,1)
#'
#'# Additive model (will work with dosages)
#'JLS_test(y,XG,XG)
#'# Genotypic model (will not work with dosages --> factor() will create many groups)
#'JLS_test(y,factor(XG),factor(XG))
#'
#'X2=round(cbind(XG==1,XG==2)) #convert to 2 columns
#'
#'# Genotypic model --> same result as results above using JLS_test(y,factor(X),factor(X))
#'# This is how genotype probabilities will be analyzed (using a 2 column design matrix)
#'JLS_test(y,X2,X2)


gL_test <-function(model, data, correlation=NULL){

  ## check if there is missing data
  if(sum(is.na(data)) > 0)  stop("missing value(s) not allowed")

  ## Obtain the p-value scale test
  if(is.null(correlation)) {
    model=terms(model)
    model0=model[-(1:2)]
    fit0<-lm(model0, data = data)
    fit<-lm(model,data=data)
    aovfit <- anova(fit0,fit)
    gS_F<-aovfit[2,5];numDF<-aovfit[2,3];denDF<-aovfit[2,1];gS_p<-aovfit[2,6]
  } else {
    fit<-gls(model,data=data,correlation=correlation,method="ML",control=lmeControl(opt = "optim"))
    aovfit<-anova(fit,Terms=2:dim(anova(fit))[1])
    gL_F<-aovfit[1,2];numDF<-aovfit[1,1];denDF<-fit$dims$N-fit$dims$p;gL_p<-aovfit[1,3]
  }
  ## return the scale p-value
  return(cbind(gL_F, numDF, denDF, gL_p))
}
