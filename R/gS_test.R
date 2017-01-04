#' Generalized Joint Location Scale (gJLS) Test
#'
#' This function performs the generalized scale (gS) test (Soave and Sun 2017).  This extension of Levene's (1960) test allows for correlated errors and group uncertainty.
#' @param model a two-sided linear formula object describing the model, with the response (y) on the left and a ~ operator separating the covariates of interest on the right, separated by + operators.
#' @param data a data frame containing the variables named in model and correlation arguments.  This is required.
#' @param correlation an optional corStruct object describing the within-group correlation structure. See the documentation of corClasses for a description of the available corStruct classes. If a grouping variable is to be used, it must be specified in the form argument to the corStruct constructor. Defaults to NULL, corresponding to uncorrelated errors.
#' @keywords gS
#' @export
#' @author David Soave
#' @import quantreg
#' @import nlme
#' @details No missing data are allowed - function will return an "error".  Absolute residuals, are estimated using least absolute deviation (LAD) regression. Outcome (phenotype) must be quantitative and covariate (genotype) may be discrete (categorical) or continuous.
#' @return gS_F the gS test statistic
#' @return numDF the gS test statistic numerator degrees of freedom
#' @return denDF the gS test statistic denominator degrees of freedom
#' @return gS_p the gS test p-value
#' @references Soave, D., Corvol, H., Panjwani, N., Gong, J., Li, W., Boelle, P.Y., Durie, P.R., Paterson, A.D., Rommens, J.M., Strug, L.J., and Sun, L. (2015). A Joint Location-Scale Test Improves Power to Detect Associated SNPs, Gene Sets, and Pathways. American journal of human genetics 97, 125-138.
#' @references Soave, D. and Sun, L. (2017). A Generalized Levene's Scale Test for Variance Heterogeneity in the Presence of Sample Correlation and Group Uncertainty. Biometrics (Accepted).
#' @examples
#' #################################################################################
#' ## Example simulating data from model [i] (Soave et al. 2015 AJHG)
#' #################################################################################
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



#Note that p_S is obtained using a stage 1 median regression (rq function, tau=0.5) where group medians are chosen to be the larger of two middle values when the group size is even.
gS_test <-function(model, data, correlation=NULL){

  ## check if there is missing data
  if(sum(is.na(data)) > 0)  stop("missing value(s) not allowed")
  model <- terms(model)

  ## Obtain the p-value scale test
  lm1 <- rq(model,tau=.5,data=data)
  data$d1 <- abs(resid(lm1))
  model2=reformulate(attr(model, "term.labels"),response="d1")

  if(is.null(correlation)) {
    aovfit <- anova(lm(model2,data=data))
    gS_F<-aovfit[1,4];numDF<-aovfit[1,1];denDF<-aovfit[2,1];gS_p<-aovfit[1,5]
  } else {
    fit<-gls(model2,data=data,correlation= correlation,method="ML",control=lmeControl(opt = "optim"))
    aovfit<-anova(fit,Terms=2:dim(anova(fit))[1])
    gS_F<-aovfit[1,2];numDF<-aovfit[1,1];denDF<-fit$dims$N-fit$dims$p;gS_p<-aovfit[1,3]
  }
  ## return the scale p-value
  return(cbind(gS_F, numDF, denDF, gS_p))
}
