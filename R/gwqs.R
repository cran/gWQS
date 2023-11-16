#' Fitting Weighted Quantile Sum regression models
#'
#' Fits Weighted Quantile Sum (WQS) regression  (Carrico et al. (2014) \doi{10.1007/s13253-014-0180-3}),
#' a random subset implementation of WQS (Curtin et al. (2019) \doi{10.1080/03610918.2019.1577971}),
#' a repeated holdout validation WQS (Tanner et al. (2019) \doi{10.1016/j.mex.2019.11.008}) and a WQS with
#' 2 indices (Renzetti et al. (2023) \doi{10.3389/fpubh.2023.1289579}) for continuous, binomial,
#' multinomial, Poisson, quasi-Poisson and negative binomial outcomes.
#'
#' @param formula An object of class \code{formula} specifying the relationship to be tested. The \code{wqs}
#' term must be included in \code{formula}, e.g. \code{y ~ wqs + ...}. To test for an interaction term with
#' a continuous variable \code{a} or for a quadratic term we can specify the \code{formula} as below:
#' \code{y ~ wqs*a + ...} and \code{y ~ wqs + I(wqs^2) + ...}, respectively.
#' @param data The \code{data.frame} containing the variables to be included in the model.
#' @param na.action \code{\link[stats]{model.frame}}. \code{na.omit} is the default.
#' @param weights An optional term containing the name of the variable in the dataset representing the weights
#' to be used in the fitting process. Should be \code{NULL} or the variable name.
#' @param mix_name A character vector listing the variables contributing to a mixture effect.
#' @param stratified The character name of the variable for which you want to stratify for.
#' It has to be a \code{factor}.
#' @param rh Number of repeated holdout validations.
#' @param b Number of bootstrap samples used in parameter estimation. No bootstrap will be performed if b = 1.
#' @param b1_pos A logical value that determines whether weights are derived from models where the beta
#' values were positive (\code{TRUE}) or negative (\code{FALSE}).
#' @param bint_cont_pos A logical value that determines whether weights are derived from models where the
#' beta parameter of the interaction term between the WQS index and a continuous variable were
#' positive (\code{TRUE}) or negative (\code{FALSE}).
#' @param bint_cat_pos A logical value or a vector of logical values that determines whether weights are
#' derived from models where the slopes of the WQS index for each level (other than the reference one)
#' of the interacting categorical variable were positive (\code{TRUE}) or negative (\code{FALSE}).
#' @param b_constr A logial value that determines whether to apply positive (if \code{b1_pos = TRUE}) or
#' negative (if \code{b1_pos = FALSE}) constraints in the optimization function for the weight estimation.
#' @param zero_infl A logical value (\code{TRUE} or \code{FALSE}) that allows to fit a zero inflated
#' model in case \code{family = "poisson"} or \code{family = "negbin"}.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). If \code{q = NULL} then
#' the values of the mixture variables are taken (these must be standardized).
#' @param validation Percentage of the dataset to be used to validate the model. If
#' \code{validation = 0} then the test dataset is used as validation dataset too.
#' @param validation_rows A list of a single (if rh=1) or multiple vectors containing the rows to be
#' considered in the validation step. When "validation_rows=NULL" (default) the function randomly choose
#' the observations to be considered in the validation step.
#' @param family A character value that allows to decide for the glm: \code{gaussian} for linear regression,
#' \code{binomial} for logistic regression, \code{poisson} for Poisson regression,
#' \code{quasipoisson} for quasi-Poisson regression, \code{"negbin"} for negative binomial regression.
#' @param signal Character identifying the signal function to be used when the average weights
#' are estimated. It can take values from \code{"one"} to apply the identity, \code{"abst"} to apply
#' the absolute value of the t-statistic, \code{"t2"} to apply the squared value of the t-statistic,
#' \code{"expt"} to apply the exponential of the t-statistic as signal function.
#' @param rs A logic value. If \code{rs = FALSE} then the bootstrap implementation of WQS is performed.
#' If \code{rs = TRUE} then the random subset implementation of WQS is applied (see the "Details" and the
#' vignette for further information).
#' @param n_vars The number of mixture components to be included at each random subset step.
#' If \code{rs = TRUE} and \code{n_vars = NULL} then the square root of the number of elements
#' in the mixture is taken.
#' @param zilink Character specification of link function in the binary zero-inflation model
#' (you can choose among \code{"logit", "probit", "cloglog", "cauchit", "log"}).
#' @param seed An \code{integer} value to fix the seed, if it is equal to \code{NULL} no seed is chosen.
#' @param wp,wn An optional set of starting weights for the positive (\code{wp}) and negative (\code{wn}) directions
#' to be passed to the optimization function. The default is \code{wp = NULL, wn = NULL} to let the \code{gwqs} function
#' set the starting values.
#' @param lambda The value of the penalization term used to shrink towards 0 the weights that are not
#' truly associated with the outcome (see the "Details" and the vignette for further information).
#' @param plan_strategy A character value that allows to choose the evaluation strategies for the
#' \code{plan} function. You can choose among "sequential", "transparent", "multisession", "multicore",
#' "multiprocess", "cluster" and "remote" (see \code{\link[future]{plan}} help page for more details).
#' @param optim.method A character identifying the method to be used by the \code{\link[stats]{optim}} function
#' (you can choose among \code{"BFGS", "Nelder-Mead", "CG", "SANN"}, \code{"BFGS"} is the default).
#' See \code{\link[stats]{optim}} for details.
#' @param control The control list of optimization parameters. See \code{\link[stats]{optim}} for details.
#' @param b1_constr The argument is deprecated, use 'b_constr' instead.
#' @param ... Additional arguments to be passed to the function
#'
#' @details
#' \code{gWQS} uses the \code{glm} function in the \bold{stats} package to fit the linear, logistic,
#' the Poisson and the quasi-Poisson regression, while the \code{glm.nb} function from the \bold{MASS}
#' package is used to fit the negative binomial regression respectively. The \code{nlm} function from
#' the \bold{stats} package was used to optimize the log-likelihood of the multinomial regression.\cr
#'
#' The \code{\link[stats]{optim}} optimization function is used to estimate the weights at each
#' bootstrap step.\cr
#'
#' The \code{seed} argument specifies a fixed seed through the \code{\link[base]{set.seed}} function.\cr
#'
#' The \code{rs} term allows to choose the type of methodology between the bootstrap implementation
#' (WQSBS) or the random subset implementation (WQSRS) of the WQS. The first method performs \code{b}
#' bootstrapped samples to estimate the weights while the second creates \code{b} randomly-selected
#' subset of the total predictor set. For further details please see the vignette
#' ("How to use gWQS package") and the references below.
#'
#' @return \code{gwqs} return the results of the WQS regression as well as many other objects and datasets.
#'
#' \item{fit}{The object that summarizes the output of the WQS model, reflecting a
#' linear, logistic, multinomial, Poisson, quasi-Poisson or negative binomial regression
#' depending on how the \code{family} parameter was specified.
#' The summary function can be used to call and print fit data (not for multinomial regression).}
#' \item{final_weights}{\code{data.frame} containing the final weights associated to each chemical.}
#' \item{conv}{Indicates whether the solver has converged (0) or not (1 or 2).}
#' \item{bres}{Matrix of estimated weights, mixture effect parameter estimates and the associated
#' standard errors, statistics and p-values estimated for each bootstrap iteration.}
#' \item{wqs}{Vector containing the wqs index for each subject.}
#' \item{pwqs}{Vector containing the positive wqs index for each subject.}
#' \item{nwqs}{Vector containing the negative wqs index for each subject.}
#' \item{qi}{List of the cutoffs used to divide in quantiles the variables in the mixture}
#' \item{bindex}{List of vectors containing the \code{rownames} of the subjects included in each
#' bootstrap dataset.}
#' \item{y_wqs_df}{\code{data.frame} containing the dependent variable values adjusted for the
#' residuals of a fitted model adjusted for covariates (original values when \code{family = binomial}
#' or \code{"multinomial"}) and the wqs index estimated values.}
#' \item{family}{The family specified.}
#' \item{call}{The matched call.}
#' \item{formula}{The formula supplied.}
#' \item{mix_name}{The vector of variable names used to identify the elements in the mixture.}
#' \item{q}{The method used to rank varibales included in the mixture.}
#' \item{n_levels}{The number of levels of the of the dependent variable when a multinomial regression is ran.}
#' \item{zero_infl}{If a zero inflated model was ran (\code{TRUE}) or not (\code{FALE})}
#' \item{zilink}{The chosen link function when a zero inflated model was ran.}
#' \item{dwqs}{A logical value whether two indices were included (\code{TRUE}) in the model or not (\code{FALSE}).}
#' \item{levelnames}{The name of each level when a multinomial regression is ran.}
#' \item{data}{The data used in the WQS analysis.}
#' \item{objfn_values}{The vector of the b values of the objective function corresponding to the optima values}
#' \item{optim_messages}{The vector of character strings giving any additional information returned by the
#' optimizer, or NULL.}
#' \item{gwqslist}{List of the output from the \code{rh} WQS models.}
#' \item{coefmat}{Matrix containing the parameter estimates from each repeated holdout WQS model.}
#' \item{wmat}{Matrix containing the weight estimates from each repeated holdout WQS model.}
#' \item{rh}{The number of repeated holdout performed.}
#'
#' @author
#' Stefano Renzetti, Paul Curtin, Allan C Just, Ghalib Bello, Chris Gennings
#'
#' @references
#' Carrico C, Gennings C, Wheeler D, Factor-Litvak P. Characterization of a weighted quantile sum
#' regression for highly correlated data in a risk analysis setting. J Biol Agricul Environ Stat.
#' 2014:1-21. ISSN: 1085-7117. \doi{10.1007/s13253-014-0180-3}.\cr
#'
#' Curtin P, Kellogg J, Cech N, Gennings C (2021). A random subset implementation of weighted quantile
#' sum (WQSRS) regression for analysis of high-dimensional mixtures, Communications in Statistics -
#' Simulation and Computation, 50:4, 1119-1134. \doi{10.1080/03610918.2019.1577971}.\cr
#'
#' Tanner EM, Bornehag CG, Gennings C. Repeated holdout validation for weighted quantile sum regression.
#' MethodsX. 2019 Nov 22;6:2855-2860. \doi{10.1016/j.mex.2019.11.008}. PMID: 31871919; PMCID: PMC6911906.\cr
#'
#' Renzetti S, Gennings C and Calza S (2023) A weighted quantile sum regression with penalized weights
#' and two indices. Front Public Health 11:1151821. \doi{10.3389/fpubh.2023.1151821}.\cr
#'
#' @seealso \link[stats]{glm}, \link[MASS]{glm.nb}, \link[nnet]{multinom}, \link[pscl]{zeroinfl}.
#'
#' @examples
#' # we save the names of the mixture variables in the variable "toxic_chems"
#' toxic_chems = names(wqs_data)[1:34]
#'
#' # To run a linear model and save the results in the variable "results". This linear model
#' # (family = gaussian) will rank/standardize variables in quartiles (q = 4), perform a
#' # 40/60 split of the data for training/validation (validation = 0.6), and estimate weights
#' # over 2 bootstrap samples (b = 2; in practical applications at least 100 bootstraps
#' # should be used). Weights will be derived from mixture effect parameters that are positive
#' # (b1_pos = TRUE). A unique seed was specified (seed = 2016) so this model will be
#' # reproducible, and plots describing the variable weights and linear relationship will be
#' # generated as output (plots = TRUE). In the end tables describing the weights values and
#' # the model parameters with the respectively statistics are generated in the plots window
#' # (tables = TRUE):
#' results = gwqs(yLBX ~ wqs, mix_name = toxic_chems, data = wqs_data, q = 4, validation = 0.6,
#'                b = 2, b1_pos = TRUE, b_constr = FALSE, family = gaussian, seed = 2016)
#'
#' # to test the significance of the covariates
#' summary(results)
#'
#' @import ggplot2
#' @import MASS
#' @import stats
#' @import knitr
#' @import kableExtra
#'
#' @importFrom broom augment
#' @importFrom rlist list.cbind
#' @importFrom reshape2 melt
#' @importFrom plotROC geom_roc style_roc calc_auc
#' @importFrom nnet multinom
#' @importFrom future plan future value
#' @importFrom future.apply future_lapply
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom pscl zeroinfl
#' @importFrom car vif
#' @importFrom Matrix bdiag
#' @importFrom utils flush.console
#' @importFrom bookdown html_document2
#'
#' @export gwqs
#' @rdname main_functions
#'
#' @usage  gwqs(formula, data, na.action, weights, mix_name, stratified, rh = 1, b = 100,
#'              b1_pos = TRUE, bint_cont_pos = NULL, bint_cat_pos = NULL, b_constr = FALSE,
#'              zero_infl = FALSE, q = 4, validation = 0.6, validation_rows = NULL,
#'              family = gaussian, signal = c("t2", "t3", "one", "abst", "expt"),
#'              rs = FALSE, n_vars = NULL,
#'              zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
#'              wp = NULL, wn = NULL, plan_strategy = "sequential", lambda = 0,
#'              optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
#'              control = list(trace = FALSE, maxit = 2000, reltol = 1e-9),
#'              b1_constr = NULL, ...)

gwqs <- function(formula, data, na.action, weights, mix_name, stratified, rh = 1, b = 100,
                 b1_pos = TRUE, bint_cont_pos = NULL, bint_cat_pos = NULL, b_constr = FALSE,
                 zero_infl = FALSE, q = 4, validation = 0.6, validation_rows = NULL,
                 family = gaussian, signal = c("t2", "t3", "one", "abst", "expt"),
                 rs = FALSE, n_vars = NULL,
                 zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
                 wp = NULL, wn = NULL, plan_strategy = "sequential", lambda = 0,
                 optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
                 control = list(trace = FALSE, maxit = 2000, reltol = 1e-9),
                 b1_constr = NULL, ...){

  wqsGaussBin <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, level_names, wqsvars, pwqsvars,
                          nwqsvars, family, dwqs, zilink, zero_infl, formula, ff, bwghts, stratified, stratlev,
                          b1_pos, b_constr, lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

    if(b_constr){
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
      if(wqsint){
        if(intcont) initp[intvars] <- sapply(1:length(intvars), function(i) ifelse(bint_cont_pos[i], abs(initp[intvars[i]]), -abs(initp[intvars[i]])))
        else initp[intvars] <- sapply(1:length(intvars), function(i) ifelse(bint_cat_pos[i] & initp[intvars[i]]-initp["wqs"]<0, abs(initp[intvars[i]]),
                                                                            ifelse(!bint_cat_pos[i] & initp[intvars[i]]-initp["wqs"]>0, -abs(initp[intvars[i]]), abs(initp[intvars[i]]))))
      }
      # if(!is.null(stratified)) initp[grepl("wqs", Xnames)] <- ifelse(b1_pos, 1, -1)*abs(initp["wqs"] - abs(initp["wqs"] + initp[grepl("wqs", Xnames)]))
    }
    w <- initp[(kx + 1):(kx + kw)]
    pen <- sum(abs(w))
    # w <- abs(w)
    w <- w^2
    if(!is.null(stratified)){
      w <- matrix(w, kw/stratlev, stratlev)
      w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
      w <- as.numeric(w)
    }
    else w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ%*%w)
    if(wqsint) bXm[,intpos] <- bXm[,"wqs"]*bXm[,intvars]
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm%*%b_covs) + boffset
    f = sum(family$dev.resids(y = bY, mu = family$linkinv(term), wt = bwghts)) + lambda*pen

    return(f)
  }

  wqsGaussBin_dwqs <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, level_names, wqsvars, pwqsvars,
                               nwqsvars, family, dwqs, zilink, zero_infl, formula, ff, bwghts, stratified, stratlev,
                               b1_pos, b_constr, lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

    initp["pwqs"] <- abs(initp["pwqs"])
    initp["nwqs"] <- -abs(initp["nwqs"])
    wp <- initp[(kx + 1):(kx + kw)]
    wn <- initp[(kx + kw + 1):(kx + 2*kw)]
    pen <- sum(abs(wp)) + sum(abs(wn))
    # wp <- abs(wp)/sum(abs(wp))
    # wn <- abs(wn)/sum(abs(wn))
    wp <- (wp^2)/sum(wp^2)
    wn <- (wn^2)/sum(wn^2)
    bXm[, "pwqs"] <- as.numeric(bQ%*%wp)
    bXm[, "nwqs"] <- as.numeric(bQ%*%wn)
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm%*%b_covs) + boffset
    f = sum(family$dev.resids(y = bY, mu = family$linkinv(term), wt = bwghts)) + lambda*pen

    return(f)
  }

  wqsPoisson <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars,
                         nwqsvars, family, dwqs, zilink, zero_infl, formula, ff, bwghts, stratified, stratlev,
                         b1_pos, b_constr, lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

    if(b_constr){
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
      if(wqsint){
        if(intcont) initp[intvars] <- sapply(1:length(intvars), function(i) ifelse(bint_cont_pos[i], abs(initp[intvars[i]]), -abs(initp[intvars[i]])))
        else initp[intvars] <- sapply(1:length(intvars), function(i) ifelse(bint_cat_pos[i] & initp[intvars[i]]-initp["wqs"]<0, abs(initp[intvars[i]]),
                                                                            ifelse(!bint_cat_pos[i] & initp[intvars[i]]-initp["wqs"]>0, -abs(initp[intvars[i]]), abs(initp[intvars[i]]))))
      }
    }
    w <- initp[(kx + 1):(kx + kw)]
    pen <- sum(abs(w))
    w <- w^2
    if(!is.null(stratified)){
      w <- matrix(w, kw/stratlev, stratlev)
      w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
      w <- as.numeric(w)
    }
    else w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ%*%w)
    if(wqsint) bXm[,intpos] <- bXm[,"wqs"]*bXm[,intvars]
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm%*%b_covs) + boffset
    f = -sum(dpois(bY, lambda = exp(term), log = TRUE)*bwghts) + lambda*pen

    return(f)
  }

  wqsPoisson_dwqs <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars,
                              nwqsvars, family, dwqs, zilink, zero_infl, formula, ff, bwghts, stratified, stratlev,
                              b1_pos, b_constr, lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

    initp["pwqs"] <- abs(initp["pwqs"])
    initp["nwqs"] <- -abs(initp["nwqs"])
    wp <- initp[(kx + 1):(kx + kw)]^2
    wn <- initp[(kx + kw + 1):(kx + 2*kw)]^2
    pen <- sum(wp) + sum(wn)
    wp <- wp/sum(wp)
    wn <- wn/sum(wn)
    bXm[, "pwqs"] <- as.numeric(bQ%*%wp)
    bXm[, "nwqs"] <- as.numeric(bQ%*%wn)
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm%*%b_covs) + boffset
    f = -sum(dpois(bY, lambda = exp(term), log = TRUE)*bwghts) + lambda*pen

    return(f)
  }

  wqsNegBin <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars,
                        nwqsvars, family, dwqs, zilink, zero_infl, formula, ff, bwghts, stratified, stratlev,
                        b1_pos, b_constr, lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

    if(b_constr){
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
      if(wqsint){
        if(intcont) initp[intvars] <- sapply(1:length(intvars), function(i) ifelse(bint_cont_pos[i], abs(initp[intvars[i]]), -abs(initp[intvars[i]])))
        else initp[intvars] <- sapply(1:length(intvars), function(i) ifelse(bint_cat_pos[i] & initp[intvars[i]]-initp["wqs"]<0, abs(initp[intvars[i]]),
                                                                            ifelse(!bint_cat_pos[i] & initp[intvars[i]]-initp["wqs"]>0, -abs(initp[intvars[i]]), abs(initp[intvars[i]]))))
      }
    }
    w <- initp[(kx + 1):(kx + kw)]
    pen <- sum(abs(w))
    w <- w^2
    if(!is.null(stratified)){
      w <- matrix(w, kw/stratlev, stratlev)
      w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
      w <- as.numeric(w)
    }
    else w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ%*%w)
    if(wqsint) bXm[,intpos] <- bXm[,"wqs"]*bXm[,intvars]
    b_covs <- initp[1:kx]
    theta <- exp(initp[length(initp)])
    term <- as.numeric(bXm%*%b_covs) + boffset
    f = -sum((suppressWarnings(dnbinom(bY, size = theta, mu = exp(term), log = TRUE)))*bwghts) + lambda*pen

    return(f)
  }

  wqsNegBin_dwqs <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars,
                             nwqsvars, family, dwqs, zilink, zero_infl, formula, ff, bwghts, stratified, stratlev,
                             b1_pos, b_constr, lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

    initp["pwqs"] <- abs(initp["pwqs"])
    initp["nwqs"] <- -abs(initp["nwqs"])
    wp <- initp[(kx + 1):(kx + kw)]^2
    wn <- initp[(kx + kw + 1):(kx + 2*kw)]^2
    pen <- sum(wp) + sum(wn)
    wp <- wp/sum(wp)
    wn <- wn/sum(wn)
    bXm[, "pwqs"] <- as.numeric(bQ%*%wp)
    bXm[, "nwqs"] <- as.numeric(bQ%*%wn)
    b_covs <- initp[1:kx]
    theta <- exp(initp[length(initp)])
    term <- as.numeric(bXm%*%b_covs) + boffset
    f = -sum((suppressWarnings(dnbinom(bY, size = theta, mu = exp(term), log = TRUE)))*bwghts) + lambda*pen

    return(f)
  }

  if(!is.null(b1_constr)){
    b_constr <- b1_constr
    warning("The argument 'b1_constr' is deprecated, use 'b_constr' instead.")
  }

  one <- function(x) rep(1, length(x))

  t2 <- function(x) x^2

  t3 <- function(x) x^3

  expt <- function(x) exp(abs(x))-1

  if(is.character(family)){
    if(family == "negbin") family <- list(family = family)
    else if(family == "multinomial") stop("'family = multinomial' is not supported by the 'gwqs' function. Please use 'gwqs_multinom' function.\n")
    else family <- get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized\n")
  }

  wqs_in_formula <- any(grepl("wqs", rownames(attr(terms(formula), "factors"))))
  pwqs_in_formula <- any(grepl("pwqs", rownames(attr(terms(formula), "factors"))))
  nwqs_in_formula <- any(grepl("nwqs", rownames(attr(terms(formula), "factors"))))
  if(!wqs_in_formula & (!pwqs_in_formula | !nwqs_in_formula))
    stop("'formula' must contain either 'wqs' term or 'pwqs' and 'nwqs' terms: e.g. y ~ wqs + ...\n y ~ pwqs + nwqs + ...\n")

  dwqs <- pwqs_in_formula & nwqs_in_formula

  if(dwqs & rs) stop("The 2iWQS does not support yet the random subset method\n")

  if(dwqs) objfn <- switch(family$family,
                           "gaussian" = wqsGaussBin_dwqs,
                           "binomial" = wqsGaussBin_dwqs,
                           "poisson" = wqsPoisson_dwqs,
                           "quasipoisson" = wqsPoisson_dwqs,
                           "negbin" = wqsNegBin_dwqs)
  else objfn <- switch(family$family,
                       "gaussian" = wqsGaussBin,
                       "binomial" = wqsGaussBin,
                       "poisson" = wqsPoisson,
                       "quasipoisson" = wqsPoisson,
                       "negbin" = wqsNegBin)

  optim.method <- match.arg(optim.method)

  signal <- match.arg(signal)
  signal <- switch(signal,
                  "one" = one,
                  "abst" = abs,
                  "t2" = t2,
                  "t3" = t3,
                  "expt" = expt)

  if(zero_infl){
    zilink <- make.link(match.arg(zilink))
    ff = formula
    if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
      formula = as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))
  }
  else zilink <- ff <- NULL

  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", "weights", "stratified"), names(mc), 0)
  dtf <- mc[c(1, m)]
  dtf[[2]] <- data
  dtf[[1]] <- as.name("selectdatavars")
  dtf <- eval(dtf, parent.frame())

  l <- list(...)
  solve_dir_issue <- ifelse(is.null(l$solve_dir_issue), FALSE, l$solve_dir_issue)

  # defining quantile variables
  if (is.null(q)) {Q = as.matrix(dtf[, mix_name]); qi = NULL}
  else {
    q_f = gwqs_rank(dtf, mix_name, q)
    Q = q_f$Q
    qi = q_f$qi
  }

  # create stratified variables if stratified is not NULL
  m <- match(c("stratified"), names(mc), 0)
  if(!is.null(l$stratified_rh)) m <- 0
  if (m){
    strtfd_out <- stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name <- strtfd_out$mix_name
  }
  else stratified <- NULL
  if(!is.null(l$stratified_rh)) stratified <- l$stratified_rh
  stratlev <- ifelse(is.null(stratified), 1, nlevels(unlist(dtf[, stratified])))
  if(dwqs & !is.null(stratified)) stop("The 2iWQS does not support yet the presence of stratified weights\n")

  if(!is.null(seed) & !is.numeric(seed)) stop("seed must be numeric or NULL\n")
  if(!is.null(seed)) set.seed(seed)
  fseed <- TRUE

  N <- nrow(dtf)

  # splitting the dataset
  if(!is.numeric(validation) | validation < 0 | validation >= 1) stop("'validation' must be numeric >= 0 and < 1\n")
  if(is.null(validation_rows)) validation_rows <- lapply(1:rh, function(i) sample(c(F, T), N, replace = T, prob = c(1-validation, validation)))
  else if(length(validation_rows) != rh | !any(unique(unlist(validation_rows)) %in% c(F, T))) stop("'validation_rows' must be a list of length 'rh' of logical vectors\n")

  # estimate starting weights
  if(dwqs){
    mc_pwqs <- mc_nwqs <- mc
    mc_pwqs$data <- mc_nwqs$data <- data
    mc_pwqs$mix_name <- mc_nwqs$mix_name <- mix_name
    mc_pwqs$validation <- mc_nwqs$validation <- 0
    mc_pwqs$plan_strategy <- mc_nwqs$plan_strategy <- "sequential"
    mc_pwqs$b_constr <- mc_nwqs$b_constr <- TRUE
    mc_pwqs$b1_pos <- TRUE
    mc_nwqs$b1_pos <- FALSE
    mc_pwqs$validation_rows <- mc_nwqs$validation_rows <- NULL
    if(rs) mc_pwqs$b <- mc_nwqs$b <- 100
    else mc_pwqs$b <- mc_nwqs$b <- 1
    mc_pwqs$rh <- mc_nwqs$rh <- 1
    mc_pwqs$formula <- mc_nwqs$formula <- as.formula(gsub("pwqs", "wqs", deparse(formula)))
    mc_pwqs$formula <- mc_nwqs$formula <- remove_terms(mc_pwqs$formula, "nwqs")
    mc_pwqs$solve_dir_issue <- mc_nwqs$solve_dir_issue <- "inverse"

    if(is.null(wp)){
      pwqs0 <- tryCatch(eval(mc_pwqs), error = function(e) NULL)
      if(is.null(pwqs0)) wp <- rep(1/length(mix_name), length(mix_name))
      else wp <- pwqs0$final_weights$mean_weight[match(mix_name, pwqs0$final_weights$mix_name)]
    }
    if(is.null(wn)){
      nwqs0 <- tryCatch(eval(mc_nwqs), error = function(e) NULL)
      if(is.null(nwqs0)) wn <- rep(1/length(mix_name), length(mix_name))
      else wn <- nwqs0$final_weights$mean_weight[match(mix_name, nwqs0$final_weights$mix_name)]
    }
  }

  if(is.null(n_vars)) n_vars = round(sqrt(length(mix_name)))

  # parameters estimation and model fitting
  m <- match(c("weights"), names(mc), 0)
  if(m[1]) dtf$wghts <- wghts <- unlist(dtf[, weights, drop = FALSE])
  else  dtf$wghts <- wghts <- rep(1, N)

  if (family$family %in% c("gaussian", "quasipoisson")) ts = "t"
  else if (family$family %in% c("binomial", "poisson", "multinomial", "negbin")) ts = "z"

  if(!is.numeric(b)) stop("'b' must be a number\n")


  Xnames <- parnames(dtf, formula, NULL)
  kx <- length(Xnames)
  kw <- ncol(Q)
  initp <- values.0(kw, Xnames, kx, formula, ff, wghts, dtf, stratified, stratlev, b1_pos, family, wp, wn, dwqs, zilink, zero_infl)
  mf <- model.frame(formula, dtf)
  Y <- model.response(mf, "any")
  if(family$family == "binomial" & any(class(Y) %in% c("factor", "character"))){
    if(inherits(Y, "character")) Y = factor(Y)
    Y <- as.numeric(Y != levels(Y)[1])
  }
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- rep(0, nrow(dtf))

  Xm <- model.matrix(formula, dtf)

  intpos <- grepl("wqs", colnames(Xm)) & grepl(":", colnames(Xm))
  intvarname <- colnames(Xm)[intpos]
  wqsint <- any(intpos)
  intvars <- gsub(":wqs", "", gsub("wqs:", "", colnames(Xm)[intpos]))
  if(wqsint){
    intvar0 <- gsub(":wqs", "", gsub("wqs:", "", attr(terms(formula), "term.labels")[grepl("wqs:", attr(terms(formula), "term.labels")) | grepl(":wqs", attr(terms(formula), "term.labels"))]))
    intcont <- is.numeric(dtf[,intvar0])
    if(wqsint & ((is.null(bint_cont_pos) & intcont) | (is.null(bint_cat_pos) & !intcont))) stop(paste0("Please specify the direction of the interaction term(s) thourgh the parameter ", ifelse(intcont, "bint_cont_pos\n", "bint_cat_pos\n")))
  }
  else bint_cont_pos <- bint_cat_pos <- intcont <- NULL

  if(dwqs & wqsint) stop("The 2iWQS does not support yet the presence of an interaction between a covariate and the WQS index\n")

  plan_strategy_rh <- ifelse(rh==1, "sequential", plan_strategy)
  plan_strategy_b <- ifelse(rh==1, plan_strategy, "sequential")

  plan(plan_strategy_rh)
  gwqslist <- future_lapply(X=1:rh, FUN = function(i){
  # gwqslist <- lapply(1:rh, function(i){
    if(control$trace) cat("start opt\n")
    plan(plan_strategy_b)
    param <- future_lapply(X = 1:b, FUN = optim.f, objfn = objfn, Y = Y[!validation_rows[[i]]],
    # param <- lapply(X = 1:b, FUN = optim.f, objfn = objfn, Y = Y[!validation_rows[[i]]],
                    Xm = Xm[!validation_rows[[i]],], Q = Q[!validation_rows[[i]],], offset = offset[!validation_rows[[i]]], wghts = wghts[!validation_rows[[i]]],
                    initp = initp, b1_pos = b1_pos, b_constr = b_constr, n_vars = n_vars, dwqs = dwqs, family = family, rs = rs,
                    zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff, kx = kx, kw = kw,
                    Xnames = Xnames, stratified = stratified, b = b, optim.method = optim.method, control = control,
                    lambda = lambda, stratlev = stratlev, intpos = intpos, wqsint = wqsint, intvars = intvars, bint_cont_pos = bint_cont_pos, bint_cat_pos = bint_cat_pos, intcont = intcont, future.seed = fseed)

    conv <- c(sapply(param, function(i) i$conv))
    counts <- c(sapply(param, function(i) i$counts))
    val <- c(sapply(param, function(i) i$val))
    mex <- lapply(param, function(i) i$mex)
    bindex <- lapply(param, function(i) i$bindex)
    slctd_vars <- lapply(param, function(i) i$slctd_vars)

    if(rs){
      plan(plan_strategy)
      param <- future_lapply(X = 1:b, FUN = set_par_names, slctd_vars, param, q_name = colnames(Q), family = family,
                             dwqs = dwqs, future.seed = FALSE)
    }


    ### NEED TO CREATE TOLERANCE FOR MULTINOMIAL REGRESSION FOR DWQS AND CHECK FOR ZERO_INFL

    build_bres <- function(wqs_var_name, interval){
        wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt[interval]))
        b1 <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, "Estimate"])
        if(wqsint){
          bint <- as.data.frame(do.call("rbind", lapply(param, function(i) summary(i$mfit$m_f)$coefficients[intpos, "Estimate"])))
          if(!intcont) bint <- b1+bint
          names(bint) <- intvarname
        }
        se <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, "Std. Error"])
        stat <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, paste0(ts, " value")])
        p_val <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, gsub("x", ts, "Pr(>|x|)")])
        if(dwqs){
          tol <- sapply(param, function(i){
            tmp <- vif(i$mfit$m_f)
            if("matrix" %in% class(tmp)) tmp <- tmp[,1]
            1/tmp[wqs_var_name]
          })
        }

      n_non_conv = sum(conv == 1)
      if(n_non_conv == 0 & control$trace) cat(paste0("The optimization function always converged\n"))
      else if(n_non_conv == b & rh == 1) stop("The optimization function never converged\n")
      else if(n_non_conv == b & rh > 1) return(NULL)
      else if(control$trace) cat(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))
      if(control$trace) cat(paste0("There are ", ifelse(b1_pos, sum(b1 >= 0, na.rm = T), sum(b1 <= 0, na.rm = T)),
                                   ifelse(b1_pos, " positive", " negative"), " bootstrapped b1 out of ", b, "\n"))
      bres <- as.data.frame(cbind(wght_matrix, b1, se, stat, p_val))
      names(bres) <- c(colnames(Q), "b1", "Std_Error", "stat", "p_val")
      if(wqsint) bres <- cbind(bres, bint)
      if(dwqs) bres$tol <- tol
      strata_names <- NULL

      return(bres)
    }

    if(dwqs){
      bresp <- build_bres("pwqs", 1:ncol(Q))
      bresn <- build_bres("nwqs", (ncol(Q)+1):(2*ncol(Q)))
      bres <- list(bresp = bresp, bresn = bresn)
    }
    else{
      bresp <- bresn <- NULL
      bres <- build_bres("wqs", 1:ncol(Q))
    }

    if(is.null(bres)) return(NULL)

    mean_weight_f <- function(bres, b1_pos, par){
        if(rs) bres[mix_name][is.na(bres[mix_name])] <- 0
        filter_direction <- apply(as.matrix(bres[,c("b1", intvarname)]>0), 1, function(i) all(i == c(b1_pos, ifelse(intcont, bint_cont_pos, bint_cat_pos))))
        mean_weight = apply(bres[filter_direction & conv == 0, mix_name], 2, weighted.mean, signal(bres[filter_direction & conv == 0, par]))
        if(all(is.nan(mean_weight))){
          if(dwqs){
            Sys.sleep(0.01)
            message(paste0("There are no ", ifelse(b1_pos, "positive", "negative"), " b1 and/or bint values in the specified direction(s) among the bootstrapped models.\nRunning the model in a single ", ifelse(!b1_pos, "positive", "negative"), " direction...\n"))
            flush.console()
            mc$formula <- as.formula(gsub(ifelse(b1_pos, "nwqs", "pwqs"), "wqs", deparse(formula)))
            mc$formula <- remove_terms(mc$formula, ifelse(b1_pos, "pwqs", "nwqs"))
            mc$b1_pos <- !b1_pos
            single_wqs <- eval(mc)
            return(single_wqs)
          }
          else if(solve_dir_issue == "inverse") mean_weight <- 1/apply(bres[, mix_name], 2, weighted.mean, signal(bres[, par]))/sum(1/apply(bres[, mix_name], 2, weighted.mean, signal(bres[, par])))
          else if(solve_dir_issue == "average") mean_weight <- rep(1/length(mix_name), length(mix_name))
          else if(!(solve_dir_issue %in% c("average", "inverse")) & rh == 1) stop("There are no b1 or bint values in the specified direction(s) among the bootstrapped models\n")
          else if(!(solve_dir_issue %in% c("average", "inverse")) & rh > 1) return(NULL)
        }

      return(mean_weight)
    }

    if(dwqs){
      mean_weight_p <- mean_weight_f(bresp, b1_pos = TRUE, par = "tol")
      if(inherits(mean_weight_p, "gwqs")) return(mean_weight_p)
      mean_weight_n <- mean_weight_f(bresn, b1_pos = FALSE, par = "tol")
      if(inherits(mean_weight_n, "gwqs")) return(mean_weight_n)
      mean_weight <- c(mean_weight_p, mean_weight_n)
    }
    else mean_weight <- mean_weight_f(bres, b1_pos = b1_pos, par = "stat")

    if(is.null(mean_weight)) return(NULL)

      # fit the final model with the estimated weights
      if(validation==0) validation_rows <- lapply(validation_rows, function(i) !i)
      wqs_model = model.fit(mean_weight, dtf[validation_rows[[i]],], Q[validation_rows[[i]],], #if(family$family == "multinomial") Y[validation_rows[[i]],] else
                            Y[validation_rows[[i]]],
                            family, dwqs, zilink, formula, ff, wghts[validation_rows[[i]]], offset[validation_rows[[i]]], initp, Xnames,
                            stratified, b1_pos, zero_infl, kx, kw, intpos, wqsint, intvars)

      # Plots
      if(dwqs){
        data_plot <- data.frame(mix_name, mean_weight_p, mean_weight_n, stringsAsFactors = TRUE)
        data_plot = data_plot[order(data_plot$mean_weight_p, decreasing = TRUE),]
        pwqs_index = as.numeric(unlist(wqs_model$pwqs))
        nwqs_index = as.numeric(unlist(wqs_model$nwqs))
        dtf$pwqs <- as.numeric(Q%*%mean_weight_p)
        dtf$nwqs <- as.numeric(Q%*%mean_weight_n)
        dtf$wqs <- wqs_index <- NULL
      }
      else{
        data_plot <- data.frame(mix_name, mean_weight, stringsAsFactors = TRUE)
        data_plot <- data_plot[order(data_plot$mean_weight, decreasing = TRUE),]
        wqs_index <- as.numeric(unlist(wqs_model$wqs))
        dtf$wqs <- as.numeric(Q%*%mean_weight)
        form_terms <- attr(terms.formula(formula), "term.labels")
        form_terms <- gsub("^`|`$", "", form_terms)
        wqsint <- any(grepl("wqs:", form_terms))
        intwqs <- any(grepl(":wqs", form_terms))
        intpos <- ifelse(wqsint, which(grepl("wqs:", form_terms)),
                         ifelse(intwqs, which(grepl(":wqs", form_terms)), NA))
        if(wqsint | intwqs){
          int_vars <- unlist(strsplit(form_terms[intpos], ":"))
          intvar <- int_vars[int_vars != "wqs"]
        }
        dtf$pwqs <- dtf$nwqs <- pwqs_index <- nwqs_index <- NULL
      }

      if(rh == 1){
        if(all(attr(terms(formula), "term.labels") %in% "wqs") | all(attr(terms(formula), "term.labels") %in% c("pwqs", "nwqs"))) y_plot <- model.response(model.frame(formula, dtf[validation_rows[[i]],]), "any")
      else{
        if(dwqs){
          formula_wo_wqs <- remove_terms(formula, "pwqs")
          formula_wo_wqs <- remove_terms(formula_wo_wqs, "nwqs")
        }
        else formula_wo_wqs <- remove_terms(formula, "wqs")
        if(zero_infl){
          if(length(ff[[3]]) > 1 && identical(ff[[3]][[1]], as.name("|"))){
            if(all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))), "term.labels")))) f1 <- as.formula(paste0(ff[[2]], " ~ ", 1))
            else f1 <- remove_terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]]))), "wqs")
            if(all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[3]])))), "term.labels")))) f2 <- as.formula(paste0(ff[[2]], " ~ ", 1))
            else f2 <- remove_terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[3]]))), "wqs")
            formula_wo_wqs <- as.formula(paste0(deparse(f1), " | ", f2[[3]]))
          }
          fit <- zeroinfl(formula_wo_wqs, dtf[validation_rows[[i]],], dist = family$family, link = zilink$name)
        }
        else{
          if(family$family == "negbin") fit = glm.nb(formula_wo_wqs, dtf[validation_rows[[i]],])
          else fit = glm(formula_wo_wqs, dtf[validation_rows[[i]],], family = family)
        }
        if(family$family == "binomial") y_plot = fit$y
        else y_plot = mean(fit$y) + resid(fit, type = "pearson")
      }

      if(dwqs){
        y_adj_wqs_df = as.data.frame(cbind(y_plot, wqs_model$pwqs, wqs_model$nwqs))
        names(y_adj_wqs_df) = c(ifelse(family$family == "binomial", "y", "y_adj"), "pwqs", "nwqs")
        y_adj_wqs_df <- melt(y_adj_wqs_df, measure.vars = c("pwqs", "nwqs"))
      }
      else{
        y_adj_wqs_df <- as.data.frame(cbind(y_plot, wqs_model$wqs))
        names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", "y", "y_adj"), "wqs")
      }
    }

    # creating the list of elements to return
    results = list(fit = wqs_model$m_f, final_weights = data_plot, conv = conv, bres = bres, wqs = wqs_index,
                   pwqs = pwqs_index, nwqs = nwqs_index, qi = qi,
                   bindex = bindex, slctd_vars = slctd_vars, validation_rows = validation_rows[[i]], vindex = validation_rows[[i]],
                   y_wqs_df = if(rh==1) y_adj_wqs_df else NULL, family = family, call = cl, formula = formula, mix_name = mix_name,
                   stratified = stratified, q = q, zero_infl = zero_infl, zilink = zilink,
                   dwqs = dwqs, data = dtf, objfn_values = val, optim_messages = mex, rh = rh)
    if(zero_infl) results$formula <- ff
    class(results) <- "gwqs"

    return(results)
  }, future.seed = fseed)

  if(rh == 1) return(gwqslist[[1]])

  gwqslist <- gwqslist[!vapply(gwqslist, is.null, logical(1))]

  if(length(gwqslist) == 0) stop("The model never converged. Try to increase the number of repeated holdout 'rh' or to change the direction of the association.\n")
  if(length(gwqslist) == 1) stop("The model converged for only a single iteration. Try to increase the number of repeated holdout 'rh'\n")

  if(control$trace) cat(paste0("The model converged ", length(gwqslist), " times out of ", rh, " repeated holdout iterations.\n"))

  coeflist <- lapply(gwqslist, function(i){
    if(zero_infl) i$fit$coefficients$count
    else i$fit$coefficients
  })
  coefmat <- do.call("rbind", coeflist)
  coefmean <- colMeans(coefmat)
  coefsd <- apply(coefmat, 2, sd)
  coefCInorm <- cbind(coefmean - 1.96*coefsd, coefmean + 1.96*coefsd)
  coefmedian <- apply(coefmat, 2, median)
  coefCIperc <- apply(coefmat, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
  coefest <- cbind(coefmean, coefsd, coefCInorm, coefmedian, t(coefCIperc))
  colnames(coefest) <- c("Mean Est.", "Std. Error", "norm 2.5 %", "norm 97.5 %", "Median Est.", "perc 2.5 %", "perc 97.5 %")
  if(zero_infl){
    countcoefmat <- coefmat
    countcoefest <- coefest
    zerocoeflist <- lapply(gwqslist, function(i) i$fit$coefficients$zero)
    zerocoefmat <- do.call("rbind", zerocoeflist)
    zerocoefmean <- colMeans(zerocoefmat)
    zerocoefsd <- apply(zerocoefmat, 2, sd)
    zerocoefCInorm <- cbind(zerocoefmean - 1.96*zerocoefsd, zerocoefmean + 1.96*zerocoefsd)
    zerocoefmedian <- apply(zerocoefmat, 2, median)
    zerocoefCIperc <- apply(zerocoefmat, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
    zerocoefest <- cbind(zerocoefmean, zerocoefsd, zerocoefCInorm, zerocoefmedian, t(zerocoefCIperc))
    colnames(zerocoefest) <- c("Mean Est.", "Std. Error", "norm 2.5 %", "norm 97.5 %", "Median Est.", "perc 2.5 %", "perc 97.5 %")
    coefest <- list(countcoefest = countcoefest, zerocoefest = zerocoefest)
    coefmat <- list(countcoefmat = countcoefmat, zerocoefmat = zerocoefmat)
  }

    wl <- lapply(gwqslist, function(i) i$final_weights)
    if(dwqs){
      for (i in 1:length(wl)) {
        if(i==1){
          wmatpos <- wl[[1]][,1:2]
          wmatneg <- wl[[1]][,c(1,3)]
        }
        else{
          wmatpos <- suppressWarnings(merge(wmatpos, wl[[i]][,1:2], by = "mix_name"))
          wmatneg <- suppressWarnings(merge(wmatneg, wl[[i]][,c(1,3)], by = "mix_name"))
        }
      }
      nameswmat <- as.character(wmatpos[,1])
      wmatpos <- t(wmatpos[,-1])
      wmatneg <- t(wmatneg[,-1])
      colnames(wmatpos) <- colnames(wmatneg) <- nameswmat
      rownames(wmatpos) <- rownames(wmatneg) <- NULL
      wmeanpos <- colMeans(wmatpos)
      wmeanneg <- colMeans(wmatneg)
      wCIpos <- apply(wmatpos, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
      wCIneg <- apply(wmatneg, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
      final_weights_pos <- as.data.frame(cbind(wmeanpos, t(wCIpos)))
      final_weights_neg <- as.data.frame(cbind(wmeanneg, t(wCIneg)))
      names(final_weights_pos) <- c("Estimate pos", "2.5% pos", "97.5% pos")
      names(final_weights_neg) <- c("Estimate neg", "2.5% neg", "97.5% neg")
      final_weights <- cbind(final_weights_pos, final_weights_neg)
      wmat <- list(wmatpos = wmatpos, wmatneg = wmatneg)
    }
    else{
      for (i in 1:length(wl)) {
        if(i==1) wmat <- wl[[1]]
        else wmat <- suppressWarnings(merge(wmat, wl[[i]], by = "mix_name"))
      }
      nameswmat <- as.character(wmat[,1])
      wmat <- t(wmat[,-1])
      colnames(wmat) <- nameswmat
      rownames(wmat) <- NULL
      wmean <- colMeans(wmat)
      wCI <- apply(wmat, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
      final_weights <- as.data.frame(cbind(wmean, t(wCI)))
      names(final_weights) <- c("Estimate", "2.5%", "97.5%")
    }

  final_weights <- cbind(rownames(final_weights), final_weights)
  names(final_weights)[1] <- "mix_name"

  # estimate wqs index
  if(dwqs){
    dtf$pwqs <- pwqs <- as.numeric(Q%*%final_weights[match(mix_name, final_weights[,1]), 2])
    dtf$nwqs <- nwqs <- as.numeric(Q%*%final_weights[match(mix_name, final_weights[,1]), 5])
    dtf$wqs <- wqs <- NULL
  }
  else{
    dtf$wqs <- wqs <- as.numeric(Q%*%final_weights[match(mix_name, final_weights[,1]), 2])
    dtf$pwqs <- pwqs <- dtf$nwqs <- nwqs <- NULL
  }

  # scatterplot dataset
  if(all(grepl("wqs", attr(terms(formula), "term.labels")))) y_plot <- model.response(model.frame(formula, dtf), "any")
  else{
    formula_wo_wqs = remove_terms(formula, "wqs")
    if(family$family != "multinomial"){
      if(zero_infl){
        if(length(ff[[3]]) > 1 && identical(ff[[3]][[1]], as.name("|"))){
          if(all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))), "term.labels")))) f1 <- as.formula(paste0(ff[[2]], " ~ ", 1))
          else f1 <- remove_terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]]))), "wqs")
          if(all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[3]])))), "term.labels")))) f2 <- as.formula(paste0(ff[[2]], " ~ ", 1))
          else f2 <- remove_terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[3]]))), "wqs")
          formula_wo_wqs <- as.formula(paste0(deparse(f1), " | ", f2[[3]]))
        }
        fit <- zeroinfl(formula_wo_wqs, dtf, dist = family$family, link = zilink$name)
      }
      else{
        if(family$family == "negbin") fit = glm.nb(formula_wo_wqs, dtf)
        else fit = glm(formula_wo_wqs, dtf, family = family)
      }
      if(family$family == "binomial") y_plot = fit$y
      else y_plot = mean(fit$y) + resid(fit, type = "pearson")
    }
  }
    final_weights <- final_weights[order(final_weights[,2], decreasing = TRUE),]
    if(dwqs){
      y_adj_wqs_df <- as.data.frame(cbind(y_plot, pwqs, nwqs))
      names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", "y", "y_adj"), "pwqs", "nwqs")
    }
    else{
      y_adj_wqs_df <- as.data.frame(cbind(y_plot, wqs))
      names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", "y", "y_adj"), "wqs")
    }

  fit <- list(coefficients = NULL, aic = NULL, dispersion = NULL, deviance = NULL, df.residual = NULL,
              null.deviance = NULL, df.null = NULL, theta = NULL, SE.theta = NULL)

  fit$coefficients <- coefest
  if(zero_infl) fit$aic <- mean(2*(nrow(coefest[[1]]) + nrow(coefest[[2]]) + ifelse(family$family == "negbin", 1, 0)) - 2*sapply(gwqslist, function(i) i$fit$loglik), na.rm = TRUE)
  else{
    fit$dispersion <- mean(sapply(gwqslist, function(i) summary(i$fit)$dispersion), na.rm = TRUE)
    fit$deviance <- mean(sapply(gwqslist, function(i) i$fit$deviance), na.rm = TRUE)
    fit$df.residual <- gwqslist[[1]]$fit$df.residual
    fit$null.deviance <- mean(sapply(gwqslist, function(i) i$fit$null.deviance), na.rm = TRUE)
    fit$df.null <- gwqslist[[1]]$fit$df.null
    fit$aic <- mean(sapply(gwqslist, function(i) i$fit$aic), na.rm = TRUE)
  }
  if(family$family == "negbin"){
    fit$theta <- mean(sapply(gwqslist, function(i) i$fit$theta), na.rm = TRUE)
    fit$SE.theta <- sd(sapply(gwqslist, function(i) i$fit$theta), na.rm = TRUE)
  }

  results <- list(fit = fit, final_weights = final_weights, wqs = wqs, qi = qi, y_wqs_df = y_adj_wqs_df,
                  family = family, call = cl, formula = formula, mix_name = mix_name, stratified = stratified,
                  q = q, zero_infl = zero_infl, zilink = zilink,
                  dwqs = dwqs, data = dtf, gwqslist = gwqslist, coefmat = coefmat, wmat = wmat, rh = rh)

  if(zero_infl) results$formula <- ff
  class(results) <- "gwqs"

  return(results)
}



#' @export gwqs_multinom
#' @rdname main_functions
#'
#' @usage  gwqs_multinom(formula, data, na.action, weights, mix_name, stratified, rh = 1, b = 100,
#'                       b1_pos = c(TRUE, TRUE), b_constr = FALSE, q = 4,
#'                       validation = 0.6, validation_rows = NULL,
#'                       signal = c("t2", "t3", "one", "abst", "expt"),
#'                       rs = FALSE, n_vars = NULL,
#'                       zilink = c("logit", "probit", "cloglog", "cauchit", "log"),
#'                       seed = NULL, wp = NULL, wn = NULL, plan_strategy = "sequential",
#'                       lambda = 0, optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
#'                       control = list(trace = FALSE, maxit = 2000, reltol = 1e-9),
#'                       b1_constr = NULL, ...)

gwqs_multinom <- function(formula, data, na.action, weights, mix_name, stratified, rh = 1, b = 100,
                          b1_pos = c(TRUE, TRUE), b_constr = FALSE, q = 4,
                          validation = 0.6, validation_rows = NULL,
                          signal = c("t2", "t3", "one", "abst", "expt"),
                          rs = FALSE, n_vars = NULL,
                          zilink = c("logit", "probit", "cloglog", "cauchit", "log"),
                          seed = NULL, wp = NULL, wn = NULL, plan_strategy = "sequential",
                          lambda = 0, optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
                          control = list(trace = FALSE, maxit = 2000, reltol = 1e-9),
                          b1_constr = NULL, ...){

  wqsMultinom <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars,
                          nwqsvars, dwqs, formula, bwghts, stratified, stratlev, b1_pos, b_constr, lambda,
                          intpos, wqsint, intvars){

    if(b_constr){
      par_pos <- which(grepl("wqs", names(initp)))
      initp[par_pos] <- sapply(1:length(b1_pos), function(i) ifelse(b1_pos[i], abs(initp[par_pos[i]]), -abs(initp[par_pos[i]])))
    }
    w <- matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels-1)
    pen <- sum(c(w))
    w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
    bXm[, wqsvars] <- bQ%*%w
    b_covs = matrix(0, kx, n_levels-1)
    i = 1:(kx/(n_levels-1))
    for (j in 0:(n_levels-2)){
      b_covs[(kx/(n_levels-1))*j+i, j+1] = initp[(kx/(n_levels-1))*j+i]
    }
    term = bXm%*%b_covs + boffset
    f = -sum((diag(bY%*%t(term)) - log(1 + rowSums(exp(term))))*bwghts) + lambda*pen

    return(f)
  }

  wqsMultinom_dwqs <- function(initp, kw, bXm, bY, boffset, bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars,
                          nwqsvars, dwqs, formula, bwghts, stratified, stratlev, b1_pos, b_constr, lambda,
                          intpos, wqsint, intvars){

    pwqs_site <- which(grepl("pwqs", names(initp)))
    nwqs_site <- which(grepl("nwqs", names(initp)))
    initp[pwqs_site] <- abs(initp[pwqs_site])
    initp[nwqs_site] <- -abs(initp[nwqs_site])
    # wp <- matrix(abs(initp[(kx + 1):(kx + kw)]), kw, n_levels-1)
    wp <- matrix(initp[(kx + 1):(kx + kw)]^2, kw, n_levels-1)
    # wn <- matrix(abs(initp[(kx + kw + 1):(kx + 2*kw)]), kw, n_levels-1)
    wn <- matrix(initp[(kx + kw + 1):(kx + 2*kw)]^2, kw, n_levels-1)
    pen <- sum(c(wp, wn))
    wp <- apply(wp, MARGIN = 2, FUN = function(i) i/sum(i))
    wn <- apply(wn, MARGIN = 2, FUN = function(i) i/sum(i))
    bXm[, pwqsvars] <- bQ%*%wp
    bXm[, nwqsvars] <- bQ%*%wn
    b_covs = matrix(0, kx, n_levels-1)
    i = 1:(kx/(n_levels-1))
    for (j in 0:(n_levels-2)){
      b_covs[(kx/(n_levels-1))*j+i, j+1] = initp[(kx/(n_levels-1))*j+i]
    }
    term = bXm%*%b_covs + boffset
    f = -sum((diag(bY%*%t(term)) - log(1 + rowSums(exp(term))))*bwghts) + lambda*pen

    return(f)
  }

  if(!is.null(b1_constr)){
    b_constr <- b1_constr
    .Deprecated("The argument 'b1_constr' is deprecated, use 'b_constr' instead.")
  }

  family <- list(family = "multinomial")

  one <- function(x) rep(1, length(x))

  t2 <- function(x) x^2

  t3 <- function(x) x^3

  expt <- function(x) exp(abs(x))-1

  wqs_in_formula <- any(grepl("wqs", rownames(attr(terms(formula), "factors"))))
  pwqs_in_formula <- any(grepl("pwqs", rownames(attr(terms(formula), "factors"))))
  nwqs_in_formula <- any(grepl("nwqs", rownames(attr(terms(formula), "factors"))))
  if(!wqs_in_formula & (!pwqs_in_formula | !nwqs_in_formula))
    stop("'formula' must contain either 'wqs' term or 'pwqs' and 'nwqs' terms: e.g. y ~ wqs + ...\n y ~ pwqs + nwqs + ...\n")

  dwqs <- pwqs_in_formula & nwqs_in_formula

  if(dwqs) stop("The 2iWQS is not supported yet by the 'gwqs_multinom' function\n")

  if(dwqs) objfn <- wqsMultinom_dwqs
  else objfn <- wqsMultinom
  optim.method <- match.arg(optim.method)

  signal <- match.arg(signal)
  signal <- switch(signal,
                   "one" = one,
                   "abst" = abs,
                   "t2" = t2,
                   "t3" = t3,
                   "expt" = expt)

  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", "weights", "stratified"), names(mc), 0)
  dtf <- mc[c(1, m)]
  dtf[[2]] <- data
  dtf[[1]] <- as.name("selectdatavars")
  dtf <- eval(dtf, parent.frame())

  l <- list(...)
  solve_dir_issue <- ifelse(is.null(l$solve_dir_issue), FALSE, l$solve_dir_issue)

  # defining quantile variables
  if (is.null(q)) {Q = as.matrix(dtf[, mix_name]); qi = NULL}
  else {
    q_f = gwqs_rank(dtf, mix_name, q)
    Q = q_f$Q
    qi = q_f$qi
  }

  # create stratified variables if stratified is not NULL
  m <- match(c("stratified"), names(mc), 0)
  if(!is.null(l$stratified_rh)) m <- 0
  if (m){
    strtfd_out <- stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name <- strtfd_out$mix_name
  }
  else stratified <- NULL
  if(!is.null(l$stratified_rh)) stratified <- l$stratified_rh
  stratlev <- ifelse(is.null(stratified), 1, nlevels(unlist(dtf[, stratified])))

  n_levels <- nlevels(eval(formula[[2]], envir = dtf))
  if(n_levels == 0) stop("y must be of class factor when 'family = \"multinomial\"'\n")

  if(!is.null(seed) & !is.numeric(seed)) stop("seed must be numeric or NULL\n")
  if(!is.null(seed)) set.seed(seed)
  fseed <- TRUE

  N <- nrow(dtf)

  # splitting the dataset
  if(!is.numeric(validation) | validation < 0 | validation >= 1) stop("'validation' must be numeric >= 0 and < 1\n")
  if(is.null(validation_rows)) validation_rows <- lapply(1:rh, function(i) sample(c(F, T), N, replace = T, prob = c(1-validation, validation)))
  else if(length(validation_rows) != rh | !any(unique(unlist(validation_rows)) %in% c(F, T))) stop("'validation_rows' values must be 0 or 1\n")

  # estimate starting weights
  if(dwqs){
    mc_pwqs <- mc_nwqs <- mc
    mc_pwqs$data <- mc_nwqs$data <- data
    mc_pwqs$mix_name <- mc_nwqs$mix_name <- mix_name
    mc_pwqs$validation <- mc_nwqs$validation <- 0
    mc_pwqs$plan_strategy <- mc_nwqs$plan_strategy <- "sequential"
    mc_pwqs$b_constr <- mc_nwqs$b_constr <- TRUE
    mc_pwqs$b1_pos <- rep(TRUE, n_levels-1)
    mc_nwqs$b1_pos <- rep(FALSE, n_levels-1)
    mc_pwqs$validation_rows <- mc_nwqs$validation_rows <- NULL
    if(rs) mc_pwqs$b <- mc_nwqs$b <- 100
    else mc_pwqs$b <- mc_nwqs$b <- 1
    mc_pwqs$formula <- mc_nwqs$formula <- as.formula(gsub("pwqs", "wqs", deparse(formula)))
    mc_pwqs$formula <- mc_nwqs$formula <- remove_terms(mc_pwqs$formula, "nwqs")
    mc_pwqs$solve_dir_issue <- mc_nwqs$solve_dir_issue <- "inverse"

    if(is.null(wp)){
      pwqs0 <- tryCatch(eval(mc_pwqs), error = function(e) NULL)
      if(is.null(pwqs0)) wp <- cbind(rep(1/length(mix_name), length(mix_name)), rep(1/length(mix_name), length(mix_name)))
      else wp <- as.matrix(pwqs0$final_weights[match(mix_name, pwqs0$final_weights$mix_name), -1])
    }
    if(is.null(wn)){
      nwqs0 <- tryCatch(eval(mc_nwqs), error = function(e) NULL)
      if(is.null(nwqs0)) wn <- cbind(rep(1/length(mix_name), length(mix_name)), rep(1/length(mix_name), length(mix_name)))
      else wn <- as.matrix(nwqs0$final_weights[match(mix_name, nwqs0$final_weights$mix_name), -1])
    }
  }

  if(is.null(n_vars)) n_vars = round(sqrt(length(mix_name)))

  # parameters estimation and model fitting
  m <- match(c("weights"), names(mc), 0)
  if(m[1]) dtf$wghts <- wghts <- unlist(dtf[, weights, drop = FALSE])
  else  dtf$wghts <- wghts <- rep(1, N)

  ts = "z"

  if(!is.numeric(b)) stop("'b' must be a number\n")

  Xnames <- parnames(dtf, formula, NULL)
  kx <- length(Xnames)
  level_names <- levels(eval(formula[[2]], envir = dtf))
  Xnames <- c(sapply(level_names[-1], function(i) paste0(Xnames[1:kx], "_", i, "_vs_", level_names[1])))
  kx <- kx*(n_levels-1)
  if(dwqs){
    pwqsvars = Xnames[grepl("^pwqs_", Xnames)]
    nwqsvars = Xnames[grepl("^nwqs_", Xnames)]
    dtf[, c(pwqsvars, nwqsvars)] <- 0
    dtf[, pwqsvars] <- Q%*%wp
    dtf[, nwqsvars] <- Q%*%wn
    wp <- as.vector(wp)
    wn <- as.vector(wn)
    wqsvars <- NULL
  }
  else{
    wqsvars = Xnames[grepl("^wqs_", Xnames)]
    dtf[, wqsvars] <- 0
    pwqsvars <- nwqsvars <- NULL
  }
  kw <- ncol(Q)
  initp <- values.0(kw, Xnames, kx, formula, ff=NULL, wghts, dtf, stratified, stratlev, b1_pos, family, wp, wn, dwqs, zilink, zero_infl=FALSE)
  mf <- model.frame(formula, dtf)
  Y <- model.response(mf, "any")
  Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- rep(0, nrow(dtf))

  if(dwqs){
    Xl = lapply(1:length(pwqsvars), function(i){
      fmp <- gsub("pwqs", pwqsvars[i], format(formula))
      fmnp <- gsub("nwqs", nwqsvars[i], fmp)
      fmi <- as.formula(fmnp)
      model.matrix(fmi, data = dtf)
    })
  }
  else{
    Xl = lapply(wqsvars, function(i){
      fmi <- as.formula(gsub("wqs", i, format(formula)))
      model.matrix(fmi, data = dtf)
    })
  }
  Xm = do.call("cbind", Xl)


  intpos <- grepl("wqs", colnames(Xm)) & grepl(":", colnames(Xm))
  intvarname <- colnames(Xm)[intpos]
  wqsint <- any(intpos)
  intvars <- gsub(":wqs", "", gsub("wqs:", "", colnames(Xm)[intpos]))

  if(wqsint) stop("The 'gwqs_multinom' function does not support yet the interaction between a covariate and the WQS index\n")

  plan_strategy_rh <- ifelse(rh==1, "sequential", plan_strategy)
  plan_strategy_b <- ifelse(rh==1, plan_strategy, "sequential")

  plan(plan_strategy_rh)
  gwqslist <- future_lapply(X = 1:rh, FUN = function(i){
  # gwqslist <- lapply(1:rh, function(i){
    if(control$trace) cat("start opt\n")
    plan(plan_strategy_b)
    param <- future_lapply(X = 1:b, FUN = optim.f_multinom, objfn = objfn, Y = Y[!validation_rows[[i]],],
    # param <- lapply(X = 1:b, FUN = optim.f_multinom, objfn = objfn, Y = Y[!validation_rows[[i]],],
                    Xm = Xm[!validation_rows[[i]],], Q = Q[!validation_rows[[i]],], offset = offset[!validation_rows[[i]]], wghts = wghts[!validation_rows[[i]]],
                    initp = initp, n_levels = n_levels, level_names = level_names, wqsvars = wqsvars,
                    pwqsvars = pwqsvars, nwqsvars = nwqsvars,
                    b1_pos = b1_pos, b_constr = b_constr, n_vars = n_vars, dwqs = dwqs, rs = rs,
                    formula = formula, kx = kx, kw = kw,
                    Xnames = Xnames, stratified = stratified, b = b, optim.method = optim.method, control = control,
                    lambda = lambda, stratlev = stratlev, intpos = intpos, wqsint = wqsint, intvars = intvars, future.seed = fseed)

    conv <- c(sapply(param, function(i) i$conv))
    counts <- c(sapply(param, function(i) i$counts))
    val <- c(sapply(param, function(i) i$val))
    mex <- lapply(param, function(i) i$mex)
    bindex <- lapply(param, function(i) i$bindex)
    slctd_vars <- lapply(param, function(i) i$slctd_vars)

    if(rs){
      plan(plan_strategy)
      param <- future_lapply(X = 1:b, FUN = set_par_names, slctd_vars, param, q_name = colnames(Q), family = family,
                             dwqs = dwqs, future.seed = FALSE)
    }

    n_levels <- dim(param[[1]]$par_opt)[2]+1

    ### NEED TO CREATE TOLERANCE FOR MULTINOMIAL REGRESSION FOR DWQS AND CHECK FOR ZERO_INFL

    build_bres <- function(wqs_var_name, interval){
      wqs_site <- which(grepl(paste0("^", wqs_var_name, "_"), rownames(param[[1]]$mfit$m_f$coefficients)))
      wght_matrix <- lapply(1:(n_levels-1), function(j) do.call("rbind", lapply(param, function(i) i$par_opt[interval,j])))
      b1 <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$nlm_out$estimate[j]))
      se <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$Standard_Error[j]))
      stat <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$stat[j]))
      p_val <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$p_value[j]))
      if(dwqs){
        if(wqs_var_name == "pwqs") formula_dwqs <- update(formula, pwqs ~ . - pwqs)
        else if(wqs_var_name == "nwqs") formula_dwqs <- update(formula, nwqs ~ . - nwqs)
        tol <- lapply(1:(n_levels-1), function(j) sapply(param, function(i){
          lmwqs <- lm(formula_dwqs, cbind(pwqs = i$mfit$pwqs[,j], nwqs = i$mfit$nwqs[,j], data[i$bindex,]))
          (1/(1-summary(lmwqs)$r.squared))^(-1)
        }))
      }

      n_non_conv = sum(conv == 1)
      if(n_non_conv == 0 & control$trace) cat(paste0("The optimization function always converged\n"))
      else if(n_non_conv == b & rh == 1) stop("The optimization function never converged\n")
      else if(n_non_conv == b & rh > 1) return(NULL)
      else if(control$trace) cat(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))
      if(control$trace) cat(paste0("There are ", ifelse(b1_pos, sum(b1 >= 0, na.rm = T), sum(b1 <= 0, na.rm = T)),
                                   ifelse(b1_pos, " positive", " negative"), " bootstrapped b1 out of ", b, "\n"))

      # estimate mean weight for each component (exclude weights from iterations with failed convergence)
      if(dwqs) bres <- Map(cbind, wght_matrix, b1, se, stat, p_val, tol)
      else bres <- Map(cbind, wght_matrix, b1, se, stat, p_val)
      bres <- lapply(bres, as.data.frame)
      if(dwqs) bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", "stat", "p_val", "tol"))
      else bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", "stat", "p_val"))
      strata_names <- gsub(paste(wqs_var_name, "_"), "", rownames(param[[1]]$mfit$m_f$coefficients)[wqs_site])
      names(bres) <- strata_names

      return(bres)
    }

    if(dwqs){
      bresp <- build_bres("pwqs", 1:ncol(Q))
      bresn <- build_bres("nwqs", (ncol(Q)+1):(2*ncol(Q)))
      bres <- list(bresp = bresp, bresn = bresn)
    }
    else{
      bresp <- bresn <- NULL
      bres <- build_bres("wqs", 1:ncol(Q))
    }
    if(dwqs) strata_names <- c(names(bresp), names(bresn))
    else strata_names <- names(bres)

    if(is.null(bres)) return(NULL)

    mean_weight_f <- function(bres, b1_pos, par){
      mean_weight <- lapply(1:(n_levels-1), function(i){
        if(rs) bres[[i]][mix_name][is.na(bres[[i]][mix_name])] <- 0
        if(b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 > 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[[i]][bres[[i]]$b1 > 0 & conv == 0, par]))
        else if(!b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 < 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[[i]][bres[[i]]$b1 < 0 & conv == 0, par]))
        if (all(is.nan(w_t))){
          if(solve_dir_issue == "inverse") w_t <- 1/apply(bres[[i]][, mix_name], 2, weighted.mean, signal(bres[[i]][, par]))/sum(1/apply(bres[[i]][, mix_name], 2, weighted.mean, signal(bres[[i]][, par])))
          else if(solve_dir_issue == "average") w_t <- rep(1/length(mix_name), length(mix_name))
          else if(!(solve_dir_issue %in% c("average", "inverse")) & rh == 1) stop(paste0("There are no ", ifelse(b1_pos[i], "positive", "negative"), " b1 in the bootstrapped models for ", strata_names[i], "\n"))
          else if(!(solve_dir_issue %in% c("average", "inverse")) & rh > 1) return(NULL)
        }
        return(w_t)
      })
      mean_weight <- list.cbind(mean_weight)

      return(mean_weight)
    }

    if(dwqs){
      mean_weight_p <- mean_weight_f(bresp, b1_pos = rep(TRUE, length(bresp)), par = "tol")
      if(inherits(mean_weight_p, "gwqs")) return(mean_weight_p)
      mean_weight_n <- mean_weight_f(bresn, b1_pos = rep(FALSE, length(bresp)), par = "tol")
      if(inherits(mean_weight_n, "gwqs")) return(mean_weight_n)
      mean_weight <- rbind(mean_weight_p, mean_weight_n)
    }
    else mean_weight <- mean_weight_f(bres, b1_pos = b1_pos, par = "stat")

    if(is.null(mean_weight)) return(NULL)

    # fit the final model with the estimated weights
    if(validation==0) validation_rows <- lapply(validation_rows, function(i) !i)
    wqs_model = model.fit_multinom(mean_weight, dtf[validation_rows[[i]],], Q[validation_rows[[i]],], Y[validation_rows[[i]],],
                          dwqs, formula, wghts[validation_rows[[i]]], offset[validation_rows[[i]]], initp, Xnames, n_levels,
                          level_names, wqsvars, pwqsvars, nwqsvars, stratified, b1_pos, kx, kw, intpos, wqsint, intvars)

    if(all(attr(terms(formula), "term.labels") %in% "wqs") | all(attr(terms(formula), "term.labels") %in% c("pwqs", "nwqs"))) y_plot <- model.response(model.frame(formula, dtf[validation_rows[[i]],]), "any")
    else{
      if(dwqs){
        formula_wo_wqs <- remove_terms(formula, "pwqs")
        formula_wo_wqs <- remove_terms(formula_wo_wqs, "nwqs")
      }
      else formula_wo_wqs <- remove_terms(formula, "wqs")
    }

    # Plots
    if(dwqs){
      data_plot <- data.frame(mix_name, mean_weight_p, mean_weight_n, stringsAsFactors = TRUE)
      pwqs_index <- wqs_model$pwqs
      nwqs_index <- wqs_model$nwqs
      wqs_index <- NULL
      dtf[, c(colnames(pwqs_index), colnames(nwqs_index))] <- cbind(Q%*%mean_weight_p, Q%*%mean_weight_n)
    }
    else{
      data_plot <- data.frame(mix_name, mean_weight, stringsAsFactors = TRUE)
      pwqs_index <- nwqs_index <- NULL
      wqs_index <- wqs_model$wqs
      dtf[, colnames(wqs_index)] <- Q%*%mean_weight
    }

    names(data_plot)[2:ncol(data_plot)] = strata_names
    data_plot <- data_plot[order(data_plot[, strata_names[1]], decreasing = TRUE),]

    if(rh==1){
      Y <- model.response(model.frame(formula, dtf[validation_rows[[i]],]), "any")
      level_names <- levels(Y)
      Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
      colnames(Y) <- gsub("pwqs_", "", strata_names[1:(n_levels-1)])
      if(dwqs){
        y_adj_wqs_df = do.call("rbind", lapply(1:length(colnames(Y)), function(i){
          data.frame(level = colnames(Y)[i],
                     y = Y[rowSums(Y[, -which(colnames(Y)==colnames(Y)[i]), drop = F]) == 0, colnames(Y)[i]],
                     pwqs = wqs_model$pwqs[rowSums(Y[, -which(colnames(wqs_model$pwqs)==colnames(wqs_model$pwqs)[i]), drop = F]) == 0, colnames(wqs_model$pwqs)[i]],
                     nwqs = wqs_model$nwqs[rowSums(Y[, -which(colnames(wqs_model$nwqs)==colnames(wqs_model$nwqs)[i]), drop = F]) == 0, colnames(wqs_model$nwqs)[i]],
                     stringsAsFactors = TRUE)
        }))
      }
      else{
        y_adj_wqs_df <- do.call("rbind", lapply(strata_names, function(i){
          data.frame(level = i,
                     y = Y[rowSums(Y[, -which(colnames(Y)==i), drop = F]) == 0, i],
                     wqs = wqs_model$wqs[rowSums(Y[, -which(colnames(wqs_model$wqs)==i), drop = F]) == 0, i],
                     stringsAsFactors = TRUE)
        }))
      }
    }

    dtf$wqs <- dtf$pwqs <- dtf$nwqs <- NULL

    # creating the list of elements to return
    results = list(fit = wqs_model$m_f, final_weights = data_plot, conv = conv, bres = bres, wqs = wqs_index,
                   pwqs = pwqs_index, nwqs = nwqs_index, qi = qi,
                   bindex = bindex, slctd_vars = slctd_vars, validation_rows = validation_rows[[i]], vindex = validation_rows[[i]],
                   y_wqs_df = if(rh==1) y_adj_wqs_df else NULL, family = family, call = cl, formula = formula, mix_name = mix_name,
                   stratified = stratified, q = q, n_levels = n_levels,
                   levelnames = strata_names, dwqs = dwqs, data = dtf, objfn_values = val, optim_messages = mex, rh = rh)
    class(results) <- "gwqs"

    return(results)
  }, future.seed = fseed)

  if(rh == 1) return(gwqslist[[1]])

  gwqslist <- gwqslist[!vapply(gwqslist, is.null, logical(1))]

  if(length(gwqslist) == 0) stop("The model never converged. Try to increase the number of repeated holdout 'rh' or to change the direction of the association.\n")
  if(length(gwqslist) == 1) stop("The model converged for only a single iteration. Try to increase the number of repeated holdout 'rh'\n")

  if(control$trace) cat(paste0("The model converged ", length(gwqslist), " times out of ", rh, " repeated holdout iterations.\n"))

  coeflist <- lapply(gwqslist, function(i){
    tmp <- i$fit$coefficients[,1]
    names(tmp) <- rownames(i$fit$coefficients)
    tmp
  })
  coefmat <- do.call("rbind", coeflist)
  coefmean <- colMeans(coefmat)
  coefsd <- apply(coefmat, 2, sd)
  coefCInorm <- cbind(coefmean - 1.96*coefsd, coefmean + 1.96*coefsd)
  coefmedian <- apply(coefmat, 2, median)
  coefCIperc <- apply(coefmat, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
  coefest <- cbind(coefmean, coefsd, coefCInorm, coefmedian, t(coefCIperc))
  colnames(coefest) <- c("Mean Est.", "Std. Error", "norm 2.5 %", "norm 97.5 %", "Median Est.", "perc 2.5 %", "perc 97.5 %")

  n_levels <- nlevels(eval(formula[[2]], envir = data))
  levelnames <- levels(eval(formula[[2]], envir = data))
  wll <- lapply(2:n_levels, function(j){
    wl <- lapply(gwqslist, function(i) i$final_weights[,c(1,j)])
    for (i in 1:length(wl)) {
      if(i==1) wmat <- wl[[1]]
      else wmat <- suppressWarnings(merge(wmat, wl[[i]], by = "mix_name"))
    }
    nameswmat <- as.character(wmat[,1])
    names(wmat) <- NULL
    wmat <- t(wmat[,-1])
    colnames(wmat) <- nameswmat
    wmean <- colMeans(wmat)
    wCI <- apply(wmat, 2, function(i) quantile(i, probs = c(0.025, 0.975)))
    final_weights <- cbind(wmean, t(wCI))
    colnames(final_weights) <- paste0(c("", "2.5 % ", "97.5% "), levelnames[j], "_vs_", levelnames[1])
    list(wmat, final_weights)
  })
  wmat <- lapply(wll, function(i) i[[1]])
  names(wmat) <- levelnames <- paste0(levelnames[-1], "_vs_", levelnames[1])
  final_weights <- as.data.frame(do.call("cbind", lapply(wll, function(i) i[[2]])))
  final_weights <- cbind(rownames(final_weights), final_weights)
  names(final_weights)[1] <- "mix_name"

  # estimate wqs index
  wqs <- Q%*%as.matrix(final_weights[match(mix_name, final_weights[,1]), levelnames])
  dtf <- cbind(dtf, wqs)

  # scatterplot dataset
  if(all(grepl("wqs", attr(terms(formula), "term.labels")))) y_plot <- model.response(model.frame(formula, dtf), "any")
  else formula_wo_wqs = remove_terms(formula, "wqs")
  final_weights <- final_weights[order(final_weights[, levelnames[1]], decreasing = TRUE),]
  Y <- model.response(model.frame(formula, dtf), "any")
  level_names <- levels(Y)
  Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
  colnames(Y) <- levelnames
  y_adj_wqs_df <- do.call("rbind", lapply(levelnames, function(i){
    if(n_levels > 3) data.frame(level = i,
                                y = Y[rowSums(Y[, -which(colnames(Y)==i)]) == 0, i],
                                wqs = wqs[rowSums(Y[, -which(colnames(wqs)==i)]) == 0, i],
                                stringsAsFactors = TRUE)
    else data.frame(level = i,
                    y = Y[Y[, -which(colnames(Y)==i)] == 0, i],
                    wqs = wqs[Y[, -which(colnames(wqs)==i)] == 0, i],
                    stringsAsFactors = TRUE)
  }))

  dtf$wqs <- NULL

  fit <- list(coefficients = NULL, aic = NULL, dispersion = NULL, deviance = NULL, df.residual = NULL,
              null.deviance = NULL, df.null = NULL, theta = NULL, SE.theta = NULL)

  fit$coefficients <- coefest
  fit$aic <- mean(sapply(gwqslist, function(i) 2*dim(i$fit$coefficients)[1] + 2*(i$fit$nlm_out$minimum)), na.rm = TRUE)
  fit$deviance <- mean(sapply(gwqslist, function(i) 2*i$fit$nlm_out$minimum), na.rm = TRUE)

  results <- list(fit = fit, final_weights = final_weights, wqs = wqs, qi = qi, y_wqs_df = y_adj_wqs_df,
                  call = cl, formula = formula, mix_name = mix_name, stratified = stratified,
                  q = q, n_levels = n_levels, levelnames = levelnames, family = family, zero_infl = FALSE,
                  dwqs = dwqs, data = dtf, gwqslist = gwqslist, coefmat = coefmat, wmat = wmat, rh = rh)

  class(results) <- "gwqs"

  return(results)
}



#' @export gwqsrh
#' @rdname main_functions
#'
#' @usage  gwqsrh(formula, data, na.action, weights, mix_name, stratified, rh = 1, b = 100,
#'                b1_pos = TRUE, bint_cont_pos = NULL, bint_cat_pos = NULL, b_constr = FALSE,
#'                zero_infl = FALSE, q = 4, validation = 0.6, validation_rows = NULL,
#'                family = gaussian, signal = c("t2", "t3", "one", "abst", "expt"),
#'                rs = FALSE, n_vars = NULL,
#'                zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
#'                wp = NULL, wn = NULL, plan_strategy = "sequential", lambda = 0,
#'                optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
#'                control = list(trace = FALSE, maxit = 2000, reltol = 1e-9), ...)

gwqsrh <- function(formula, data, na.action, weights, mix_name, stratified, rh = 1, b = 100,
                   b1_pos = TRUE, bint_cont_pos = NULL, bint_cat_pos = NULL, b_constr = FALSE,
                   zero_infl = FALSE, q = 4, validation = 0.6, validation_rows = NULL,
                   family = gaussian, signal = c("t2", "t3", "one", "abst", "expt"),
                   rs = FALSE, n_vars = NULL,
                   zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
                   wp = NULL, wn = NULL, plan_strategy = "sequential", lambda = 0,
                   optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
                   control = list(trace = FALSE, maxit = 2000, reltol = 1e-9), ...){

  .Defunct("gwqs", package = "gWQS", "The 'gwqsrh' function is deprecated, please specify the parameter 'rh' within the 'gwqs' function.")

}

#' Methods for gwqs objects
#'
#' Methods for extracting information from fitted Weighted Quantile Sum (WQS) regression model objects
#' of class "gwqs".
#'
#' @param object,x An object of class "gwqs" as returned by \link[gWQS]{gwqs}.
#' @param sumtype Type of summary statistic to be used: "norm" takes the mean of the estimated parameters on the
#' validation sets and the 95% CI assuming a normal distribution of the parameters, while "perc" uses the median
#' as the parameters estimates and the 2.5, 97.5 percentiles as CI. This option is only available for objects of
#' class \code{gwqsrh}.
#' @param digits The number of significant digits to use when printing.
#' @param newdata Optionally, a data frame in which to look for variables with which to predict.
#' If omitted, the original observations are used.
#' @param type Character specifying the type of predictions, fitted values or residuals, respectively.
#' For details see below.
#' @param model Character specifying for which component of the model the varance-covariance matrix
#' should be extracted when zero_infl = TRUE.
#' @param ... Further arguments to be passed.
#'
#' @details
#' A set of standard extractor functions for fitted model objects is available for objects of class "gwqs",
#' including methods to the generic functions print and summary which print the estimated coefficients
#' along with some further information. As usual, the summary method returns an object of class
#' "summary.gwqs" containing the relevant summary statistics which can subsequently be printed using the
#' associated print method.
#'
#' The methods for \link[stats]{coef} and \link[stats]{vcov} by default return a single vector of
#' coefficients (a matrix when family = "multinomial") and their associated covariance matrix, respectively.
#' By setting the model argument, the estimates for the corresponding model components can be extracted.
#'
#' Both the \link[stats]{fitted} and \link[stats]{predict} methods can compute fitted responses. The latter
#' sets the default on the scale of the linear predictors; the alternative "response" is on the scale of the
#' response variable. Thus for a default binomial model the default predictions are of log-odds
#' (probabilities on logit scale) and type = "response" gives the predicted probabilities. Type can be equal
#' to "prob", "count" or "zero" when zero_infl = T to estimate the predicted density (i.e., probabilities for
#' the observed counts), the predicted mean from the count component (without zero inflation) and the
#' predicted probability for the zero component. Type = "class" allow to predict the dependent variable
#' categories when family = "multinomial". The "terms" option returns a matrix giving the fitted values of
#' each term in the model formula on the linear predictor scale.
#'
#' The \link[stats]{residuals} method allows to extracts model residuals from the objects of class "gwqs".
#'
#' @return All these methods return the classic output as for the corresponding glm, glm.nb, multinom and
#' zeroinfl classes. Only the predict method gives a different output made of the following values.
#'
#' \item{df_pred}{A data.frame containing the dependent varible and the predicted values.}
#' \item{Q}{The matrix containing the new dataset quantiled variables of the elements included in the mixture.}
#' \item{qi}{A list of vectors containing the cut points used to determine the quantiled variables.}
#' \item{wqs}{The vetor containing the wqs index built on the new dataset.}
#'
#' @author
#' Stefano Renzetti, Paul Curtin, Allan C Just, Ghalib Bello, Chris Gennings
#'
#' @examples
#' toxic_chems = names(wqs_data)[1:34]
#' set.seed(1234)
#' rws <- sample(1:500, 150)
#' results = gwqs(yLBX ~ wqs, mix_name = toxic_chems, data = wqs_data[-rws,], q = 4, validation = 0.6,
#'                b = 2, b1_pos = TRUE, b_constr = FALSE, family = gaussian)
#'
#' # to test the significance of the covariates
#' summary(results)
#'
#' # extract regression coefficients
#' coef(results)
#'
#' # estimate variance-covariance matrix
#' vcov(results)
#'
#' # estimate fitted values
#' fitted(results)
#'
#' # estimate regression residuals
#' residuals(results)
#'
#' # estimate predicted values on the left part of wqs_data
#' pred_res <- predict(results, wqs_data[rws,])
#' pred_res$df_pred
#'
#' @importFrom utils tail
#'
#' @rawNamespace S3method(summary, gwqs)
#' @rdname methods

summary.gwqs <- function(object, sumtype = c("norm", "perc"), ...){
  if(object$rh == 1){
    if(object$family$family == "multinomial"){
      ans <- vector("list", 7)
      ans$call <- object$call
      ans$is.binomial <- FALSE
      ans$digits <- options()$digits
      ans$coefficients <- matrix(object$fit$coefficients$Estimate, nrow = dim(object$fit$coefficients)[1]/2, byrow = F)
      ans$standard.errors <- matrix(object$fit$coefficients$Standard_Error, nrow = dim(object$fit$coefficients)[1]/2, byrow = F)
      colnames(ans$coefficients) <- colnames(ans$standard.errors) <- levels(unlist(object$data[, all.vars(object$formula)[1]]))[-1]
      if(object$dwqs) object$data$pwqs <-object$data$nwqs <- 0
      else object$data$wqs <- 0
      rownames(ans$coefficients) <- rownames(ans$standard.errors) <- colnames(model.matrix(object$formula, object$data))
      ans$AIC <- 2*dim(object$fit$coefficients)[1] + 2*(object$fit$nlm_out$minimum)
      ans$deviance <- 2*object$fit$nlm_out$minimum
    }
    else{
      object$fit$call <- object$call
      if(object$zero_infl)
        ans <- summary(object$fit)
      else if(object$family$family == "negbin") ans <- summary(object$fit)
      else ans <- summary(object$fit)
    }
    ans$zero_infl <- object$zero_infl
    ans$family <- object$family
    ans$rh <- object$rh
    class(ans) <- "summary.gwqs"
    return(ans)
  }

  else{
    sumtype <- match.arg(sumtype)
    if(object$zero_infl){
      if(sumtype == "norm"){
        countcoefest <- object$fit$coefficients$countcoefest[,c(1,2,3,4)]
        colnames(countcoefest) <- c("Estimate", "Std. Error", "2.5 %", "97.5 %")
        zerocoefest <- object$fit$coefficients$zerocoefest[,c(1,2,3,4)]
        colnames(zerocoefest) <- c("Estimate", "Std. Error", "2.5 %", "97.5 %")
      }
      else if(sumtype == "perc"){
        countcoefest <- object$fit$coefficients$countcoefest[,c(5,6,7)]
        colnames(countcoefest) <- c("Estimate", "2.5 %", "97.5 %")
        zerocoefest <- object$fit$coefficients$zerocoefest[,c(5,6,7)]
        colnames(zerocoefest) <- c("Estimate", "2.5 %", "97.5 %")
      }
      coefest <- list(countcoefest = countcoefest, zerocoefest = zerocoefest)
    }
    else{
      if(sumtype == "norm"){
        coefest <- object$fit$coefficients[,c(1,2,3,4)]
        colnames(coefest) <- c("Estimate", "Std. Error", "2.5 %", "97.5 %")
      }
      else if(sumtype == "perc"){
        coefest <- object$fit$coefficients[,c(5,6,7)]
        colnames(coefest) <- c("Estimate", "2.5 %", "97.5 %")
      }
    }

    ans <- list(coefficients = coefest, call = object$call, family = object$family, zero_infl = object$zero_infl,
                dispersion = object$fit$dispersion, deviance = object$fit$deviance, df.residual = object$fit$df.residual,
                null.deviance = object$fit$null.deviance, df.null = object$fit$df.null, aic = object$fit$aic, theta = object$fit$theta,
                SE.theta = object$fit$SE.theta, zilink = object$zilink, rh = object$rh)
    ans$is.binomial <- FALSE
    class(ans) <- "summary.gwqs"
    return(ans)
  }
}

#' @rawNamespace S3method(print, gwqs)
#' @rdname methods

print.gwqs <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  if(x$rh==1){
    if(x$family$family == "multinomial"){
      if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
      }
      cat("\nCoefficients:\n")
      print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
      cat("\nResidual Deviance:", format(2*x$fit$nlm_out$minimum), "\n")
      cat("AIC:", format(2*dim(x$fit$coefficients)[1] + 2*(x$fit$nlm_out$minimum)), "\n")
      invisible(x)
    }
    else{
      x$fit$call <- x$call
      print(x$fit, digits, ...)
    }
  }
  else{
    if(!is.null(cl <- x$call)) {
      cat("Call:\n")
      dput(cl, control = NULL)
    }
    cat("\nMean Coefficients:\n")
    print.default(format(coef(x, sumtype = "norm"), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nMedian Coefficients:\n")
    print.default(format(coef(x, sumtype = "perc"), digits = digits), print.gap = 2, quote = FALSE)
    if(x$family$family == "multinomial") cat("\nMean Residual Deviance:", format(2*x$fit$deviance), "\n")
    else if(!x$zero_infl){
      cat("\n")
      cat(apply(cbind(paste(format(c("Mean null","Mean residual"), justify="right"),
                            "deviance:"),
                      format(c(x$fit$null.deviance, x$fit$deviance),
                             digits = max(5L, digits + 1L)), " on",
                      format(c(x$fit$df.null, x$fit$df.residual)),
                      " degrees of freedom\n"),
                1L, paste, collapse = " "), sep = "")
    }
    cat("\nMean AIC: ", format(x$fit$aic, digits = max(4L, digits + 1L)),"\n")
    cat("\n")
    invisible(x)
    if(x$family$family == "negbin"){
      dp <- max(2 - floor(log10(x$fit$SE.theta)), 0)
      cat("Mean Theta: ", format(round(x$fit$theta, dp), nsmall = dp), "\nStd. Err.: ",
          format(round(x$fit$SE.theta, dp), nsmall = dp), "\n")
      invisible(x)
    }
  }
}

#' @rawNamespace S3method(print, summary.gwqs)
#' @rdname methods

print.summary.gwqs <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  if(x$rh==1){
    if(x$family$family == "multinomial"){
      if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
      }
      cat("\nCoefficients:\n")
      if(x$is.binomial) {
        print(cbind(Values = x$coefficients,
                    "Std. Err." = x$standard.errors,
                    "Value/SE" = x$Wald.ratios),
              digits = digits)
      } else {
        print(x$coefficients, digits = digits)
        cat("\nStd. Errors:\n")
        print(x$standard.errors, digits = digits)
        if(!is.null(x$Wald.ratios)) {
          cat("\nValue/SE (Wald statistics):\n")
          print(x$coefficients/x$standard.errors, digits = digits)
        }
      }
      cat("\nResidual Deviance:", format(x$deviance), "\n")
      cat("AIC:", format(x$AIC), "\n")
      if(!is.null(correl <- x$correlation)) {
        p <- dim(correl)[2L]
        if(p > 1) {
          cat("\nCorrelation of Coefficients:\n")
          ll <- lower.tri(correl)
          correl[ll] <- format(round(correl[ll], digits))
          correl[!ll] <- ""
          print(correl[-1L, -p], quote = FALSE, ...)
        }
      }
      invisible(x)
    }
    else{
      if(x$zero_infl){
        cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

        if(!x$converged) {
          cat("model did not converge\n")
        } else {

          cat("Pearson residuals:\n")
          print(structure(quantile(x$residuals),
                          names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)

          cat(paste("\nCount model coefficients (", x$dist, " with log link):\n", sep = ""))
          printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)

          cat(paste("\nZero-inflation model coefficients (binomial with ", x$link, " link):\n", sep = ""))
          printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)

          if(getOption("show.signif.stars") & any(rbind(x$coefficients$count, x$coefficients$zero)[,4] < 0.1))
            cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

          if(x$dist == "negbin") cat(paste("\nTheta =", round(x$theta, digits), "\n")) else cat("\n")
          cat(paste("Number of iterations in", x$method, "optimization:", tail(na.omit(x$optim$count), 1), "\n"))
          cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
        }

        invisible(x)
      }
      else{
        cat("\nCall:\n",
            paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
        cat("Deviance Residuals: \n")
        if(x$df.residual > 5) {
          x$deviance.resid <- setNames(quantile(x$deviance.resid, na.rm = TRUE),
                                       c("Min", "1Q", "Median", "3Q", "Max"))
        }
        xx <- zapsmall(x$deviance.resid, digits + 1L)
        print.default(xx, digits = digits, na.print = "", print.gap = 2L)

        if(length(x$aliased) == 0L) {
          cat("\nNo Coefficients\n")
        } else {
          df <- if ("df" %in% names(x)) x[["df"]] else NULL
          if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
          else cat("\nCoefficients:\n")
          coefs <- x$coefficients
          if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L,
                            dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
          }
          printCoefmat(coefs, digits = digits, signif.stars = getOption("show.signif.stars"),
                       na.print = "NA", ...)
        }
        cat("\n(Dispersion parameter for ", x$family$family,
            " family taken to be ", format(x$dispersion), ")\n\n",
            apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                              "deviance:"),
                        format(unlist(x[c("null.deviance","deviance")]),
                               digits = max(5L, digits + 1L)), " on",
                        format(unlist(x[c("df.null","df.residual")])),
                        " degrees of freedom\n"),
                  1L, paste, collapse = " "), sep = "")
        if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
        cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n\n",
            "Number of Fisher Scoring iterations: ", x$iter,
            "\n", sep = "")

        correl <- x$correlation
        if(!is.null(correl)) {
          p <- NCOL(correl)
          if(p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if(is.logical(x$symbolic.cor) && x$symbolic.cor) {
              print(symnum(correl, abbr.colnames = NULL))
            } else {
              correl <- format(round(correl, 2L), nsmall = 2L,
                               digits = digits)
              correl[!lower.tri(correl)] <- ""
              print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
          }
        }
        cat("\n")
        invisible(x)
        if(x$family$family == "negbin"){
          dp <- max(2 - floor(log10(x$SE.theta)), 0)
          cat("\n              Theta: ", format(round(x$theta, dp), nsmall = dp), "\n          Std. Err.: ",
              format(round(x$SE.theta, dp), nsmall = dp), "\n")
          if (!is.null(x$th.warn))
            cat("Warning while fitting theta:", x$th.warn, "\n")
          cat("\n 2 x log-likelihood: ", format(round(x$twologlik, 3), nsmall = dp), "\n")
          invisible(x)
        }
      }
    }
  }
  else{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if(x$zero_infl){
      cat(paste("\nCount model coefficients (", x$family$family, " with log link):\n", sep = ""))
      printCoefmat(x$coefficients$countcoefest, digits = digits, signif.legend = FALSE)
      cat(paste("\nZero-inflation model coefficients (binomial with ", x$zilink, " link):\n", sep = ""))
      printCoefmat(x$coefficients$zerocoefest, digits = digits, signif.legend = FALSE)
    }
    else{
      cat("\nCoefficients:\n")
      coefs <- x$coefficients
      printCoefmat(coefs, digits = digits, signif.stars = getOption("show.signif.stars"),
                   na.print = "NA", ...)
    }
    if(x$family$family == "multinomial"){
      cat("\n(Mean dispersion parameter for ", x$family$family,
          " family taken to be ", format(x$dispersion), ")\n\n",
          apply(cbind(paste("Mean residual deviance:"),
                      format(x["deviance"], digits = max(5L, digits + 1L)), "\n"),
                1L, paste, collapse = " "), sep = "")
    }
    else if(!x$zero_infl){
      cat("\n(Mean dispersion parameter for ", x$family$family,
          " family taken to be ", format(x$dispersion), ")\n\n",
          apply(cbind(paste(format(c("Mean null","Mean residual"), justify="right"),
                            "deviance:"),
                      format(unlist(x[c("null.deviance","deviance")]),
                             digits = max(5L, digits + 1L)), " on",
                      format(unlist(x[c("df.null","df.residual")])),
                      " degrees of freedom\n"),
                1L, paste, collapse = " "), sep = "")
    }
    cat("\nMean AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n")
    cat("\n")
    invisible(x)
    if(x$family$family == "negbin"){
      dp <- max(2 - floor(log10(x$SE.theta)), 0)
      cat("Mean Theta: ", format(round(x$theta, dp), nsmall = dp), "\nStd. Err.: ",
          format(round(x$SE.theta, dp), nsmall = dp), "\n")
      if (!is.null(x$th.warn))
        cat("Warning while fitting theta:", x$th.warn, "\n")
      invisible(x)
    }
  }
}

#' @rawNamespace S3method(predict, gwqs)
#' @rdname methods

predict.gwqs <- function(object, newdata, sumtype = c("norm", "perc"), type=c("link", "response", "prob", "count", "zero", "class", "probs", "terms"), ...){
  if(object$rh==1){
    type <- match.arg(type)
    if(missing(newdata)){
      data <- object$data[object$validation_rows,]
      data <- cbind(data, object$wqs)
      if(object$family$family == "multinomial") predictmultinom(object, data, type)$pred
      else predict(object$fit, type = type)
    }
    else{
      data <- newdata
      if (is.null(object$q)) {Q <- as.matrix(data[, object$mix_name]); qi <- NULL}
      else {
        Ql <- lapply(1:length(object$mix_name), function(i) cut(unlist(data[,object$mix_name[i]]), breaks = object$qi[[i]], labels = FALSE, include.lowest = TRUE) - 1)
        Q <- do.call("cbind", Ql)
      }
      qi <- object$qi
      if(!is.null(object$stratified)){
        strtfd_out <- stratified_f(Q, data, object$stratified, object$mix_name)
        Q <- strtfd_out$Q
        object$mix_name <- strtfd_out$mix_name
      }
      if(object$family$family == "multinomial"){
        wqs <- Q%*%as.matrix(object$final_weights[match(object$mix_name, object$final_weights[,1]),-1])
        data <- cbind(data, wqs)
        predm <- predictmultinom(object, data, NULL, type)
        pred <- predm$pred
        y <- predm$y
      }
      else{
        wqs <- data$wqs <- as.numeric(Q%*%object$final_weights$mean_weight[match(object$mix_name, as.character(object$final_weights$mix_name))])
        pred <- predict(object$fit, newdata = data, type = type)
        y <- model.response(model.frame(object$formula, data), "any")
      }
      df_pred <- data.frame(y = y, ypred = pred)
      return(list(df_pred = df_pred, Q = Q, qi = qi, wqs = wqs))
    }
  }
  else{
    type <- match.arg(type)
    sumtype <- match.arg(sumtype)
    if(missing(newdata)){
      data <- object$data
      if(!object$dwqs) data <- cbind(data, object$wqs)
      if(object$family$family == "multinomial") predictmultinom(object, data, sumtype, type)$pred
      else{
        if(object$zero_infl){
          object$gwqslist[[1]]$fit$coefficients$count <- switch(sumtype, norm = object$fit$coefficients$countcoefest[,1], perc = object$fit$coefficients$countcoefest[,5])
          object$gwqslist[[1]]$fit$coefficients$zero <- switch(sumtype, norm = object$fit$coefficients$zerocoefest[,1], perc = object$fit$coefficients$zerocoefest[,5])
        }
        else object$gwqslist[[1]]$fit$coefficients <- coef(object, sumtype)
        predict(object$gwqslist[[1]]$fit, newdata = object$data, type = type)
      }
    }
    else{
      data <- newdata
      if (is.null(object$q)) {Q <- as.matrix(data[, object$mix_name]); qi <- NULL}
      else {
        Ql <- lapply(1:length(object$mix_name), function(i) cut(unlist(data[,object$mix_name[i]]), breaks = object$qi[[i]], labels = FALSE, include.lowest = TRUE) - 1)
        Q <- do.call("cbind", Ql)
      }
      qi <- object$qi
      if(!is.null(object$stratified)){
        strtfd_out <- stratified_f(Q, data, object$stratified, object$mix_name)
        Q <- strtfd_out$Q
        object$mix_name <- strtfd_out$mix_name
      }
      if(object$family$family == "multinomial"){
        wqs <- Q%*%as.matrix(object$final_weights[,-1][match(object$mix_name, object$final_weights[,1]), c(3*0:(object$n_levels-2)+1)])
        data <- cbind(data, wqs)
        predm <- predictmultinom(object, data, sumtype, type)
        pred <- predm$pred
        y <- predm$y
      }
      else{
        wqs <- data$wqs <- as.numeric(Q%*%object$final_weights$Estimate[match(object$mix_name, as.character(object$final_weights$mix_name))])
        object$gwqslist[[1]]$fit$coefficients <- coef(object, sumtype)
        pred <- predict(object$gwqslist[[1]]$fit, newdata = data, type = type)
        y <- model.response(model.frame(object$formula, data), "any")
      }
      df_pred <- data.frame(y = y, ypred = pred)
      return(list(df_pred = df_pred, Q = Q, qi = qi, wqs = wqs))
    }
  }
}

#' @rawNamespace S3method(coef, gwqs)
#' @rdname methods

coef.gwqs <- function(object, sumtype = c("norm", "perc"), ...){
  if(object$rh==1){
    if(object$family$family == "multinomial"){
      coef <- matrix(object$fit$coefficients$Estimate, nrow = dim(object$fit$coefficients)[1]/2, byrow = T)
      rownames(coef) <- levels(unlist(object$data[, all.vars(object$formula)[1]]))[-1]
      object$data$wqs <- 0
      colnames(coef) <- colnames(model.matrix(object$formula, object$data))
      coef
    }
    else coef(object$fit)
  }
  else{
    sumtype <- match.arg(sumtype)
    i <- ifelse(sumtype == "norm", 1, ifelse(sumtype == "perc", 5, NA))
    if(is.na(i)) stop("sumtype can be equal to norm or perc.\n")
    if(object$family$family == "multinomial"){
      coef <- matrix(object$fit$coefficients[,i], nrow = dim(object$fit$coefficients)[1]/2, byrow = T)
      rownames(coef) <- levels(unlist(object$data[, all.vars(object$formula)[1]]))[-1]
      object$data$wqs <- 0
      colnames(coef) <- colnames(model.matrix(object$formula, object$data))
      coef
    }
    else if(object$zero_infl){
      countcoef <- object$fit$coefficients$countcoefest[,i]
      names(countcoef) <- paste0("count_", names(countcoef))
      zerocoef <- object$fit$coefficients$zerocoefest[,i]
      names(zerocoef) <- paste0("zero_", names(zerocoef))
      c(countcoef, zerocoef)
    }
    else object$fit$coefficients[,i]
  }
}

#' @rawNamespace S3method(vcov, gwqs)
#' @rdname methods

vcov.gwqs <- function(object, model = c("full", "count", "zero"), ...){
  if(object$rh==1){
    if(object$family$family == "multinomial"){
      ans <- solve(object$fit$nlm_out$hessian)
      colnames(ans) <- rownames(ans) <- rownames(object$fit$coefficients)
      ans
    }
    else if(object$zero_infl){
      model <- match.arg(model)
      vcov(object$fit, model)
    }
    else vcov(object$fit)
  }
  else{
    if(object$zero_infl){
      model <- match.arg(model)
      if(model == "full"){
        vcovmat <- var(cbind(object$coefmat$countcoefmat, object$coefmat$zerocoefmat))
        rownames(vcovmat) <- colnames(vcovmat) <- c(paste0("count_", colnames(object$coefmat$countcoefmat)),
                                                    paste0("zero_", colnames(object$coefmat$zerocoefmat)))
        vcovmat
      }
    }
    else var(object$coefmat)
  }
}

#' @rawNamespace S3method(fitted, gwqs)
#' @rdname methods

fitted.gwqs <- function(object, sumtype = c("norm", "perc"), type = c("link", "response", "prob", "count", "zero", "class", "probs", "terms"), ...){
  type <- match.arg(type)
  sumtype <- match.arg(sumtype)
  predict.gwqs(object, sumtype = sumtype, type = type)
}

#' @rawNamespace S3method(residuals, gwqs)
#' @rdname methods

residuals.gwqs <- function(object, sumtype = c("norm", "perc"), type = c("deviance", "pearson", "working", "response", "partial"), ...){
  if(object$rh==1){
    type <- match.arg(type)
    if(object$family$family == "multinomial"){
      if(!(type %in% c("pearson", "response"))) stop("If family is \"multinomial\" then residuals type must be \"response\" or \"pearson\"\n")
      else{
        y <- unlist(object$data[object$validation_rows, as.character(object$formula[[2]])])
        res <- mapply(function(i) i == y, levels(y)) - predict.gwqs(object, type = "prob")
        if(type == "pearson") res <- res/sqrt(predict.gwqs(object, type = "prob"))
        res
      }
    }
    else residuals(object$fit, type)
  }
  else{
    type <- match.arg(type)
    sumtype <- match.arg(sumtype)
    y <- unlist(object$data[, as.character(object$formula[[2]])])
    if(object$family$family == "multinomial"){
      if(!(type %in% c("pearson", "response"))) stop("If family is \"multinomial\" then residuals type must be \"response\" or \"pearson\"\n")
      else{
        res <- mapply(function(i) i == y, levels(y)) - predict.gwqs(object, sumtype = sumtype, type = "prob")
        if(type == "pearson") res <- res/sqrt(predict.gwqs(object, sumtype = sumtype, type = "prob"))
      }
    }
    else{
      res <- y - predict.gwqs(object, sumtype = sumtype, type = "response")
      if(object$zero_infl & type == "pearson"){
        mu <- predict(object, sumtype = sumtype, type = "count")
        phi <- predict(object, sumtype = sumtype, type = "zero")
        theta1 <- switch(object$family$family, poisson = 0, negbin = 1/object$fit$theta)
        vv <- fitted(object, sumtype = sumtype, type = "response")*(1 + (phi + theta1) * mu)
        res <- res/sqrt(vv)
      }
      else if(type == "pearson") res <- res/sqrt(object$family$variance(predict.gwqs(object, sumtype = sumtype, type = "response")))
    }
    res
  }
}

#' Plots and tables functions
#'
#' Functions that allow to generate plots and tables helping in visualizing and summarise Weighted Quantile Sum (WQS) regression results.
#'
#' @param object An object of class "gwqs" as returned by \link[gWQS]{gwqs}.
#' @param tau A number identifying the cutoff for the significant weights. Is tau is missing then reciprocal of
#' the number of elements in the mixture is considered. To avoid printing the threshold line set \code{tau = NULL}.
#' @param sumtype Type of summary statistic to be used: "norm" takes the mean of the estimated parameters on the
#' validation sets and the 95% CI assuming a normal distribution of the parameters, while "perc" uses the median
#' as the parameters estimates and the 2.5, 97.5 percentiles as CI. This option is only available for objects of
#' class \code{gwqsrh}.
#' @param newdata A data frame in which to look for variables with which to predict and generate the
#' ROC curve.
#' @param data Dataset from which you want to select the variables you are interested in.
#' @param na.action Allows to choose what action has to be taken to deal with NAs.
#' @param formula Formula used in the model to specify the dependent and independent variables.
#' @param mix_name Vector containing element names included in the mixture.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}).
#' @param ... Further arguments to be passed.
#'
#' @details
#' The \code{gwqs_barplot}, \code{gwqs_scatterplot}, \code{gwqs_fitted_vs_resid}, \code{gwqs_levels_scatterplot},
#' \code{gwqs_ROC} and \code{gwqs_boxplot} functions produce five figures through the \code{\link[ggplot2]{ggplot}} function.
#'
#' The \code{gwqs_summary_tab} and \code{gwqs_weights_tab} functions produce two tables in the viewr pane
#' through the use of the \code{\link[knitr]{kable}} and \code{\link[kableExtra]{kable_styling}} functions.
#'
#' The \code{gwqs_barplot}, \code{gwqs_scatterplot} plots are available for all family types while
#' \code{gwqs_fitted_vs_resid} is not available when \code{family = binomial} or \code{"multinomial"}.
#' \code{gwqs_levels_scatterplot} plot is only available when \code{family = "multinomial"} and \code{gwqs_ROC}
#' when \code{family = binomial}. The \code{gwqs_boxplot} can be used when the parameter \code{rh} within
#' the \code{gwqs} function is set greater than 1.
#'
#' The \code{gwqs_rank} function allows to split the variables selected through the vector \code{mix_name}
#' in quantiles (depending by the value assigned to \code{q}).
#'
#' @return All the plot functions print the output in the Plots pane while the table functions print
#' the output in the Viewer pane.
#'
#' \item{Qm}{The matrix containing the quantiled variables of the elements included in the mixture.}
#' \item{qi}{A list of vectors containing the cut points used to determine the quantiled variables.}
#'
#' @author
#' Stefano Renzetti, Paul Curtin, Allan C Just, Ghalib Bello, Chris Gennings
#'
#' @examples
#' toxic_chems = names(wqs_data)[1:34]
#' results = gwqs(yLBX ~ wqs, mix_name = toxic_chems, data = wqs_data, q = 4, validation = 0.6,
#'                b = 2, b1_pos = TRUE, b_constr = FALSE, family = gaussian)
#'
#' # barplot
#' gwqs_barplot(results)
#'
#' # scatterplot
#' gwqs_scatterplot(results)
#'
#' # fitted values vs rediduals scatterplot
#' gwqs_fitted_vs_resid(results)
#'
#' @export gwqs_barplot
#' @rdname secondary_functions

# Functions to generate the plots
gwqs_barplot <- function(object, tau, ...){
  if(object$family$family == "multinomial"){
    if(object$rh>1) object$final_weights <- object$final_weights[,c(1, 3*0:(object$n_levels-2)+2)]
    data_plot <- object$final_weights[order(object$final_weights[, object$levelnames[1]]),]
    pos <- match(data_plot$mix_name, sort(object$mix_name))
    data_plot$mix_name <- factor(data_plot$mix_name, levels(data_plot$mix_name)[pos])
    data_plot_l <- melt(data_plot, id.vars = "mix_name")
    bar_plot_h <- ggplot(data_plot_l, aes_string(x = "mix_name", y = "value")) +
      facet_wrap(~ variable)
  }
  else {
    if(object$rh>1){
      if(object$dwqs) object$final_weights <- melt(object$final_weights[,c(1,2,5)])[,c(1,3,2)]
      else object$final_weights <- object$final_weights[,c(1, 2)]
    }
    else if(object$rh==1 & object$dwqs) object$final_weights <- melt(object$final_weights)[,c(1,3,2)]
    names(object$final_weights)[1:2] <- c("mix_name", "mean_weight")
    data_plot <- object$final_weights[order(object$final_weights$mean_weight),]
    if(!is.null(object$stratified)){
      data_plot$strata <- sapply(data_plot$mix_name, function(i){
        l <- levels(object$data[,object$stratified])
        l[sapply(l, function(j) grepl(paste0("_", j, "$"), i))]
      })
      data_plot$mix_name <- sapply(1:nrow(data_plot), function(i) sub(paste0("_", data_plot$strata[i], "$"), "", data_plot$mix_name[i]))
      data_plot <- data_plot[order(paste0(data_plot$strata, data_plot$mean_weight)),]
    }
    if(object$dwqs){
      data_plot <- data_plot[order(paste0(data_plot$variable, data_plot$mean_weight), decreasing = T),]
      data_plot$mix_name <- factor(data_plot$mix_name, levels = unique(data_plot$mix_name)[(nrow(data_plot)/2):1])
    }
    else data_plot$mix_name <- factor(data_plot$mix_name, levels = unique(data_plot$mix_name))
    bar_plot_h <- ggplot(data_plot, aes_string(x = "mix_name", y = "mean_weight"))
  }

  bar_plot_h <- bar_plot_h + geom_bar(stat = "identity", color = "black") + theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(color='black'),
          legend.position = "none") + coord_flip()
  if(missing(tau)) tau <- 1/length(object$mix_name)
  if(!is.null(tau)) bar_plot_h <- bar_plot_h + geom_hline(yintercept = tau, linetype="dashed", color = "red")
  if(object$dwqs) bar_plot_h <- bar_plot_h + facet_wrap(~ variable, nrow = 1)
  if(!is.null(object$stratified)) bar_plot_h <- bar_plot_h + facet_wrap(~ strata, nrow = 1)

  print(bar_plot_h)
}

#' @export gwqs_scatterplot
#' @rdname secondary_functions

gwqs_scatterplot <- function(object, ...){
  y_labs = ifelse(object$family$family %in% c("multinomial", "binomial"), "y", "y_adj")

  if(object$dwqs) object$y_wqs_df <- melt(object$y_wqs_df, measure.vars = c("pwqs", "nwqs"), value.name = "wqs")

  yadj_vs_wqs = ggplot(object$y_wqs_df, aes_string("wqs", y_labs)) +
    geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()

  if(object$dwqs) yadj_vs_wqs <- yadj_vs_wqs + facet_wrap(~ variable)

  if(object$family$family == "multinomial") yadj_vs_wqs = yadj_vs_wqs + facet_wrap(~ level)

  print(yadj_vs_wqs)
}

#' @export gwqs_fitted_vs_resid
#' @rdname secondary_functions

gwqs_fitted_vs_resid <- function(object, sumtype = c("norm", "perc"), ...){
  sumtype <- match.arg(sumtype)
  if(!(object$family$family %in% c("binomial", "multinomial"))){
    # if(object$zero_infl) fit_df = data.frame(fitted = fitted(object), resid = residuals(object, type = "response"))
    # if(object$zero_infl) fit_df = data.frame(.fitted = object$fit$fitted.values, .resid = object$fit$residuals)
    # else fit_df = augment(object$fit)
    if(object$rh==1) fit_df = data.frame(fitted = fitted(object), resid = residuals(object, type = "response"))
    else if(object$rh>1) fit_df = data.frame(fitted = fitted(object, sumtype, type = "response"), resid = residuals(object, sumtype = sumtype, type = "response"))
    res_vs_fitted = ggplot(fit_df, aes_string(x = "fitted", y = "resid")) + geom_point() + theme_bw() +
      xlab("Fitted values") + ylab("Residuals")
    print(res_vs_fitted)
  }
  else stop("This plot is not available for binomial or multinomial family\n")
}

#' @export gwqs_levels_scatterplot
#' @rdname secondary_functions

gwqs_levels_scatterplot <- function(object, ...){
  if(object$n_levels == 3){
    if(object$rh>1) object$final_weights <- object$final_weights[,c(1, 3*0:(object$n_levels-2)+2)]
    dataplot_names <- names(object$final_weights)
    w1_vs_w2 = ggplot(object$final_weights, aes_string(dataplot_names[2], dataplot_names[3])) + geom_point() +
      theme_bw() + xlab(dataplot_names[2]) + ylab(dataplot_names[3]) + geom_abline(linetype = 2) +
      geom_text_repel(aes(label = object$mix_name))
    print(w1_vs_w2)
  }
  else stop("This plot is only available if family is multinomial and the dependent variable has 3 levels\n")
}

#' @export gwqs_ROC
#' @rdname secondary_functions

gwqs_ROC <- function(object, newdata, sumtype = c("norm", "perc"), ...){
  sumtype <- match.arg(sumtype)
  if(missing(newdata)) stop("argument \"newdata\" is missing, with no default\n")
  if(object$family$family == "binomial"){
    if(object$rh==1) wqs_pred <- predict(object, newdata, type = "response")
    else if(object$rh>1) wqs_pred <- predict(object, newdata, sumtype, type = "response")
    df_roc <- wqs_pred$df_pred
    if(inherits(df_roc$y, "character")) df_roc$y = factor(df_roc$y)
    if(inherits(df_roc$y, "factor")) df_roc$y <- as.numeric(df_roc$y != levels(df_roc$y)[1])
    gg_roc = suppressWarnings(ggplot(df_roc, aes_string(d="y", m="ypred")) + geom_roc(n.cuts = 0) +
                                style_roc(xlab = "1 - Specificity", ylab = "Sensitivity"))
    auc_est = calc_auc(gg_roc)
    gg_roc = gg_roc + annotate("text", x=0.75, y=0.25, label=paste0("AUC = ", round(auc_est[, "AUC"], 3)))

    print(gg_roc)
  }
  else stop("The ROC curve is only available for binomial family\n")
}

#' @export gwqs_boxplot
#' @rdname secondary_functions

gwqs_boxplot <- function(object, tau, ...){
  if(object$rh>1){
    if(object$family$family == "multinomial"){
      wboxplotl <- lapply(object$levelnames, function(i){
        tmp <- melt(object$wmat[[i]], varnames = c("rh", "mix_name"))
        tmp$level <- i
        return(tmp)
      })
      wboxplot <- do.call("rbind", wboxplotl)
    }
    else wboxplot <- melt(object$wmat, varnames = c("rh", "mix_name"))
    if(!is.null(object$stratified)){
      wboxplot$strata <- sapply(wboxplot$mix_name, function(i){
        l <- levels(object$data[,object$stratified])
        l[sapply(l, function(j) grepl(paste0("_", j, "$"), i))]
      })
      wboxplot$mix_name <- sapply(1:nrow(wboxplot), function(i) sub(paste0("_", wboxplot$strata[i], "$"), "", wboxplot$mix_name[i]))
      tmp <- object$final_weights
      tmp$strata <- sapply(tmp$mix_name, function(i){
        l <- levels(object$data[,object$stratified])
        l[sapply(l, function(j) grepl(paste0("_", j, "$"), i))]
      })
      tmp$mix_name <- sapply(1:nrow(tmp), function(i) sub(paste0("_", tmp$strata[i], "$"), "", tmp$mix_name[i]))
      tmp <- tmp[order(paste0(tmp$strata, tmp$mean_weight)),]
      wboxplot$mix_name <- factor(wboxplot$mix_name, levels = unique(tmp$mix_name))
    }
    else wboxplot$mix_name <- factor(wboxplot$mix_name, levels = object$final_weights$mix_name)
    if(object$dwqs) wboxplot$L1 <- factor(wboxplot$L1, levels = c("wmatpos", "wmatneg"), labels = c("pwqs", "nwqs"))
    box_plot <- ggplot(wboxplot, aes_string(x = "mix_name", y = "value")) +
      geom_boxplot(outlier.shape = " ") + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Weight") +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3) + geom_jitter(alpha = 0.3)
    if(object$family$family == "multinomial") box_plot <- box_plot + facet_wrap(~level)
    if(missing(tau)) tau <- 1/length(object$mix_name)
    if(!is.null(tau)) box_plot <- box_plot + geom_hline(yintercept = tau, linetype="dashed", color = "red")
    if(object$dwqs) box_plot <- box_plot + facet_wrap(~ L1, ncol = 1)
    if(!is.null(object$stratified)) box_plot <- box_plot + facet_wrap(~ strata, ncol = 1)

    print(box_plot)
  }
  else stop("The function gwqs_boxplot is only available when applying a repeated holdout validation procedure (setting the argument rh > 1)\n")
}

#' @export gwqs_summary_tab
#' @rdname secondary_functions

# functions to generate tables
gwqs_summary_tab <- function(object, sumtype = c("norm", "perc"), ...){
  sumtype <- match.arg(sumtype)
  if(object$rh==1){
    if(object$family$family == "multinomial") mf_df = signif(object$fit$coefficients, 3)
    else{
      if(object$zero_infl){
        mf_df_count = as.data.frame(signif(coef(summary(object$fit))$count, 3))
        mf_df_zero = as.data.frame(signif(coef(summary(object$fit))$zero, 3))
        rownames(mf_df_zero) = paste0("z_", rownames(mf_df_zero))
        mf_df = rbind(mf_df_count, mf_df_zero)
      }
      else mf_df = as.data.frame(signif(coef(summary(object$fit)), 3))
    }
  }
  else if(object$rh>1){
    if(object$zero_infl){
      mf_df_count = as.data.frame(signif(summary(object, sumtype = sumtype)$coefficients$count, 3))
      mf_df_zero = as.data.frame(signif(summary(object, sumtype = sumtype)$coefficients$zero, 3))
      rownames(mf_df_zero) = paste0("z_", rownames(mf_df_zero))
      mf_df = rbind(mf_df_count, mf_df_zero)
    }
    else mf_df = as.data.frame(signif(summary(object, sumtype = sumtype)$coefficients, 3))
  }
  if(object$zero_infl) print(kable_styling(kable(mf_df, row.names = TRUE, ...)) %>%
                        group_rows("Count model", 1, dim(mf_df_count)[1]) %>%
                        group_rows("Zero-inflation model", dim(mf_df_count)[1]+1, dim(mf_df_count)[1]+dim(mf_df_zero)[1]))
  else print(kable_styling(kable(mf_df, row.names = TRUE, ...)))
}

#' @export gwqs_weights_tab
#' @rdname secondary_functions

gwqs_weights_tab <- function(object, ...){
  final_weight <- object$final_weights
  final_weight[, -1] <- signif(final_weight[, -1], 3)
  print(kable_styling(kable(final_weight, row.names = FALSE, ...)))
}

#' @export selectdatavars
#' @rdname secondary_functions

# function to select variables
selectdatavars <- function(data, na.action, formula, mix_name, ...){
  allvars = all.vars(formula)
  other_vars <- c(...)
  data$wqs = 0
  data$pwqs = 0
  data$nwqs = 0
  data <- data[, c(allvars, mix_name, other_vars)]
  if(missing(na.action)) na.action <- na.omit
  dtf <- na_action(data, na.action)
  return(dtf)
}

#' @export gwqs_rank
#' @rdname secondary_functions

# function to create variables with quantile of the components
gwqs_rank <- function(data, mix_name, q){

  if(!is.numeric(q)) stop("'q' must be a number\n")
  Ql <- lapply(1:length(mix_name), function(i){
    qi <- unique(quantile(data[[mix_name[i]]], probs = seq(0, 1, by = 1/q), na.rm = TRUE))
    if(length(qi) == 1) qi = c(-Inf, qi)
    else{
      qi[1] <- -Inf
      qi[length(qi)] <- Inf
    }
    q <- cut(data[[mix_name[i]]], breaks = qi, labels = FALSE, include.lowest = TRUE) - 1
    return(list(qi, q))
  })
  qi <- lapply(Ql, function(x) x[[1]])
  Qm <- matrix(unlist(lapply(Ql, function(x) x[[2]])), ncol = length(Ql))
  colnames(Qm) <- names(qi) <- mix_name

  qf_out <- list(Qm = Qm, qi = qi)
  return(qf_out)
}


