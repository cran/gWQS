#' Fitting Weighted Quantile Sum regression models
#'
#' Fits Weighted Quantile Sum (WQS) regression, a random subset inplementation of WQS and a repeated holdout validation WQS for continuous, binomial, multinomial, poisson, quasi-poisson and negative binomial outcomes.
#'
#' @param formula An object of class \code{formula} specifying the relationship to be tested. The \code{wqs}
#' term must be included in \code{formula}, e.g. \code{y ~ wqs + ...}. To test for an interaction term with
#' a continuous variable \code{a} or for a quadratic term we can specify the \code{formula} as below:
#' \code{y ~ wqs*a + ...} and \code{y ~ wqs + I(wqs^2) + ...}, respectively.
#' @param data The \code{data.frame} containing the variables to be included in the model.
#' @param na.action \code{\link[stats]{model.frame}}. \code{na.omit} is the default.
#' @param weights An optional vector of weights to be used in the fitting process.
#' Should be \code{NULL} or a numeric vector.
#' @param mix_name A character vector listing the variables contributing to a mixture effect.
#' @param stratified The character name of the variable for which you want to stratify for.
#' It has to be a \code{factor}.
#' @param valid_var A character value containing the name of the variable that identifies the validation
#' and the training dataset. You previously need to create a variable in the dataset which is equal to 1
#' for the observations you want to include in the validation dataset, equal to 0 for the observation
#' you want to include in the training dataset (use 0 also for the validation dataset if you want to train and
#' validate the model on the same data) and equal to 2 if you want to keep part of the data for the
#' predictive model.
#' @param rh Number of repeated holdout validations. This option is only available for \code{gwqsrh} function.
#' @param b Number of bootstrap samples used in parameter estimation.
#' @param b1_pos A logical value that determines whether weights are derived from models where the beta
#' values were positive or negative.
#' @param b1_constr A logial value that determines whether to apply positive (if \code{b1_pos = TRUE}) or
#' negative (if \code{b1_pos = FALSE}) constraints in the optimization function for the weight estimation.
#' @param zero_infl A logical value (\code{TRUE} or \code{FALSE}) that allows to fit a zero inflated
#' model in case \code{family = "poisson"} or \code{family = "negbin"}.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). If \code{q = NULL} then
#' the values of the mixture variables are taken (these must be standardized).
#' @param validation Percentage of the dataset to be used to validate the model. If
#' \code{validation = 0} then the test dataset is used as validation dataset too.
#' @param family A character value that allows to decide for the glm: \code{gaussian} for linear regression,
#' \code{binomial} for logistic regression \code{"multinomial"} for multinomial regression,
#' \code{poisson} for Poisson regression, \code{quasipoisson} for quasi-Poisson regression,
#' \code{"negbin"} for negative binomial regression.
#' @param signal Character identifying the signal function to be used when the average weights
#' are estimated. It can take values from \code{"one"} to apply the identity, \code{"abst"} to apply
#' the absolute value of the t-statistic, \code{"t2"} to apply the squared value of the t-statistic,
#' \code{"expt"} to apply the exponential of the t-statistic as signal function.
#' @param rs A logic value. If \code{rs = FALSE} then the bootstrap implementation of WQS is performed.
#' If \code{rs = TRUE} then the random subset implementation of WQS is applied (see the "Details" and the
#' vignette for further infromation).
#' @param n_vars The number of mixture components to be included at each random subset step.
#' If \code{rs = TRUE} and \code{n_vars = NULL} then the square root of the number of elements
#' in the mixture is taken.
#' @param zilink Character specification of link function in the binary zero-inflation model
#' (you can choose among \code{"logit", "probit", "cloglog", "cauchit", "log"}).
#' @param seed An \code{integer} value to fix the seed, if it is equal to \code{NULL} no seed is chosen.
#' @param plan_strategy A character value that allows to choose the evaluation strategies for the
#' \code{plan} function. You can choose among "sequential", "transparent", "multisession", "multicore",
#' "multiprocess", "cluster" and "remote" (see \code{\link[future]{plan}} help page for more details).
#' @param optim.method A character identifying the method to be used by the \code{\link[stats]{optim}} function
#' (you can choose among \code{"BFGS", "Nelder-Mead", "CG", "SANN"}, \code{"BFGS"} is the default).
#' See \code{\link[stats]{optim}} for details.
#' @param control The control list of optimization parameters. See \code{\link[stats]{optim}} for details.
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
#' The \code{seed} argument specifies a fixed seed through the \code{set.seed} function.\cr
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
#' \item{qi}{List of the cutoffs used to divide in quantiles the variables in the mixture}
#' \item{bindex}{List of vectors containing the \code{rownames} of the subjects included in each
#' bootstrap dataset.}
#' \item{tindex}{Vector containing the rows used to estimate the weights in each bootstrap.}
#' \item{vindex}{Vector containing the rows used to estimate the parameters of the final model.}
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
#' \item{levelnames}{The name of each level when a multinomial regression is ran.}
#' \item{data}{The data used in the WQS analysis.}
#' \item{objfn_values}{The vector of the b values of the objective function corresponding to the optima values}
#' \item{optim_messages}{The vector of character strings giving any additional information returned by the
#' optimizer, or NULL.}
#' \item{gwqslist}{List of the output from the \code{rh} WQS models.}
#' \item{coefmat}{Matrix containing the parameter estimates from each repeated holdout WQS model.}
#' \item{wmat}{Matrix containing the weight estimates from each repeated holdout WQS model.}
#'
#' @author
#' Stefano Renzetti, Paul Curtin, Allan C Just, Ghalib Bello, Chris Gennings
#'
#' @references
#' Carrico C, Gennings C, Wheeler D, Factor-Litvak P. Characterization of a weighted quantile sum
#' regression for highly correlated data in a risk analysis setting. J Biol Agricul Environ Stat.
#' 2014:1-21. ISSN: 1085-7117. DOI: 10.1007/ s13253-014-0180-3.
#' \url{http://dx.doi.org/10.1007/s13253-014-0180-3}.\cr
#'
#' Czarnota J, Gennings C, Colt JS, De Roos AJ, Cerhan JR, Severson RK, Hartge P, Ward MH,
#' Wheeler D. 2015. Analysis of environmental chemical mixtures and non-Hodgkin lymphoma risk in the
#' NCI-SEER NHL study. Environmental Health Perspectives, DOI:10.1289/ehp.1408630.\cr
#'
#' Czarnota J, Gennings C, Wheeler D. 2015. Assessment of weighted quantile sum regression for modeling
#' chemical mixtures and cancer risk. Cancer Informatics,
#' 2015:14(S2) 159-171 DOI: 10.4137/CIN.S17295.\cr
#'
#' Brunst KJ, Sanchez Guerra M, Gennings C, et al. Maternal Lifetime Stress and Prenatal Psychological
#' Functioning and Decreased Placental Mitochondrial DNA Copy Number in the PRISM Study.
#' Am J Epidemiol. 2017;186(11):1227-1236. doi:10.1093/aje/kwx183.\cr
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
#' results = gwqs(y ~ wqs, mix_name = toxic_chems, data = wqs_data, q = 4, validation = 0.6,
#'                b = 2, b1_pos = TRUE, b1_constr = FALSE, family = gaussian, seed = 2016)
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
#' @importFrom plotROC geom_roc style_roc calc_auc
#' @importFrom nnet multinom
#' @importFrom future plan future value
#' @importFrom future.apply future_lapply
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom pscl zeroinfl
#'
#' @export gwqs
#' @rdname main_functions
#'
#' @usage  gwqs(formula, data, na.action, weights, mix_name, stratified, valid_var, b = 100,
#'              b1_pos = TRUE, b1_constr = FALSE, zero_infl = FALSE, q = 4,
#'              validation = 0.6, family = gaussian, signal = c("t2", "one", "abst", "expt"),
#'              rs = FALSE, n_vars = NULL,
#'              zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
#'              plan_strategy = "sequential",
#'              optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
#'              control = list(trace = FALSE, maxit = 1000, reltol = 1e-8), ...)

gwqs <- function(formula, data, na.action, weights, mix_name, stratified, valid_var,
                 b = 100, b1_pos = TRUE, b1_constr = FALSE, zero_infl = FALSE, q = 4,
                 validation = 0.6, family = gaussian,
                 signal = c("t2", "one", "abst", "expt"), rs = FALSE, n_vars = NULL,
                 zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
                 plan_strategy = "sequential",
                 optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
                 control = list(trace = FALSE, maxit = 1000, reltol = 1e-8), ...){

  wqsGaussBin <- function(initp, kw, bdtf, Y, offset, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){

    if(b1_constr) initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
    w <- initp[(kx + 1):(kx + kw)]^2
    w <- w/sum(w)
    bdtf$wqs <- as.numeric(Q%*%w)
    X <- model.matrix(formula, bdtf)
    b_covs <- initp[1:kx]
    term <- as.numeric(X%*%b_covs) + offset
    f = sum(family$dev.resids(y = Y, mu = family$linkinv(term), wt = wghts))

    return(f)
  }

  wqsPoisson <- function(initp, kw, bdtf, Y, offset, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){

    if(b1_constr) initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
    w <- initp[(kx + 1):(kx + kw)]^2
    w <- w/sum(w)
    bdtf$wqs <- as.numeric(Q%*%w)
    X <- model.matrix(formula, bdtf)
    b_covs <- initp[1:kx]
    term <- as.numeric(X%*%b_covs) + offset
    f = -sum(dpois(Y, lambda = exp(term), log = TRUE))

    return(f)
  }

  wqsNegBin <- function(initp, kw, bdtf, Y, offset, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){

    if(b1_constr) initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), -abs(initp["wqs"]))
    w <- initp[(kx + 1):(kx + kw)]^2
    w <- w/sum(w)
    bdtf$wqs <- as.numeric(Q%*%w)
    X <- model.matrix(formula, bdtf)
    b_covs <- initp[1:kx]
    theta <- exp(initp[length(initp)])
    term <- as.numeric(X%*%b_covs) + offset
    f = -sum((suppressWarnings(dnbinom(Y, size = theta, mu = exp(term), log = TRUE)))*wghts)

    return(f)
  }

  wqsMultinom <- function(initp, kw, bdtf, Y, offset, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, wghts, stratified, b1_pos, b1_constr){

    if(b1_constr){
      par_pos <- which(grepl("wqs", names(initp)))
      initp[par_pos] <- sapply(1:length(b1_pos), function(i) ifelse(b1_pos[i], abs(initp[par_pos[i]]), -abs(initp[par_pos[i]])))
    }
    w <- matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels-1)
    w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
    bdtf[, wqsvars] <- Q%*%w
    Xl = lapply(wqsvars, function(i){
      fmi <- as.formula(gsub("wqs", i, format(formula)))
      model.matrix(fmi, data = bdtf)
    })
    X = do.call("cbind", Xl)
    b_covs = matrix(0, kx, n_levels-1)
    i = 1:(kx/(n_levels-1))
    for (j in 0:(n_levels-2)){
      b_covs[(kx/(n_levels-1))*j+i, j+1] = initp[(kx/(n_levels-1))*j+i]
    }
    term = X%*%b_covs + offset
    f = -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*wghts)

    return(f)
  }

  one <- function(x){
    rep(1, length(x))
  }

  t2 <- function(x){
    x^2
  }

  if(is.character(family)){
    if(family %in% c("multinomial", "negbin")) family <- list(family = family)
    else family <- get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized\n")
  }

  objfn <- switch(family$family,
                  "gaussian" = wqsGaussBin,
                  "binomial" = wqsGaussBin,
                  "poisson" = wqsPoisson,
                  "quasipoisson" = wqsPoisson,
                  "negbin" = wqsNegBin,
                  "multinomial" = wqsMultinom)
  optim.method <- match.arg(optim.method)

  signal <- match.arg(signal)
  signal <- switch(signal,
                  "one" = one,
                  "abst" = abs,
                  "t2" = t2,
                  "expt" = exp)

  if(!any(grepl("wqs", rownames(attr(terms(formula), "factors"))))) stop("'formula' must contain 'wqs' term: e.g. y ~ wqs + ...\n")

  if(zero_infl){
    zilink <- make.link(match.arg(zilink))
    ff = formula
    if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
      formula = as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))
  }
  else zilink <- ff <- NULL

  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", "weights", "stratified", "valid_var"), names(mc), 0)
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
  if (m){
    strtfd_out = stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name = strtfd_out$mix_name
  }
  else stratified <- NULL

  if(!is.null(seed) & !is.numeric(seed)) stop("seed must be numeric or NULL\n")
  if(!is.null(seed)){
    set.seed(seed)
    fseed <- TRUE
  }
  else fseed <- FALSE

  N <- nrow(dtf)

  # splitting the dataset
  m <- match(c("valid_var"), names(mc), 0)
  rindex = create_rindex(dtf, N, validation, valid_var, m, family)

  if(is.null(n_vars)) n_vars = round(sqrt(length(mix_name)))

  # parameters estimation and model fitting
  m <- match(c("weights"), names(mc), 0)
  if(m[1]) dtf$wghts <- wghts <- unlist(dtf[, weights, drop = FALSE])
  else  dtf$wghts <- wghts <- rep(1, N)

  if (family$family %in% c("gaussian", "quasipoisson")) ts = "t"
  else if (family$family %in% c("binomial", "poisson", "multinomial", "negbin")) ts = "z"

  if(!is.numeric(b)) stop("'b' must be a number\n")

  plan(plan_strategy)
  if(control$trace) cat("start opt\n")
  param <- future_lapply(X = 1:b, FUN = optim.f, objfn = objfn, dtf = dtf[rindex$it,], bQ = Q[rindex$it,],
                         b1_pos = b1_pos, b1_constr = b1_constr, n_vars = n_vars, family = family, rs = rs,
                         zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff, wghts = wghts[rindex$it],
                         stratified = stratified, optim.method = optim.method, control = control, future.seed = FALSE)

  conv <- c(sapply(param, function(i) i$conv))
  counts <- c(sapply(param, function(i) i$counts))
  val <- c(sapply(param, function(i) i$val))
  mex <- lapply(param, function(i) i$mex)
  bindex <- lapply(param, function(i) i$bindex)
  slctd_vars <- lapply(param, function(i) i$slctd_vars)

  if(rs){
    plan(plan_strategy)
    param <- future_lapply(X = 1:b, FUN = set_par_names, slctd_vars, param, q_name = colnames(Q), family = family,
                           future.seed = FALSE)
  }

  if(family$family == "multinomial"){
    n_levels <- dim(param[[1]]$par_opt)[2]+1
    wqs_site <- which(grepl("^wqs_", rownames(param[[1]]$mfit$m_f$coefficients)))
    wght_matrix <- lapply(1:(n_levels-1), function(j) do.call("rbind", lapply(param, function(i) i$par_opt[,j])))
    b1 <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$nlm_out$estimate[j]))
    se <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$Standard_Error[j]))
    stat <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$stat[j]))
    p_val <- lapply(wqs_site, function(j) sapply(param, function(i) i$mfit$m_f$coefficients$p_value[j]))
  }
  else{
    wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt))
    if(zero_infl){
      b1_count <- sapply(param, function(i) i$mfit$m_f$coefficients$count["wqs"])
      se_count <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$count["wqs", "Std. Error"])
      stat_count <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$count["wqs", paste0(ts, " value")])
      p_val_count <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$count["wqs", gsub("x", ts, "Pr(>|x|)")])
      if("wqs" %in% names(param[[1]]$mfit$m_f$coefficients$zero)){
        b1_zero <- sapply(param, function(i) i$mfit$m_f$coefficients$zero["wqs"])
        se_zero <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$zero["wqs", "Std. Error"])
        stat_zero <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$zero["wqs", paste0(ts, " value")])
        p_val_zero <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients$zero["wqs", gsub("x", ts, "Pr(>|x|)")])
      }
      else b1_zero <- se_zero <- stat_zero <- p_val_zero <- NULL
    }
    else{
      b1 <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", "Estimate"])
      se <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", "Std. Error"])
      stat <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", paste0(ts, " value")])
      p_val <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")])
    }
    n_levels <- 1
  }

  n_non_conv = sum(conv == 1)
  if(n_non_conv == 0 & control$trace) cat(paste0("The optimization function always converged\n"))
  else if(n_non_conv == b) stop("The optimization function never converged\n")
  else if(control$trace) cat(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))
  if(control$trace) cat(paste0("There are ", ifelse(b1_pos, sum(b1 >= 0, na.rm = T), sum(b1 <= 0, na.rm = T)),
                               ifelse(b1_pos, " positive", " negative"), " bootstrapped b1 out of ", b, "\n"))

  # estimate mean weight for each component (exclude weights from iterations with failed convergence)
  if (family$family == "multinomial"){
    bres <- Map(cbind, wght_matrix, b1, se, stat, p_val)
    bres <- lapply(bres, as.data.frame)
    bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", "stat", "p_val"))
    strata_names <- gsub("wqs_", "", rownames(param[[1]]$mfit$m_f$coefficients)[wqs_site])
    names(bres) <- strata_names
  }
  else {
    if(zero_infl){
      if(is.null(b1_zero)){
        bres <- as.data.frame(cbind(wght_matrix, b1_count, se_count, stat_count, p_val_count))
        names(bres) <- c(colnames(Q), "b1_count", "Std_Error_count", "stat_count", "p_val_count")
      }
      else{
        bres <- as.data.frame(cbind(wght_matrix, b1_count, se_count, stat_count, p_val_count, b1_zero, se_zero, stat_zero, p_val_zero))
        names(bres) <- c(colnames(Q), "b1_count", "Std_Error_count", "stat_count", "p_val_count", "b1_zero", "Std_Error_zero", "stat_zero", "p_val_zero")
      }
    }
    else{
      bres <- as.data.frame(cbind(wght_matrix, b1, se, stat, p_val))
      names(bres) <- c(colnames(Q), "b1", "Std_Error", "stat", "p_val")
    }
    strata_names <- NULL
  }

  if (family$family == "multinomial"){
    mean_weight <- lapply(1:(n_levels-1), function(i){
      if(rs) bres[[i]][mix_name][is.na(bres[[i]][mix_name])] <- 0
      if(b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 > 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[[i]][bres[[i]]$b1 > 0 & conv == 0, "stat"]))
      else if(!b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 < 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[[i]][bres[[i]]$b1 < 0 & conv == 0, "stat"]))
      if (all(is.nan(w_t))){
        w_t <- 1/apply(bres[[i]][, mix_name], 2, weighted.mean, signal(bres[[i]][, "stat"]))/sum(1/apply(bres[[i]][, mix_name], 2, weighted.mean, signal(bres[[i]][, "stat"])))
        if(solve_dir_issue) warning(paste0("No models converged in ", ifelse(b1_pos[i], "positive", "negative"), " direction for ", strata_names[i], "\n"))
        else stop(paste0("There are no ", ifelse(b1_pos[i], "positive", "negative"), " b1 in the bootstrapped models for ", strata_names[i], "\n"))
      }
      return(w_t)
    })
    mean_weight <- list.cbind(mean_weight)
  }
  else{
    if(rs) bres[mix_name][is.na(bres[mix_name])] <- 0
    if(zero_infl){
      if(b1_pos) mean_weight = apply(bres[bres$b1_count > 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[bres$b1_count > 0 & conv == 0, "stat_count"]))
      else mean_weight = apply(bres[bres$b1_count < 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[bres$b1_count < 0 & conv == 0, "stat_count"]))
    }
    else{
      if(b1_pos) mean_weight = apply(bres[bres$b1 > 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[bres$b1 > 0 & conv == 0, "stat"]))
      else mean_weight = apply(bres[bres$b1 < 0 & conv == 0, mix_name], 2, weighted.mean, signal(bres[bres$b1 < 0 & conv == 0, "stat"]))
    }
    if(all(is.nan(mean_weight))){
      mean_weight <- 1/apply(bres[, mix_name], 2, weighted.mean, signal(bres[, "stat"]))/sum(1/apply(bres[, mix_name], 2, weighted.mean, signal(bres[, "stat"])))
      if(solve_dir_issue) warning(paste0("No models converged in ", ifelse(b1_pos, "positive", "negative"), " direction\n"))
      else stop("There are no ", ifelse(b1_pos, "positive", "negative"), " b1 in the bootstrapped models\n")
    }
  }

  # fit the final model with the estimated weights
  wqs_model = model.fit(mean_weight, dtf[rindex$iv,], Q[rindex$iv,], family, zilink, formula, ff, wghts[rindex$iv], stratified, b1_pos, zero_infl)

  if(all(grepl("wqs", attr(terms(formula), "term.labels")))) y_plot <- model.response(model.frame(formula, dtf[rindex$iv,]), "any")
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
        fit <- zeroinfl(formula_wo_wqs, dtf[rindex$iv,], dist = family$family, link = zilink$name)
      }
      else{
        if(family$family == "negbin") fit = glm.nb(formula_wo_wqs, dtf[rindex$iv,])
        else fit = glm(formula_wo_wqs, dtf[rindex$iv,], family = family)
      }
      if(family$family == "binomial") y_plot = fit$y
      else y_plot = mean(fit$y) + resid(fit, type = "pearson")
    }
  }

  # Plots
  data_plot <- data.frame(mix_name, mean_weight, stringsAsFactors = TRUE)

  if(family$family == "multinomial"){
    Y <- model.response(model.frame(formula, dtf[rindex$iv,]), "any")
    level_names <- levels(Y)
    names(data_plot)[2:n_levels] = strata_names
    data_plot <- data_plot[order(data_plot[, strata_names[1]], decreasing = TRUE),]
    wqs_index <- wqs_model$wqs
    Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
    colnames(Y) <- strata_names
    y_adj_wqs_df <- do.call("rbind", lapply(strata_names, function(i){
      if(n_levels > 3) data.frame(level = i,
                                  y = Y[rowSums(Y[, -which(colnames(Y)==i)]) == 0, i],
                                  wqs = wqs_model$wqs[rowSums(Y[, -which(colnames(wqs_model$wqs)==i)]) == 0, i],
                                  stringsAsFactors = TRUE)
      else data.frame(level = i,
                      y = Y[Y[, -which(colnames(Y)==i)] == 0, i],
                      wqs = wqs_model$wqs[Y[, -which(colnames(wqs_model$wqs)==i)] == 0, i],
                      stringsAsFactors = TRUE)
    }))
    dtf <- cbind(dtf, Q%*%mean_weight)
    names(dtf)[(ncol(dtf)-ncol(mean_weight)+1):ncol(dtf)] <- colnames(wqs_index)
    dtf$wqs <- NULL
  }
  else{
    y_adj_wqs_df <- as.data.frame(cbind(y_plot, wqs_model$wqs))
    names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", "y", "y_adj"), "wqs")
    data_plot <- data_plot[order(data_plot$mean_weight, decreasing = TRUE),]
    wqs_index <- as.numeric(unlist(wqs_model$wqs))
    dtf$wqs <- as.numeric(Q%*%mean_weight)
  }

  # creating the list of elements to return
  results = list(fit = wqs_model$m_f, final_weights = data_plot, conv = conv, bres = bres, wqs = wqs_index, qi = qi,
                 bindex = bindex, slctd_vars = slctd_vars, tindex = rindex$it, vindex = rindex$iv,
                 y_wqs_df = y_adj_wqs_df, family = family, call = cl, formula = formula, mix_name = mix_name,
                 stratified = stratified, q = q, n_levels = n_levels, zero_infl = zero_infl, zilink = zilink,
                 levelnames = strata_names, data = dtf, objfn_values = val, optim_messages = mex)
  if(zero_infl) results$formula <- ff

  class(results) <- "gwqs"

  return(results)
}



#' @export gwqsrh
#' @rdname main_functions
#'
#' @usage  gwqsrh(formula, data, na.action, weights, mix_name, stratified, valid_var, rh = 100,
#'                b = 100, b1_pos = TRUE, b1_constr = FALSE, zero_infl = FALSE, q = 4,
#'                validation = 0.6, family = gaussian,
#'                signal = c("t2", "one", "abst", "expt"), rs = FALSE, n_vars = NULL,
#'                zilink = c("logit", "probit", "cloglog", "cauchit", "log"), seed = NULL,
#'                plan_strategy = "sequential",
#'                optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
#'                control = list(trace = FALSE, maxit = 1000, reltol = 1e-8), ...)

gwqsrh <- function(formula, data, na.action, weights, mix_name, stratified, valid_var,
                   rh = 100, b = 100, b1_pos = TRUE, b1_constr = FALSE, zero_infl = FALSE,
                   q = 4, validation = 0.6, family = gaussian,
                   signal = c("t2", "one", "abst", "expt"), rs = FALSE, n_vars = NULL,
                   zilink = c("logit", "probit", "cloglog", "cauchit", "log"),
                   seed = NULL, plan_strategy = "sequential",
                   optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"),
                   control = list(trace = FALSE, maxit = 1000, reltol = 1e-8), ...){

  if(is.character(family)){
    if(family %in% c("multinomial", "negbin")) family <- list(family = family)
    else family <- get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized\n")
  }

  if(zero_infl){
    zilink <- match.arg(zilink)
    ff = formula
    if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
      formula = as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))
  }
  else zilink <- ff <- NULL

  mc <- cl <- match.call()
  # mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", "weights", "stratified", "valid_var"), names(mc), 0)
  dtf <- mc[c(1, m)]
  dtf[[2]] <- data
  dtf[[1]] <- as.name("selectdatavars")
  dtf <- eval(dtf, parent.frame())

  # defining quantile variables
  if (is.null(q)) {Q = as.matrix(dtf[, mix_name]); qi = NULL}
  else {
    q_f = gwqs_rank(dtf, mix_name, q)
    Q = q_f$Q
    qi = q_f$qi
  }
  dtf[,mix_name] <- Q

  # create stratified variables if stratified is not NULL
  m <- match(c("stratified"), names(mc), 0)
  if(m){
    strtfd_out <- stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name <- strtfd_out$mix_name
    dtf <- cbind(dtf, Q)
  }
  else stratified <- NULL

  # parameters estimation and model fitting
  m <- match(c("weights"), names(mc), 0)
  if(m[1]) dtf$wghts <- wghts <- unlist(dtf[, weights, drop = FALSE])
  else  dtf$wghts <- wghts <- rep(1, nrow(dtf))

  mc[[1]] <- as.name("gwqs")
  mc$data <- dtf
  mc$weights <- "wghts"
  mc$rh <- mc$seed <- mc$valid_var <- mc$stratified <- NULL
  mc[which(names(mc)=="q")] <- list(NULL)
  mc$mix_name <- mix_name
  if(is.numeric(rh)) rh <- 1:rh
  else if(is.list(rh)) mc$valid_var <- "group"
  else stop("rh has to be either numeric or a list")

  set.seed(seed)
  gwqslist <- lapply(rh, function(i){
    if(is.list(rh)){
      mc$data$group <- 0
      mc$data$group[i] <- 1
    }
    gwqsout <- eval(mc)
    return(gwqsout)
  })
  dtf$group <- NULL

  coeflist <- lapply(gwqslist, function(i){
    if(family$family == "multinomial"){
      tmp <- i$fit$coefficients[,1]
      names(tmp) <- rownames(i$fit$coefficients)
      tmp
    }
    else if(zero_infl) i$fit$coefficients$count
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

  if(family$family == "multinomial"){
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
  }
  else{
    n_levels <- levelnames <- NULL
    wl <- lapply(gwqslist, function(i) i$final_weights)
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
    names(final_weights) <- c("Estimate", "2.5 %", "97.5%")
  }
  final_weights <- cbind(rownames(final_weights), final_weights)
  names(final_weights)[1] <- "mix_name"

  # estimate wqs index
  if(family$family == "multinomial"){
    wqs <- Q%*%as.matrix(final_weights[match(mix_name, final_weights[,1]), levelnames])
    dtf <- cbind(dtf, wqs)
  }
  else dtf$wqs <- wqs <- as.numeric(Q%*%final_weights[match(mix_name, final_weights[,1]), 2])

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
  if(family$family == "multinomial"){
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
  }
  else{
    final_weights <- final_weights[order(final_weights$Estimate, decreasing = TRUE),]
    y_adj_wqs_df <- as.data.frame(cbind(y_plot, wqs))
    names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", "y", "y_adj"), "wqs")
  }

  if(family$family == "multinomial") dtf$wqs <- NULL

  fit <- list(coefficients = NULL, aic = NULL, dispersion = NULL, deviance = NULL, df.residual = NULL,
              null.deviance = NULL, df.null = NULL, theta = NULL, SE.theta = NULL)

  fit$coefficients <- coefest
  if(family$family == "multinomial"){
    fit$aic <- mean(sapply(gwqslist, function(i) 2*dim(i$fit$coefficients)[1] + 2*(i$fit$nlm_out$minimum)), na.rm = TRUE)
    fit$deviance <- mean(sapply(gwqslist, function(i) 2*i$fit$nlm_out$minimum), na.rm = TRUE)
  }
  else{
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
  }

  results <- list(fit = fit, final_weights = final_weights, wqs = wqs, qi = qi, y_wqs_df = y_adj_wqs_df,
                  family = family, call = cl, formula = formula, mix_name = mix_name, stratified = stratified,
                  q = q, n_levels = n_levels, zero_infl = zero_infl, zilink = zilink, levelnames = levelnames,
                  data = dtf, gwqslist = gwqslist, coefmat = coefmat, wmat = wmat)

  if(zero_infl) results$formula <- ff
  class(results) <- "gwqsrh"

  return(results)

}

#' Methods for gwqs objects
#'
#' Methods for extracting information from fitted Weighted Quantile Sum (WQS) regression model objects
#' of class "gwqs".
#'
#' @param object,x An object of class "gwqs" as returned by \code{gwqs}.
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
#' Both the \code{fitted} and \link[stats]{predict} methods can compute fitted responses. The latter
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
#'                b = 1, b1_pos = TRUE, b1_constr = FALSE, family = gaussian)
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

summary.gwqs <- function(object, ...){
  if(object$family$family == "multinomial"){
    ans <- vector("list", 7)
    ans$call <- object$call
    ans$is.binomial <- FALSE
    ans$digits <- options()$digits
    ans$coefficients <- matrix(object$fit$coefficients$Estimate, nrow = dim(object$fit$coefficients)[1]/2, byrow = T)
    ans$standard.errors <- matrix(object$fit$coefficients$Standard_Error, nrow = dim(object$fit$coefficients)[1]/2, byrow = T)
    rownames(ans$coefficients) <- rownames(ans$standard.errors) <- levels(unlist(object$data[, all.vars(object$formula)[1]]))[-1]
    object$data$wqs <- 0
    colnames(ans$coefficients) <- colnames(ans$standard.errors) <- colnames(model.matrix(object$formula, object$data))
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
  class(ans) <- "summary.gwqs"
  return(ans)
}

#' @rawNamespace S3method(summary, gwqsrh)
#' @rdname methods

summary.gwqsrh <- function(object, sumtype = c("norm", "perc"), ...){
  sumtype <- match.arg(sumtype)
  if(object$zero_infl){
    if(sumtype == "norm"){
      countcoefest <- object$fit$coefficients$countcoefest[,c(1,2,3,4)]
      colnames(countcoefest) <- c("Estimate", "Std. Error", "2.5 %", "97.5 %")
      zerocoefest <- object$fit$coefficients$zerocoefest[,c(1,2,3,4)]
      colnames(zerocoefest) <- c("Estimate", "Std. Error", "2.5 %", "97.5 %")
    }
    if(sumtype == "perc"){
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
    if(sumtype == "perc"){
      coefest <- object$fit$coefficients[,c(5,6,7)]
      colnames(coefest) <- c("Estimate", "2.5 %", "97.5 %")
    }
  }

  ans <- list(coefficients = coefest, call = object$call, family = object$family, zero_infl = object$zero_infl,
              dispersion = object$fit$dispersion, deviance = object$fit$deviance, df.residual = object$fit$df.residual,
              null.deviance = object$fit$null.deviance, df.null = object$fit$df.null, aic = object$fit$aic, theta = object$fit$theta,
              SE.theta = object$fit$SE.theta, zilink = object$zilink)
  ans$is.binomial <- FALSE
  class(ans) <- "summary.gwqsrh"
  return(ans)
}

#' @rawNamespace S3method(print, gwqs)
#' @rdname methods

print.gwqs <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
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

#' @rawNamespace S3method(print, gwqsrh)
#' @rdname methods

print.gwqsrh <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
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

#' @rawNamespace S3method(print, summary.gwqs)
#' @rdname methods

print.summary.gwqs <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
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

#' @rawNamespace S3method(print, summary.gwqsrh)
#' @rdname methods

print.summary.gwqsrh <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
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
      # cat("\n Mean 2 x log-likelihood: ", format(round(x$twologlik, 3), nsmall = dp), "\n")
      invisible(x)
  }
}

#' @rawNamespace S3method(predict, gwqs)
#' @rdname methods

predict.gwqs <- function(object, newdata, type=c("link", "response", "prob", "count", "zero", "class", "probs", "terms"), ...){
  type <- match.arg(type)
  if(missing(newdata)){
    data <- object$data[object$vindex,]
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

#' @rawNamespace S3method(predict, gwqsrh)
#' @rdname methods

predict.gwqsrh <- function(object, newdata, sumtype = c("norm", "perc"), type=c("link", "response", "prob", "count", "zero", "class", "probs", "terms"), ...){
  type <- match.arg(type)
  sumtype <- match.arg(sumtype)
  if(missing(newdata)){
    data <- object$data
    data <- cbind(data, object$wqs)
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

#' @rawNamespace S3method(coef, gwqs)
#' @rdname methods

coef.gwqs <- function(object, ...){
  if(object$family$family == "multinomial"){
    coef <- matrix(object$fit$coefficients$Estimate, nrow = dim(object$fit$coefficients)[1]/2, byrow = T)
    rownames(coef) <- levels(unlist(object$data[, all.vars(object$formula)[1]]))[-1]
    object$data$wqs <- 0
    colnames(coef) <- colnames(model.matrix(object$formula, object$data))
    coef
  }
  else coef(object$fit)
}

#' @rawNamespace S3method(coef, gwqsrh)
#' @rdname methods

coef.gwqsrh <- function(object, sumtype = c("norm", "perc"), ...){
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

#' @rawNamespace S3method(vcov, gwqs)
#' @rdname methods

vcov.gwqs <- function(object, model = c("full", "count", "zero"), ...){
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

#' @rawNamespace S3method(vcov, gwqsrh)
#' @rdname methods

vcov.gwqsrh <- function(object, model = c("full", "count", "zero"), ...){
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

#' @rawNamespace S3method(fitted, gwqs)
#' @rdname methods

fitted.gwqs <- function(object, type = c("prob", "response"), ...){
  type <- match.arg(type)
  if(object$family$family == "multinomial") predict.gwqs(object, type = type)
  else fitted(object$fit)
}

#' @rawNamespace S3method(fitted, gwqsrh)
#' @rdname methods

fitted.gwqsrh <- function(object, sumtype = c("norm", "perc"), type = c("prob", "response"), ...){
  type <- match.arg(type)
  sumtype <- match.arg(sumtype)
  predict.gwqsrh(object, sumtype = sumtype, type = type)
}

#' @rawNamespace S3method(residuals, gwqs)
#' @rdname methods

residuals.gwqs <- function(object, type = c("deviance", "pearson", "working", "response", "partial"), ...){
  type <- match.arg(type)
  if(object$family$family == "multinomial"){
    if(!(type %in% c("pearson", "response"))) stop("If family is \"multinomial\" then residuals type must be \"response\" or \"pearson\"\n")
    else{
      y <- unlist(object$data[object$vindex, as.character(object$formula[[2]])])
      res <- mapply(function(i) i == y, levels(y)) - predict.gwqs(object, type = "prob")
      if(type == "pearson") res <- res/sqrt(predict.gwqs(object, type = "prob"))
      res
    }
  }
  else residuals(object$fit, type)
}

#' @rawNamespace S3method(residuals, gwqsrh)
#' @rdname methods

residuals.gwqsrh <- function(object, sumtype = c("norm", "perc"), type = c("pearson", "response"), ...){
  type <- match.arg(type)
  sumtype <- match.arg(sumtype)
  y <- unlist(object$data[, as.character(object$formula[[2]])])
  if(object$family$family == "multinomial"){
    if(!(type %in% c("pearson", "response"))) stop("If family is \"multinomial\" then residuals type must be \"response\" or \"pearson\"\n")
    else{
      res <- mapply(function(i) i == y, levels(y)) - predict.gwqsrh(object, sumtype = sumtype, type = "prob")
      if(type == "pearson") res <- res/sqrt(predict.gwqsrh(object, sumtype = sumtype, type = "prob"))
    }
  }
  else{
    res <- y - predict.gwqsrh(object, sumtype = sumtype, type = "response")
    if(object$zero_infl & type == "pearson"){
      mu <- predict(object, sumtype = sumtype, type = "count")
      phi <- predict(object, sumtype = sumtype, type = "zero")
      theta1 <- switch(object$family$family, poisson = 0, negbin = 1/object$fit$theta)
      vv <- fitted(object, sumtype = sumtype, type = "response")*(1 + (phi + theta1) * mu)
      res <- res/sqrt(vv)
    }
    else if(type == "pearson") res <- res/sqrt(object$family$variance(predict.gwqsrh(object, sumtype = sumtype, type = "response")))
  }
  res
}

#' Plots and tables functions
#'
#' Functions that allow to generate plots and tables helping in visualizing and summarise Weighted Quantile Sum (WQS) regression results.
#'
#' @param object An object of class "gwqs" as returned by \code{gwqs}.
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
#' @param ... Further arguments to be passed to the function.
#'
#' @details
#' The \code{gwqs_barplot}, \code{gwqs_scatterplot}, \code{gwqs_fitted_vs_resid}, \code{gwqs_levels_scatterplot},
#' \code{gwqs_ROC} and \code{gwqsrh_boxplot} functions produce five figures through the \code{\link[ggplot2]{ggplot}} function.
#'
#' The \code{gwqs_summary_tab} and \code{gwqs_weights_tab} functions produce two tables in the viewr pane
#' through the use of the \code{\link[knitr]{kable}} and \code{\link[kableExtra]{kable_styling}} functions.
#'
#' The \code{gwqs_barplot}, \code{gwqs_scatterplot} plots are available for all family types while
#' \code{gwqs_fitted_vs_resid} is not available when \code{family = binomial} or \code{"multinomial"}.
#' \code{gwqs_levels_scatterplot} plot is only available when \code{family = "multinomial"} and \code{gwqs_ROC}
#' when \code{family = binomial}. All these plots can also be applied to the objects of class \code{gwqsrh}.
#' For these objects an additional plot is available through the function \code{gwqs_boxplot}.
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
#'                b = 1, b1_pos = TRUE, b1_constr = FALSE, family = gaussian)
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
#' @importFrom reshape2 melt
#'
#' @export gwqs_barplot
#' @rdname secondary_functions

# Functions to generate the plots
gwqs_barplot <- function(object, tau, ...){
  if(object$family$family == "multinomial"){
    if(class(object) == "gwqsrh") object$final_weights <- object$final_weights[,c(1, 3*0:(object$n_levels-2)+2)]
    data_plot <- object$final_weights[order(object$final_weights[, object$levelnames[1]]),]
    pos <- match(data_plot$mix_name, sort(object$mix_name))
    data_plot$mix_name <- factor(data_plot$mix_name, levels(data_plot$mix_name)[pos])
    data_plot_l <- melt(data_plot, id.vars = "mix_name")
    bar_plot_h <- ggplot(data_plot_l, aes_string(x = "mix_name", y = "value")) +
      facet_wrap(~ variable)
  }
  else {
    if(class(object) == "gwqsrh"){
      object$final_weights <- object$final_weights[,c(1, 2)]
      names(object$final_weights) <- c("mix_name", "mean_weight")
    }
    data_plot <- object$final_weights[order(object$final_weights$mean_weight),]
    data_plot$mix_name <- factor(data_plot$mix_name, levels = data_plot$mix_name)
    bar_plot_h <- ggplot(data_plot, aes_string(x = "mix_name", y = "mean_weight"))
  }

  bar_plot_h <- bar_plot_h + geom_bar(stat = "identity", color = "black") + theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(color='black'),
          legend.position = "none") + coord_flip()
  if(missing(tau)) tau <- 1/length(object$mix_name)
  if(!is.null(tau)) bar_plot_h <- bar_plot_h + geom_hline(yintercept = tau, linetype="dashed", color = "red")

  print(bar_plot_h)
}

#' @export gwqs_scatterplot
#' @rdname secondary_functions

gwqs_scatterplot <- function(object, ...){
  y_labs = ifelse(object$family$family %in% c("multinomial", "binomial"), "y", "y_adj")

  yadj_vs_wqs = ggplot(object$y_wqs_df, aes_string("wqs", y_labs)) +
    geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()

  if(object$family$family == "multinomial") yadj_vs_wqs = yadj_vs_wqs + facet_wrap(~ level)

  print(yadj_vs_wqs)
}

#' @export gwqs_fitted_vs_resid
#' @rdname secondary_functions

gwqs_fitted_vs_resid <- function(object, sumtype = c("norm", "perc"), ...){
  sumtype <- match.arg(sumtype)
  if(!(object$family$family %in% c("binomial", "multinomial"))){
    if(class(object) == "gwqs") fit_df = data.frame(fitted = fitted(object), resid = residuals(object, type = "response"))
    else if(class(object) == "gwqsrh") fit_df = data.frame(fitted = fitted(object, sumtype, type = "response"), resid = residuals(object, sumtype = sumtype, type = "response"))
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
    if(class(object) == "gwqsrh") object$final_weights <- object$final_weights[,c(1, 3*0:(object$n_levels-2)+2)]
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
    if(class(object) == "gwqs") wqs_pred <- predict(object, newdata, type = "response")
    else if(class(object) == "gwqsrh") wqs_pred <- predict(object, newdata, sumtype, type = "response")
    df_roc <- wqs_pred$df_pred
    if(class(df_roc$y) == "character") df_roc$y = factor(df_roc$y)
    if(class(df_roc$y) == "factor") df_roc$y <- as.numeric(df_roc$y != levels(df_roc$y)[1])
    gg_roc = suppressWarnings(ggplot(df_roc, aes_string(d="y", m="ypred")) + geom_roc(n.cuts = 0) +
                                style_roc(xlab = "1 - Specificity", ylab = "Sensitivity"))
    auc_est = calc_auc(gg_roc)
    gg_roc = gg_roc + annotate("text", x=0.75, y=0.25, label=paste0("AUC = ", round(auc_est[, "AUC"], 3)))

    print(gg_roc)
  }
  else stop("The ROC curve is only available for binomial family\n")
}

#' @export gwqsrh_boxplot
#' @rdname secondary_functions

gwqsrh_boxplot <- function(object, tau, ...){
  if(class(object) == "gwqsrh"){
    if(object$family$family == "multinomial"){
      wboxplotl <- lapply(object$levelnames, function(i){
        tmp <- melt(object$wmat[[i]], varnames = c("rh", "mix_name"))
        tmp$level <- i
        return(tmp)
      })
      wboxplot <- do.call("rbind", wboxplotl)
    }
    else wboxplot <- melt(object$wmat, varnames = c("rh", "mix_name"))
    wboxplot$mix_name <- factor(wboxplot$mix_name, levels = object$final_weights$mix_name)
    box_plot <- ggplot(wboxplot, aes_string(x = "mix_name", y = "value")) +
      geom_boxplot(outlier.shape = " ") + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Weight (%)") +
      stat_summary(fun.y = mean, geom = "point", shape = 18, size = 3) + geom_jitter(alpha = 0.3)
    if(object$family$family == "multinomial") box_plot <- box_plot + facet_wrap(~level)
    if(missing(tau)) tau <- 1/length(object$mix_name)
    if(!is.null(tau)) box_plot <- box_plot + geom_hline(yintercept = tau, linetype="dashed", color = "red")

    print(box_plot)
  }
  else stop("The function gwqsrh_boxplot is only available for objects of class gwqsrh\n")
}

#' @export gwqs_summary_tab
#' @rdname secondary_functions

# functions to generate tables
gwqs_summary_tab <- function(object, sumtype = c("norm", "perc"), ...){
  sumtype <- match.arg(sumtype)
  if(class(object) == "gwqs"){
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
  else if(class(object) == "gwqsrh"){
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


