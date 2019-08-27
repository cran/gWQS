#' Fitting Weighted Quantile Sum regression models
#'
#' Fits Weighted Quantile Sum (WQS) regressions for continuous, binomial, multinomial, poisson, quasi-poisson and negative binomial outcomes.
#'
#' @param formula An object of class \code{formula} specifying the relationship to be tested. The \code{wqs}
#' term must be included in \code{formula}, e.g. \code{y ~ wqs + ...}. To test for an interaction term with
#' a continuous variable \code{a} or for a quadratic term we can specify the \code{formula} as below:
#' \code{y ~ wqs*a + ...} and \code{y ~ wqs + I(wqs^2) + ...}, respectively.
#' @param data The \code{data.frame} containing the variables to be included in the model.
#' @param na.action \code{\link[stats]{model.frame}}. \code{na.omit} is the default.
#' @param weights an optional vector of weights to be used in the fitting process.
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
#' @param zilink character specification of link function in the binary zero-inflation model
#' (you can choose among "logit", "probit", "cloglog", "cauchit", "log").
#' @param seed An \code{integer} value to fix the seed, if it is equal to \code{NULL} no seed is chosen.
#' @param pred Percentage of the dataset to be used for the predictive model. If \code{pred = 0} then no
#' predicitve model is going to be built.
#' @param plots A logical value indicating whether plots should be generated with the output
#' (\code{plots = TRUE}) or not (\code{plots = FALSE}).
#' @param tables A logical value indicating whether tables should be generated in the output
#' (\code{tables = TRUE}) or not (\code{tables = FALSE}).
#' @param plan_strategy A character value that allows to choose the evaluation strategies for the
#' \code{plan} function. You can choose among "sequential", "transparent", "multisession", "multicore",
#' "multiprocess", "cluster" and "remote" (see \code{\link[future]{plan}} help page for more details).
#' @param control The control list of optimization parameters. See \code{\link[Rsolnp]{solnp}} for details.
#' @param lp The lambda parameter that add a penlization term when we want to constrain in the negative direction.
#' This is an alternative to \code{b1_constr = TRUE}.
#' @param ln The lambda parameter that add a penlization term when we want to constrain in the positive direction.
#' This is an alternative to \code{b1_constr = TRUE}.
#'
#' @details
#' \code{gWQS} uses the \code{glm} function in the \bold{stats} package to fit the linear, logistic,
#' the Poisson and the quasi-Poisson regression, while the \code{glm.nb} function from the \bold{MASS}
#' package is used to fit the negative binomial regression respectively. The \code{nlm} function from
#' the \bold{stats} package was used to optimize the log-likelihood of the multinomial regression.\cr
#'
#' The \code{\link[Rsolnp]{solnp}} optimization function is used to estimate the weights at each
#' bootstrap step.\cr
#'
#' The \code{seed} argument specifies a fixed seed through the \code{\link[base]{set.seed}} function.\cr
#'
#' The \code{plots} argument produces three figures (two if \code{family = binomial} or \code{"multinomial"})
#' through the \code{\link[ggplot2]{ggplot}} function. One more plot will be printed if \code{pred > 0} and
#' \code{family = binomial}.\cr
#'
#' The \code{tables} argument produces two tables in the viewr pane through the use of the functions
#' \code{\link[knitr]{kable}} and \code{\link[kableExtra]{kable_styling}}.\cr
#'
#' @return \code{gwqs} return the results of the WQS regression as well as many other objects and datasets.
#'
#' \item{fit}{The object that summarizes the output of the WQS model, reflecting a
#' linear, logistic, multinomial, Poisson, quasi-Poisson or negative binomial regression
#' depending on how the \code{family} parameter was specified.
#' The summary function can be used to call and print fit data (not for multinomial regression).}
#' \item{conv}{Indicates whether the solver has converged (0) or not (1 or 2).}
#' \item{bres}{Matrix of estimated weights, mixture effect parameter estimates and the associated
#' standard errors, statistics and p-values estimated for each bootstrap iteration.}
#' \item{wqs}{Vector containing the wqs index for each subject.}
#' \item{q_i}{List of the cutoffs used to divide in quantiles the variables in the mixture}
#' \item{bindex}{List of vectors containing the \code{rownames} of the subjects included in each
#' bootstrap dataset.}
#' \item{tindex}{Vector containing the rows used to estimate the weights in each bootstrap.}
#' \item{vindex}{Vector containing the rows used to estimate the parameters of the final model.}
#' \item{final_weights}{\code{data.frame} containing the final weights associated to each chemical.}
#' \item{y_wqs_df}{\code{data.frame} containing the dependent variable values adjusted for the
#' residuals of a fitted model adjusted for covariates (original values when \code{family = binomial}
#' or \code{"multinomial"}) and the wqs index estimated values.}
#' \item{df_pred}{\code{data.frame} containing the variables to print the ROC curve. It is generated only
#' when \code{pred > 0}}
#' \item{pindex}{Vector containing the subjects used for prediction. It is generated only when \code{pred > 0}}
#'
#' @author
#' Stefano Renzetti, Paul Curtin, Allan C Just, Ghalib Bello, Chris Gennings
#'
#' @references
#' Renzetti S, Gennings C, Curtin PC. 2019. gWQS: An R Package for Linear and Generalized Weighted
#' Quantile Sum (WQS) Regression. Journal of Statistical Software.\cr
#'
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
#' @examples
#' # we save the names of the mixture variables in the variable "toxic_chems"
#' toxic_chems = c("log_LBX074LA", "log_LBX099LA", "log_LBX105LA", "log_LBX118LA",
#' "log_LBX138LA", "log_LBX153LA", "log_LBX156LA", "log_LBX157LA", "log_LBX167LA",
#' "log_LBX170LA", "log_LBX180LA", "log_LBX187LA", "log_LBX189LA", "log_LBX194LA",
#' "log_LBX196LA", "log_LBX199LA", "log_LBXD01LA", "log_LBXD02LA", "log_LBXD03LA",
#' "log_LBXD04LA", "log_LBXD05LA", "log_LBXD07LA", "log_LBXF01LA", "log_LBXF02LA",
#' "log_LBXF03LA", "log_LBXF04LA", "log_LBXF05LA", "log_LBXF06LA", "log_LBXF07LA",
#' "log_LBXF08LA", "log_LBXF09LA", "log_LBXPCBLA", "log_LBXTCDLA", "log_LBXHXCLA")
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
#'                b = 2, b1_pos = TRUE, b1_constr = FALSE, family = gaussian, seed = 2016,
#'                plots = TRUE, tables = TRUE)
#'
#' # to test the significance of the covariates
#' summary(results$fit)
#'
#' @import ggplot2
#' @import Rsolnp
#' @import stats
#' @import knitr
#' @import kableExtra
#'
#' @importFrom broom augment
#' @importFrom rlist list.cbind
#' @importFrom MASS glm.nb
#' @importFrom reshape2 melt
#' @importFrom plotROC geom_roc style_roc calc_auc
#' @importFrom nnet multinom
#' @importFrom future plan future value
#' @importFrom future.apply future_lapply
#' @importFrom ggrepel geom_text_repel
#' @importFrom pscl zeroinfl
#'
#' @export

gwqs <- function(formula, data, na.action, weights, mix_name, stratified, valid_var, b = 100,
                 b1_pos = TRUE, b1_constr = FALSE, zero_infl = FALSE, q = 4, validation = 0.6,
                 family = gaussian, zilink = c("logit", "probit", "cloglog", "cauchit", "log"),
                 seed = NULL, pred = 0, plots = FALSE, tables = FALSE, plan_strategy = "sequential",
                 control = list(rho = 1, outer.iter = 400, inner.iter = 800, delta = 1.0e-7, tol = 1e-8, trace = 0),
                 lp = 0, ln = 0){

  if(is.character(family)){
    if(family %in% c("multinomial", "negbin")) family <- list(family = family)
    else family <- get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if(!any(grepl("wqs", rownames(attr(terms(formula), "factors"))))) stop("'formula' must contain 'wqs' term: e.g. y ~ wqs + ...")

  if(zero_infl){
    zilink <- make.link(match.arg(zilink))
    ff = formula
    if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
      formula = as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))
  }
  else zilink <- ff <- NULL

  # cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", "weights", "stratified", "valid_var"), names(mc), 0)
  dtf <- mc[c(1, m)]
  dtf[[2]] <- data
  dtf[[1]] <- as.name("select_vars")
  dtf <- eval(dtf, environment(gwqs))

  # select rows of the original dataset based on rows selected with na.action

  # defining quantile variables
  if (is.null(q)) {Q = as.matrix(dtf[, mix_name]); q_i = NULL}
  else {
    q_f = quantile_f(dtf, mix_name, q)
    Q = q_f$Q
    q_i = q_f$q_i
  }

  # create stratified variables if stratified is not NULL
  m <- match(c("stratified"), names(mc), 0)
  if (m){
    if(is.null(q)) stop("q must be different from NULL if you want to stratify for a categorical variable")
    strtfd_out = stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name = strtfd_out$mix_name
  }

  if(!is.null(seed) & !is.numeric(seed)) stop("seed must be numeric or NULL")
  if(!is.null(seed)) set.seed(seed)

  N <- nrow(dtf)

  # splitting the dataset
  m <- match(c("valid_var"), names(mc), 0)
  rindex = create_rindex(dtf, N, validation, pred, valid_var, m, family)

  # parameters estimation and model fitting
  m <- match(c("weights"), names(mc), 0)
  if(m[1]) weights <- dtf[, weights, drop = FALSE]
  else  dtf$weights <- weights <- rep(1, N)
  par_model <- par.modl.est(dtf[rindex$it,], Q[rindex$it,], formula, ff, weights[rindex$it], b, b1_pos, b1_constr, family, zilink, zero_infl, plan_strategy, control, lp, ln)

  bres = par_model$bres
  conv = par_model$conv
  bindex = par_model$bindex
  nfuneval = par_model$n_funeval
  n_levels = par_model$n_levels
  strata_names = par_model$strata_names

  mean_weight = mean_weight_f(mix_name, bres, conv, b1_pos, family, n_levels, strata_names, zero_infl)

  # fit the final model with the estimated weights
  wqs_model = model.fit(mean_weight, dtf[rindex$iv,], Q[rindex$iv,], family, zilink, formula, ff, weights[rindex$iv], b1_pos, zero_infl)

  if(all(grepl("wqs", attr(terms(formula), "term.labels")))) y_plot <- model.response(model.frame(formula, dtf[rindex$iv,]), "any")
  else{
    formula_wo_wqs = remove_terms(formula, "wqs")
    if(family$family != "multinomial"){
      if(zero_infl){
        if(length(ff[[3]]) > 1 && identical(ff[[3]][[1]], as.name("|"))){
          f1 <- remove_terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]]))), "wqs")
          f2 <- remove_terms(as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[3]]))), "wqs")
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

  # prediction
  if(!is.null(rindex$ip)) df_pred = predict_f(Q[rindex$ip,], dtf[rindex$ip,], mean_weight, wqs_model$m_f, formula)
  else df_pred = NULL

  # Plots
  data_plot <- data.frame(mix_name, mean_weight)

  if(family$family == "multinomial"){
    Y <- model.response(model.frame(formula, dtf[rindex$iv,]), "any")
    level_names = levels(Y)
    names(data_plot)[2:n_levels] = strata_names
    Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
    colnames(Y) = strata_names
    y_adj_wqs_df = do.call("rbind", lapply(strata_names, function(i){
      if(n_levels > 3) data.frame(level = i,
                                  y = Y[rowSums(Y[, -which(colnames(Y)==i)]) == 0, i],
                                  wqs = wqs_model$wqs[rowSums(Y[, -which(colnames(wqs_model$wqs)==i)]) == 0, i])
      else data.frame(level = i,
                      y = Y[Y[, -which(colnames(Y)==i)] == 0, i],
                      wqs = wqs_model$wqs[Y[, -which(colnames(wqs_model$wqs)==i)] == 0, i])
    }))
  }
  else{
    y_adj_wqs_df = as.data.frame(cbind(y_plot, wqs_model$wqs))
    names(y_adj_wqs_df) = c(ifelse(family$family == "binomial", "y", "y_adj"), "wqs")
  }

  if(plots) plots(data_plot, y_adj_wqs_df, q, mix_name, mean_weight, wqs_model$m_f, family,
                          n_levels, strata_names, df_pred, zero_infl)

  if(family$family == "multinomial"){
    data_plot = data_plot[order(data_plot[, strata_names[1]], decreasing = TRUE),]
    wqs_index = wqs_model$wqs
  }
  else{
    data_plot = data_plot[order(data_plot$mean_weight, decreasing = TRUE),]
    wqs_index = as.numeric(unlist(wqs_model$wqs))
  }

  # Tables
  if (tables) tables(data_plot, wqs_model$m_f, family, n_levels, zero_infl)

  # creating the list of elements to return
  results = list(wqs_model$m_f, conv, bres, wqs_index, q_i, bindex, rindex$it, rindex$iv,
                 data_plot, y_adj_wqs_df, df_pred, rindex$ip)
  names(results) = c("fit", "conv", "bres", "wqs", "q_i", "bindex", "tindex", "vindex",
                     "final_weights", "y_wqs_df", "df_pred", "pindex")

  return(results)
}

