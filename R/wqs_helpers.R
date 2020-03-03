# starting message
.onAttach <- function(...) packageStartupMessage("Welcome to Weighted Quantile Sum (WQS) Regression.
If you are using a Mac you have to install XQuartz.
You can download it from: https://www.xquartz.org/\n")


# function to remove terms from formula
remove_terms <- function(form, term) {
  fterms <- terms(form)
  fac <- attr(fterms, "factors")
  if(term %in% rownames(fac)){
    fac_wqs <- fac[grep(term, rownames(fac)), ]
    if(NCOL(fac_wqs) == 1) idx <- which(as.logical(fac[term, ]))
    else idx <- which(apply(fac_wqs, 2, function(i) any(i==1)))
    new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
    return(formula(new_fterms))
  }
  else return(form)
}


# function to select variables
select_vars <- function(data, na.action, formula, mix_name, ...){
  allvars = all.vars(formula)
  other_vars <- c(...)
  data$wqs = 0
  data <- data[, c(allvars, mix_name, other_vars)]
  if(missing(na.action)) na.action <- na.omit
  dtf <- na_action(data, na.action)
  return(dtf)
}


# function to manage NAs
na_action <- function(data, na.action){
  dtf <- match.call(expand.dots = FALSE)
  m <- match(c("na.action", "data"), names(dtf), 0)
  dtf <- dtf[m]
  names(dtf)[2] <- "object"
  dtf <- eval(dtf, parent.frame())
  return(dtf)
}


# function to create stratified elements in the mixture by levels of factors
stratified_f = function(Q, dtf, stratified, mix_name){
  ls = levels(unlist(dtf[, stratified, drop = FALSE]))
  if(is.null(ls)) stop("'stratified' must be factor")
  llsm = lapply(ls, function(ls){
    mat = diag(as.numeric(dtf[, stratified] == ls))
    sub_dt = mat%*%Q[, mix_name]
    colnames(sub_dt) = paste(mix_name, ls, sep = "_")
    return(sub_dt)
  })
  Q = do.call("cbind", llsm)
  mix_name = colnames(Q)
  strtfd_out = list(Q, mix_name)
  names(strtfd_out) = c("Q", "mix_name")

  return(strtfd_out)
}


# function to create variables with quantile of the components
quantile_f <- function(dtf, mix_name, q){

  if(!is.numeric(q)) stop("'q' must be a number")
  Ql <- lapply(1:length(mix_name), function(i){
    q_i <- unique(quantile(dtf[[mix_name[i]]], probs = seq(0, 1, by = 1/q), na.rm = TRUE))
    if(length(q_i) == 1) q_i = c(-Inf, q_i)
    q <- cut(dtf[[mix_name[i]]], breaks = q_i, labels = FALSE, include.lowest = TRUE) - 1
    return(list(q_i, q))
  })
  q_i <- lapply(Ql, function(x) x[[1]])
  Q <- matrix(unlist(lapply(Ql, function(x) x[[2]])), ncol = length(Ql))
  colnames(Q) <- names(q_i) <- mix_name

  qf_out <- list(Q, q_i)
  names(qf_out) <- c("Q", "q_i")
  return(qf_out)
}


# function to split the dataset
create_rindex <- function(dtf, N, validation, pred, valid_var, m, family){

  if(!m){
    if(!is.numeric(validation) | validation < 0 | validation >= 1) stop("'validation' must be numeric >= 0 and < 1")
    if(!is.numeric(pred) | pred < 0 | pred >= 1) stop("'pred' must be numeric >= 0 and < 1")
    if(pred + validation >= 1) stop("the sum of pred and validation must be between 0 and 1")
    if(pred>0 & family$family == "multinomial") stop("The predictive model is not available for multinomial regression")
    groups <- rep(0, N)
    if(validation > 0) groups[sample(1:N, round(N*validation))] <-1
    if(pred > 0) groups[sample(which(groups!=1), round(N*pred))] <-2
  }
  else{
    groups = unlist(dtf[, valid_var, drop = FALSE])
    if(!any(unique(groups) %in% c(0, 1, 2))) stop("valid_var values must be 0, 1 or 2")
    if(!(0 %in% unique(groups))) stop(("0 must identify test dataset"))
  }
  it = which(groups == 0)
  iv = which(groups == 1)
  ip = which(groups == 2)
  if(length(iv) == 0) iv = it
  if(length(ip) == 0) ip = NULL

  indexl = list(it, iv, ip)
  names(indexl) = c("it", "iv", "ip")
  return(indexl)
}


# parameter names in model matrix
parnames <- function(df, formula, form2){
  if(!is.null(form2)){
    mf <- model.frame(form2, df)
    Y <- model.response(mf, "any")
    df$yz <- ifelse(Y == 0, 0, 1)
  }
  mm <- model.matrix(formula, df)
  colnames(mm)
}


# functtion to define the objective function
objfn <- function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, weights, b1_constr, lp, ln){

  mf <- model.frame(formula, bdtf)
  Y <- model.response(mf, "any")
  if(family$family == "binomial" & class(Y) %in% c("factor", "character")){
    if(class(Y) == "character") Y = factor(Y)
    Y <- as.numeric(Y != levels(Y)[1])
  }
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0
  if (family$family == "multinomial"){
    w <- matrix(initp[(kx + 1):length(initp)], kw, n_levels-1)
    bdtf[, wqsvars] <- Q%*%w
    Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
    fm_l = sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
    Xl = lapply(fm_l, function(i) model.matrix(i, data = bdtf))
    X = do.call("cbind", Xl)
    b_covs = matrix(0, kx, n_levels-1)
    i = 1:(kx/(n_levels-1))
    for (j in 0:(n_levels-2)){
      b_covs[(kx/(n_levels-1))*j+i, j+1] = initp[(kx/(n_levels-1))*j+i]
    }
    term = X%*%b_covs + offset
  }
  else{
    w <- initp[(kx + 1):(kx + kw)]
    bdtf$wqs <- as.numeric(Q%*%w)
    X <- model.matrix(formula, bdtf)
    b_covs <- initp[1:kx]
    if(family$family == "negbin") theta <- exp(initp[length(initp)])
    term <- as.numeric(X%*%b_covs) + offset
  }

  if(family$family == "multinomial") f = -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*weights)
  else if(family$family == "negbin") f = -sum((suppressWarnings(dnbinom(Y, size = theta, mu = exp(term), log = TRUE)))*weights)
  else if(family$family %in% c("poisson", "quasipoisson")) f = -sum(dpois(Y, lambda = exp(term), log = TRUE))
  else f = sum(family$dev.resids(y = Y, mu = family$linkinv(term), wt = weights)) + lp*max(0, initp["wqs"]) - ln*min(0, initp["wqs"])

  return(f)
}


# function to define the equality constraint
linconst = function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, weights, b1_constr, lp, ln){

  if (family$family == "multinomial"){
    w = matrix(initp[(kx + 1):length(initp)], kw, (n_levels-1))
    wsum = colSums(w)
  }
  else{
    w = initp[(kx + 1):(kx + kw)]
    wsum = sum(w)
  }

  return(wsum)
}


# function to determine bounded prameters
bounded_param = function(initp, kw, bdtf, Q, kx, Xnames, n_levels, level_names, wqsvars, family, zilink, zero_infl, formula, ff, weights, b1_constr, lp, ln){

  bp = initp[(kx + 1):(kx + kw*(n_levels-1))]
  if (b1_constr){
    wqs_site = which(Xnames == "wqs")
    bp = c(initp[wqs_site], bp)
  }

  return(bp)
}


# function to define the lower bounds
LBound = function(kw, n_levels, family, b1_pos, b1_constr){

  LB = rep(0, kw)
  if(family$family == "multinomial") LB = rep(LB, n_levels-1)
  if(b1_constr) LB = c(sapply(b1_pos, function(i) ifelse(i, 0, -Inf)), LB)

  return(LB)
}


# function to define the upper bounds
UBound = function(kw, n_levels, family, b1_pos, b1_constr){

  UB = rep(1, kw)
  if(family$family == "multinomial") UB = rep(UB, n_levels-1)
  if(b1_constr) UB = c(sapply(b1_pos, function(i) ifelse(i, Inf, 0)), UB)

  return(UB)
}


# function to define the parameters initial values
values.0 = function(kw, Xnames, kx, n_levels, formula, ff, weights, bdtf, b1_pos, family, zilink, zero_infl){

  w = rep(1/kw, kw)

  if(family$family == "multinomial"){
    fit = multinom(formula, bdtf, trace = F, weights = weights)
    bj = c(sapply(1:(n_levels-1), function(i) coef(fit)[i,]))
    w = rep(w, (n_levels-1))
    val.0 = c(bj, w)
  }
  else{
    if(family$family == "negbin") fit = suppressWarnings(glm.nb(formula, bdtf, weights = weights))
    else fit = glm(formula, bdtf, family = family, weights = weights)
    bj = coef(fit)
    val.0 = c(bj, w)
    if(family$family == "negbin"){
      if(length(attr(terms(formula), "term.labels")) == 1) val.0 = c(val.0, 1)
      else val.0 = c(val.0, log(fit$theta))
    }
  }

  wqs_site <- which(grepl("wqs", Xnames))
  val.0[wqs_site] = sapply(b1_pos, function(i) ifelse(i, 0.0001, -0.0001))

  return(val.0)
}


# optimization function to estimate the weights
optim.f <- function(bdtf, bQ, b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, weights, control, lp, ln){

  Xnames <- parnames(bdtf, formula, NULL)
  kx <- length(Xnames)
  if(family$family == "multinomial"){
    n_levels <- nlevels(eval(formula[[2]], envir = bdtf))
    if(n_levels == 0) stop("When 'family = \"multinomial\"' y must be of class factor")
    level_names <- levels(eval(formula[[2]], envir = bdtf))
    Xnames <- c(sapply(level_names[-1], function(i) paste0(Xnames[1:kx], "_", i, "_vs_", level_names[1])))
    kx <- kx*(n_levels-1)
    wqsvars = Xnames[grepl("^wqs_", Xnames)]
    bdtf[, wqsvars] <- 0
  }
  else {n_levels <- 2; level_names <- wqsvars <- NULL}
  kw <- dim(bQ)[2]
  initp <- values.0(kw, Xnames, kx, n_levels, formula, ff, weights, bdtf, b1_pos, family, zilink, zero_infl)
  LowB <- LBound(kw, n_levels, family, b1_pos, b1_constr)
  UpB <- UBound(kw, n_levels, family, b1_pos, b1_constr)
  eq_b <- rep(1, length(b1_pos))

  if(control$trace) cat("start opt")

  opt_res <- tryCatch(solnp(pars = initp, fun = objfn, eqfun = linconst, eqB = eq_b, ineqfun = bounded_param,
                            ineqLB = LowB, ineqUB = UpB, control = control, kw = kw, bdtf = bdtf, Q = bQ, kx = kx,
                            Xnames = Xnames, n_levels = n_levels, level_names = level_names,
                            wqsvars = wqsvars, family = family, zilink = zilink, zero_infl = zero_infl,
                            formula = formula, ff = ff, weights = weights, b1_constr = b1_constr, lp = lp, ln = ln),
                      error = function(e) NULL)

  if(!is.null(opt_res)) {
    if(family$family == "multinomial") par_opt <- matrix(opt_res$pars[(kx + 1):length(initp)], kw, n_levels-1)
    else par_opt <- opt_res$pars[(kx + 1):(kx + kw)]
    conv <- opt_res$convergence
    nfuneval <- opt_res$nfuneval
  }
  else {
    if(family$family == "multinomial") par_opt <- matrix(initp[(kx + 1):length(initp)], kw, n_levels-1)
    else par_opt <- initp[(kx + 1):(kx + kw)]
    conv <- 2
    nfuneval <- 0
  }
  out <- list(par_opt, conv, nfuneval)
  names(out) <- c("par_opt", "conv", "nfuneval")

  return(out)
}


# function that fit the wqs model
model.fit <- function(w, bdtf, bQ, family, zilink, formula, ff, weights, b1_pos, zero_infl){

  if (family$family == "multinomial"){
    mf <- model.frame(formula, bdtf)
    Y <- model.response(mf, "any")
    n_levels <- nlevels(Y)
    level_names <- levels(Y)
    Y <- cbind(sapply(2:n_levels, function(i) ifelse(Y == level_names[i], 1, 0)))
    offset <- model.offset(mf)
    if(is.null(offset)) offset <- 0
    Xnames <- parnames(bdtf, formula, NULL)
    Xnames <- c(sapply(level_names[-1], function(i) paste0(Xnames[1:length(Xnames)], "_", i, "_vs_", level_names[1])))
    kx <- length(Xnames)
    kw <- dim(bQ)[2]
    w <- matrix(w, kw, n_levels-1)
    wqsvars <- Xnames[grepl("^wqs_", Xnames)]
    wqs <- bdtf[, wqsvars] <- bQ%*%w
    colnames(wqs) <- paste0(level_names[-1], "_vs_", level_names[1])

    fm_l <- sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
    Xl <- lapply(fm_l, function(i) model.matrix(i, data = bdtf))
    X <- do.call("cbind", Xl)
    initp <- values.0(kw, Xnames, kx, n_levels, formula, ff, weights, bdtf, b1_pos, family, zilink, zero_infl)
    initp <- initp[1:kx]
  }
  else wqs <- bdtf$wqs <- as.numeric(bQ%*%w)

  if(zero_infl) m_f <- pscl::zeroinfl(ff, bdtf, dist = family$family, link = zilink$name)
  else{
    if(family$family == "multinomial") {
      LL <- function(p){
        b_covs = matrix(0, kx, n_levels-1)
        i = 1:(kx/(n_levels-1))
        for (j in 0:(n_levels-2)){
          b_covs[(kx/(n_levels-1))*j+i, j+1] = p[(kx/(n_levels-1))*j+i]
        }
        term = X%*%b_covs + offset
        -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*weights)
      }
      nlm_out = nlm(LL, initp, hessian = TRUE)

      Estimate = nlm_out$estimate
      Standard_Error = sqrt(diag(solve(nlm_out$hessian)))
      stat = Estimate/Standard_Error
      p_value = 2*pnorm(-abs(stat))
      sum_stat <- data.frame(Estimate = Estimate, Standard_Error = Standard_Error, stat = stat, p_value = p_value, stringsAsFactors = TRUE)
      rownames(sum_stat) = Xnames
      m_f = list(nlm_out, sum_stat)
      names(m_f) = c("nlm_out", "sum_stat")
    }
    else if(family$family == "negbin") m_f = glm.nb(formula, data = bdtf)
    else m_f = glm(formula, data = bdtf, family = family)
  }

  mf_out = list(wqs, m_f)
  names(mf_out) = c("wqs", "m_f")

  return(mf_out)
}

# function to sample from data rownames
sample_f = function(i, dtf){

  bindex = sample(1:nrow(dtf), nrow(dtf), replace=TRUE)

  return(bindex)
}

# function that call optim.f to estimate parameters for each bootstrap sample
estimate_param = function(bindex, dtf, Q, b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, weights, control, lp, ln){

  param = optim.f(dtf[bindex,], Q[bindex,], b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, weights[bindex], control, lp, ln)

  return(param)
}


# function to be passed to future_lapply to fit the models
model.fit_f = function(i, param, bindex, dtf, Q, b1_pos, family, zilink, formula, ff, weights, zero_infl){

  b_fit <- model.fit(param[[i]]$par_opt, dtf[bindex[[i]],], Q[bindex[[i]],], family, zilink, formula, ff, weights, b1_pos, zero_infl)

  return(b_fit)
}

# function that calls the optimization function and the function to fit the model for each bootstrap sample
par.modl.est <- function(dtf, Q, formula, ff, weights, b, b1_pos, b1_constr, family, zilink, zero_infl, plan_strategy, control, lp, ln){

  if (family$family %in% c("gaussian", "quasipoisson")) ts = "t"
  else if (family$family %in% c("binomial", "poisson", "multinomial", "negbin")) ts = "z"

  if(!is.numeric(b)) stop("'b' must be a number")

  bindex = lapply(X=1:b, FUN = sample_f, dtf = dtf)

  plan(plan_strategy)
  param <- future_lapply(X = bindex, FUN = estimate_param, dtf = dtf, Q = Q, b1_pos = b1_pos, b1_constr = b1_constr,
                         family = family, zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff,
                         weights = weights, control = control, lp = lp, ln = ln, future.seed = FALSE)

  plan(plan_strategy)
  b_fit <- future_lapply(X = 1:b, FUN = model.fit_f, param = param, bindex = bindex, dtf = dtf, Q = Q, b1_pos = b1_pos,
                         family = family, zilink = zilink, formula = formula, ff = ff, weights = weights,
                         zero_infl = zero_infl, future.seed = FALSE)

  conv <- c(sapply(param, function(i) i$conv))
  nfuneval <- c(sapply(param, function(i) i$nfuneval))
  if(family$family == "multinomial"){
    n_levels <- dim(param[[1]]$par_opt)[2]+1
    wqs_site <- which(grepl("^wqs_", rownames(b_fit[[1]]$m_f$sum_stat)))
    wght_matrix <- lapply(1:(n_levels-1), function(j) do.call("rbind", lapply(param, function(i) i$par_opt[,j])))
    b1 <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$nlm_out$estimate[j]))
    se <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$sum_stat$Standard_Error[j]))
    stat <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$sum_stat$stat[j]))
    p_val <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$sum_stat$p_value[j]))
  }
  else{
    wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt))
    if(zero_infl){
      b1_count <- sapply(b_fit, function(i) i$m_f$coefficients$count["wqs"])
      se_count <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$count["wqs", "Std. Error"])
      stat_count <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$count["wqs", paste0(ts, " value")])
      p_val_count <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$count["wqs", gsub("x", ts, "Pr(>|x|)")])
      if("wqs" %in% names(b_fit[[1]]$m_f$coefficients$zero)){
        b1_zero <- sapply(b_fit, function(i) i$m_f$coefficients$zero["wqs"])
        se_zero <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$zero["wqs", "Std. Error"])
        stat_zero <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$zero["wqs", paste0(ts, " value")])
        p_val_zero <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$zero["wqs", gsub("x", ts, "Pr(>|x|)")])
      }
      else b1_zero <- se_zero <- stat_zero <- p_val_zero <- NULL
    }
    else{
      b1 <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", "Estimate"])
      se <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", "Std. Error"])
      stat <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", paste0(ts, " value")])
      p_val <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")])
    }
    n_levels <- 1
  }

  n_non_conv = sum(conv == 2)
  if(n_non_conv == 0 & control$trace) cat(paste0("The optimization function always converged\n"))
  else if(n_non_conv == b) stop("The optimization function never converged\n")
  else if(control$trace) cat(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))

  # estimate mean weight for each component (exclude weights from iterations with failed convergence)
  if (family$family == "multinomial"){
    bres <- Map(cbind, wght_matrix, b1, se, stat, p_val)
    bres <- lapply(bres, as.data.frame)
    bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", "stat", "p_val"))
    strata_names <- gsub("wqs_", "", rownames(b_fit[[1]]$m_f$sum_stat)[wqs_site])
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

  par_model_out <- list(bres, conv, bindex, nfuneval, n_levels, strata_names)
  names(par_model_out) <- c("bres", "conv", "bindex", "nfuneval", "n_levels", "strata_names")

  return(par_model_out)
}


# function to estimate mean weights for each component
mean_weight_f = function(mix_name, bres, conv, b1_pos, family, n_levels, strata_names, zero_infl){

  if (family$family == "multinomial"){
    mean_weight <- lapply(1:(n_levels-1), function(i){
      if(b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[[i]][bres[[i]]$b1 > 0 & conv!=2, "stat"]))
      else if(!b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[[i]][bres[[i]]$b1 < 0 & conv!=2, "stat"]))
      if (all(is.nan(w_t)))
        stop(paste0("There are no ", ifelse(b1_pos[i], "positive", "negative"), " b1 in the bootstrapped models for ", strata_names[i]))
      return(w_t)
    })
    mean_weight <- list.cbind(mean_weight)
  }
  else{
    if(zero_infl){
      if(b1_pos) mean_weight = apply(bres[bres$b1_count > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1_count > 0 & conv!=2, "stat_count"]))
      else mean_weight = apply(bres[bres$b1_count < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1_count < 0 & conv!=2, "stat_count"]))
    }
    else{
      if(b1_pos) mean_weight = apply(bres[bres$b1 > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1 > 0 & conv!=2, "stat"]))
      else mean_weight = apply(bres[bres$b1 < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1 < 0 & conv!=2, "stat"]))
    }
    if(all(is.nan(mean_weight)))
      stop("There are no ", ifelse(b1_pos, "positive", "negative"), " b1 in the bootstrapped models")
  }

  return(mean_weight)
}


# function to build predictive model in case of binomial dist
predict_f = function(Q, dtf, w, m_f, formula){

  dtf$wqs = as.numeric(Q%*%w)
  pred = predict(m_f, newdata = dtf, type = "response")
  y = model.response(model.frame(formula, dtf), "any")
  df_pred = data.frame(y = y, p_y = pred, stringsAsFactors = TRUE)

  return(df_pred)
}

# Function for plotting the graphs
plots <- function(data_plot, y_adj_wqs_df, q, mix_name, mean_weight, fit, family, n_levels, strata_names, df_roc, zero_infl){

  if(family$family == "multinomial"){
    data_plot = data_plot[order(data_plot[, strata_names[1]]),]
    pos = match(data_plot$mix_name, sort(mix_name))
    data_plot$mix_name = factor(data_plot$mix_name, levels(data_plot$mix_name)[pos])
    data_plot_l = melt(data_plot, id.vars = "mix_name")
    bar_plot_h = ggplot(data_plot_l, aes(x = mix_name, y = value, fill = mix_name)) +
      facet_wrap(~ variable)
  }
  else {
    data_plot <- data_plot[order(data_plot$mean_weight),]
    data_plot$mix_name <- factor(data_plot$mix_name, levels = data_plot$mix_name)
    bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name))
  }

  bar_plot_h = bar_plot_h + geom_bar(stat = "identity", color = "black") + theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(color='black'),
          legend.position = "none") + coord_flip()

  print(bar_plot_h)

  y_labs = ifelse(family$family %in% c("multinomial", "binomial"), "y", "y_adj")

  yadj_vs_wqs = ggplot(y_adj_wqs_df, aes_string("wqs", y_labs)) +
    geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()

  if(family$family == "multinomial") yadj_vs_wqs = yadj_vs_wqs + facet_wrap(~ level)

  print(yadj_vs_wqs)

  if(!(family$family %in% c("binomial", "multinomial"))){
    if(zero_infl) fit_df = data.frame(.fitted = fit$fitted.values, .resid = fit$residuals, stringsAsFactors = TRUE)
    else fit_df = augment(fit)
    res_vs_fitted = ggplot(fit_df, aes_string(x = ".fitted", y = ".resid")) + geom_point() + theme_bw() +
      xlab("Fitted values") + ylab("Residuals")
    print(res_vs_fitted)
  }

  if(n_levels == 3){
    w1_vs_w2 = ggplot(data_plot, aes_string(names(data_plot)[2], names(data_plot)[3])) + geom_point() +
      theme_bw() + xlab(names(data_plot)[2]) + ylab(names(data_plot)[3]) + geom_abline(linetype = 2) +
      geom_text_repel(aes(label=mix_name))
    print(w1_vs_w2)
  }

  if(!is.null(df_roc) & family$family == "binomial"){
    if(class(df_roc$y) == "character") df_roc$y = factor(df_roc$y)
    if(class(df_roc$y) == "factor") df_roc$y <- as.numeric(df_roc$y != levels(df_roc$y)[1])
    gg_roc = suppressWarnings(ggplot(df_roc, aes_string(d="y", m="p_y")) + geom_roc(n.cuts = 0) +
                                style_roc(xlab = "1 - Specificity", ylab = "Sensitivity"))
    auc_est = calc_auc(gg_roc)
    gg_roc = gg_roc + annotate("text", x=0.75, y=0.25, label=paste0("AUC = ", round(auc_est[, "AUC"], 3)))

    print(gg_roc)
  }
}

# Function for creating the file containing the tables
tables <- function(final_weight, mf, family, n_levels, zero_infl){

  if(family$family == "multinomial"){
    final_weight[, c(2:n_levels)] = signif(final_weight[, c(2:n_levels)], 3)
    mf_df = signif(mf$sum_stat, 3)
  }
  else{
    final_weight <- data.frame(Mix_name = final_weight$mix_name, Final_weight = signif(final_weight$mean_weight, 3), stringsAsFactors = TRUE)
    if(zero_infl){
      mf_df_count = as.data.frame(signif(coef(summary(mf))$count, 3))
      mf_df_zero = as.data.frame(signif(coef(summary(mf))$zero, 3))
      rownames(mf_df_zero) = paste0("z_", rownames(mf_df_zero))
      mf_df = rbind(mf_df_count, mf_df_zero)
    }
    else mf_df = as.data.frame(signif(coef(summary(mf)), 3))
  }

  if(zero_infl) print(kable_styling(kable(mf_df, row.names = TRUE)) %>%
                        group_rows("Count model", 1, dim(mf_df_count)[1]) %>%
                        group_rows("Zero-inflation model", dim(mf_df_count)[1]+1, dim(mf_df_count)[1]+dim(mf_df_zero)[1]))
  else print(kable_styling(kable(mf_df, row.names = TRUE)))
  print(kable_styling(kable(final_weight, row.names = FALSE)))
}
