
# starting message
.onAttach <- function(...) packageStartupMessage("Welcome to Weighted Quantile Sum (WQS) Regression.\nIf you are using a Mac you have to install XQuartz.\nYou can download it from: https://www.xquartz.org/\n")


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
  if(is.null(ls)) stop("'stratified' must be factor\n")
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


# function to split the dataset
create_rindex <- function(dtf, N, validation, valid_var, m, family){

  if(!m){
    if(!is.numeric(validation) | validation < 0 | validation >= 1) stop("'validation' must be numeric >= 0 and < 1\n")
    groups <- rep(0, N)
    if(validation > 0) groups[sample(1:N, round(N*validation))] <-1
  }
  else{
    groups = unlist(dtf[, valid_var, drop = FALSE])
    if(!any(unique(groups) %in% c(0, 1))) stop("valid_var values must be 0 or 1\n")
    if(!(0 %in% unique(groups))) stop(("0 must identify test dataset\n"))
  }
  it = which(groups == 0)
  iv = which(groups == 1)
  if(length(iv) == 0) iv = it

  indexl = list(it, iv)
  names(indexl) = c("it", "iv")
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


# function to define the parameters initial values
values.0 = function(kw, Xnames, kx, n_levels, formula, ff, wghts, bdtf, stratified, b1_pos, family, zilink, zero_infl){

  w = rep(1/kw*ifelse(is.null(stratified), 1, n_levels), kw)
  w <- sqrt(w)

  if(family$family == "multinomial"){
    fit = multinom(formula, bdtf, trace = F, weights = wghts)
    bj = c(sapply(1:(n_levels-1), function(i) coef(fit)[i,]))
    names(bj) <- Xnames
    w = rep(w, (n_levels-1))
    val.0 = c(bj, w)
  }
  else{
    if(family$family == "negbin") fit = suppressWarnings(glm.nb(formula, bdtf, weights = wghts))
    else fit = glm(formula, bdtf, family = family, weights = wghts)
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
optim.f <- function(i, objfn, Y, Xm, Q, offset, wghts, initp, n_levels, level_names, wqsvars, b1_pos, b1_constr,
                    n_vars, family, rs, zilink, zero_infl, formula, ff, kx, kw, Xnames, stratified, b, optim.method, control){

  if(rs){
    bindex <- 1:nrow(Xm)
    slctd_vars <- sample(colnames(Q), n_vars, replace=FALSE)
    initp <- initp[c(rep(T, ncol(Xm)), if(family$family == "multinomial") rep(colnames(Q) %in% slctd_vars, (n_levels-1)) else colnames(Q) %in% slctd_vars)]
    kw <- length(slctd_vars)
  }
  else{
    if(b == 1) bindex <- 1:nrow(Xm)
    else bindex <- sample(1:nrow(Xm), nrow(Xm), replace=TRUE)
    slctd_vars <- colnames(Q)
  }
  bXm <- Xm[bindex,]
  bQ <- Q[bindex, slctd_vars]
  bY <- if(family$family == "multinomial")  Y[bindex,] else Y[bindex]
  bwghts <- wghts[bindex]
  boffset <- offset[bindex]
  opt_res <- tryCatch(optim(par = initp, fn = objfn, method = optim.method, control = control,
                            kw = kw, bXm = bXm, bY = bY, boffset = boffset, bQ = bQ, kx = kx,
                            Xnames = Xnames, n_levels = n_levels, level_names = level_names,
                            wqsvars = wqsvars, family = family, zilink = zilink, zero_infl = zero_infl,
                            formula = formula, ff = ff, bwghts = bwghts, stratified = stratified, b1_pos = b1_pos,
                            b1_constr = b1_constr), error = function(e) NULL)

  if(!is.null(opt_res)) {
    if(family$family == "multinomial") par_opt <- apply(matrix(opt_res$par[(kx + 1):length(initp)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
    else par_opt <- (opt_res$par[(kx + 1):(kx + kw)]^2)/sum(opt_res$par[(kx + 1):(kx + kw)]^2)
    conv <- opt_res$convergence
    counts <- opt_res$counts
    val <- opt_res$val
    mex <- opt_res$message
  }
  else {
    if(family$family == "multinomial") par_opt <- apply(matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
    else par_opt <- (initp[(kx + 1):(kx + kw)]^2)/sum(initp[(kx + 1):(kx + kw)]^2)
    conv <- 1
    counts <- val <- mex <- NULL
  }
  if(any(is.infinite(par_opt))){
    if(family$family == "multinomial") par_opt[is.infinite(par_opt)] <- 1/colSums(is.infinite(par_opt))
    else par_opt[is.infinite(par_opt)] <- 1/sum(is.infinite(par_opt))
    par_opt[!(is.infinite(par_opt))] <- 0
  }
  mfit <- model.fit(w = par_opt, dt = bXm, bQ = bQ, Y = bY, family = family, zilink = zilink, formula = formula, ff = ff,
                    wghts = bwghts, offset = boffset, initp = initp, Xnames = Xnames, n_levels = n_levels,
                    level_names = level_names, wqsvars = wqsvars, stratified = stratified, b1_pos = b1_pos,
                    zero_infl = zero_infl, kx = kx, kw = kw)

  out <- list(par_opt = par_opt, conv = conv, counts = counts, val = val, mex = mex, mfit = mfit, bindex = bindex, slctd_vars = slctd_vars)
  return(out)
}


# function that fit the wqs model
model.fit <- function(w, dt, bQ, Y, family, zilink, formula, ff, wghts, offset, initp, Xnames, n_levels, level_names,
                      wqsvars, stratified, b1_pos, zero_infl, kx, kw){

  if(is.matrix(dt)){
    dtf <- as.data.frame(dt)
    formula <- ff <- Y ~ 0 + .
    zero_infl <- FALSE
  }
  else{
    if(family$family == "multinomial"){
      fm_l <- sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
      dtfl <- lapply(fm_l, function(i) model.matrix(i, data = dt))
      dtf <- do.call("cbind", dtfl)
    }
    else dtf <- dt
  }

  if (family$family == "multinomial"){
    w <- matrix(w, kw, n_levels-1)
    wqs <- dtf[, wqsvars] <- bQ%*%w
    colnames(wqs) <- paste0(level_names[-1], "_vs_", level_names[1])
    initp <- initp[1:kx]
  }
  else wqs <- dtf[,"wqs"] <- as.numeric(bQ%*%w)

  if(zero_infl) m_f <- zeroinfl(ff, dtf, dist = family$family, link = zilink$name)
  else{
    if(family$family == "multinomial") {
      LL <- function(p){
        tmp <- lapply(1:(n_levels-1), function(i) as.matrix(p[((i-1)*kx/(n_levels-1)+1):(i*kx/(n_levels-1))]))
        b_covs <- as.matrix(bdiag(tmp))
        term = as.matrix(dtf)%*%b_covs + offset
        -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*wghts)
      }
      nlm_out = nlm(LL, initp, hessian = TRUE)

      Estimate = nlm_out$estimate
      Standard_Error = sqrt(diag(solve(nlm_out$hessian)))
      stat = Estimate/Standard_Error
      p_value = 2*pnorm(-abs(stat))
      coefficients <- data.frame(Estimate = Estimate, Standard_Error = Standard_Error, stat = stat, p_value = p_value)
      rownames(coefficients) = Xnames
      m_f = list(nlm_out, coefficients)
      names(m_f) = c("nlm_out", "coefficients")
    }
    else if(family$family == "negbin") m_f = glm.nb(formula, data = dtf, weights = wghts)
    else m_f = glm(formula, data = dtf, family = family, weights = wghts)
  }

  mf_out = list(wqs = wqs, m_f = m_f)
  return(mf_out)
}


predictmultinom <- function(object, data, sumtype, type){
  fm_l <- sapply(colnames(object$wqs), function(i) as.formula(gsub("wqs", i, format(object$formula))))
  Xl <- lapply(fm_l, function(i) model.matrix(i, data = data))
  predl <- lapply(1:length(Xl), function(i) exp(Xl[[i]]%*%coef(object, sumtype)[i,]))
  predl <- lapply(predl, function(i) i/(1+Reduce("+", predl)))
  pred <- do.call("cbind", predl)
  pred <- cbind(1-rowSums(pred), pred)
  y <- model.response(model.frame(fm_l[[1]], data), "any")
  if(type == "response"){
    Ylevels <- levels(y)
    pred <- factor(Ylevels[apply(pred, 1, which.max)], levels = Ylevels)
  }
  else if(!(type %in% c("response", "prob"))) stop("If family is \"multinomial\" then predict type must be \"response\" or \"prob\"\n")
  return(list(pred = pred, y = y))
}

set_par_names <- function(i, slctd_vars, param, q_name, family){

  temp <- param[[i]]$par_opt
  if(family$family == "multinomial"){
    param[[i]]$par_opt <- matrix(NA, length(q_name), dim(temp)[2])
    param[[i]]$par_opt[which(q_name %in% slctd_vars[[i]]),] <- temp
    rownames(param[[i]]$par_opt) <- q_name
  }
  else{
    param[[i]]$par_opt <- rep(NA, length(q_name))
    names(param[[i]]$par_opt) <- q_name
    param[[i]]$par_opt[slctd_vars[[i]]] <- temp
  }

  return(param[[i]])
}

