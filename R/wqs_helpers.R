
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
values.0 = function(kw, Xnames, kx, formula, ff, wghts, bdtf, stratified, stratlev, b1_pos, family, wp, wn, dwqs, zilink, zero_infl){

  w = rep(1/kw*stratlev, kw)
  w <- sqrt(w)

  if(family$family == "multinomial"){
    n_levels <- nlevels(eval(formula[[2]], envir = bdtf))
    fit = multinom(formula, bdtf, trace = F, weights = wghts)
    bj = c(sapply(1:(n_levels-1), function(i) coef(fit)[i,]))
    names(bj) <- Xnames
    if(dwqs) val.0 = c(bj, sqrt(wp), sqrt(wn))
    else{
      w = rep(w, (n_levels-1))
      val.0 = c(bj, w)
    }
  }
  else{
    if(family$family == "negbin") fit = suppressWarnings(glm.nb(formula, bdtf, weights = wghts))
    else fit = glm(formula, bdtf, family = family, weights = wghts)
    bj = coef(fit)
    if(dwqs) val.0 = c(bj, sqrt(wp), sqrt(wn))
    else val.0 = c(bj, w)
    if(family$family == "negbin"){
      if(length(attr(terms(formula), "term.labels")) == 1) val.0 = c(val.0, 1)
      else val.0 = c(val.0, log(fit$theta))
    }
  }

  if(dwqs){
    pwqs_site <- which(grepl("pwqs", Xnames))
    nwqs_site <- which(grepl("nwqs", Xnames))
    val.0[pwqs_site] = 0.0001
    val.0[nwqs_site] = -0.0001
  }
  else{
    wqs_site <- which(grepl("wqs", Xnames))
    val.0[wqs_site] = sapply(b1_pos, function(i) ifelse(i, 0.0001, -0.0001))
  }

  return(val.0)
}


# optimization function to estimate the weights
optim.f <- function(i, objfn, Y, Xm, Q, offset, wghts, initp,
                    b1_pos, b_constr, n_vars, dwqs, family, rs, zilink, zero_infl, formula, ff, kx, kw, Xnames,
                    stratified, b, optim.method, control, lambda, stratlev, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont){

  if(rs){
    bindex <- 1:nrow(Xm)
    slctd_vars <- sample(colnames(Q), n_vars, replace=FALSE)
    if(dwqs) initp <- initp[c(rep(T, ncol(Xm)), rep(colnames(Q) %in% slctd_vars, 2))]
    else initp <- initp[c(rep(T, ncol(Xm)), colnames(Q) %in% slctd_vars)]
    kw <- length(slctd_vars)
  }
  else{
    if(b == 1) bindex <- 1:nrow(Xm)
    else bindex <- sample(1:nrow(Xm), nrow(Xm), replace=TRUE)
    slctd_vars <- colnames(Q)
  }
  bXm <- Xm[bindex,]
  bQ <- Q[bindex, slctd_vars]
  bY <- Y[bindex]
  bwghts <- wghts[bindex]
  boffset <- offset[bindex]
  opt_res <- optim(par = initp, fn = objfn, method = optim.method, control = control,
                            kw = kw, bXm = bXm, bY = bY, boffset = boffset, bQ = bQ, kx = kx,
                            Xnames = Xnames, family = family, dwqs = dwqs,
                            zilink = zilink, zero_infl = zero_infl, formula = formula, ff = ff, bwghts = bwghts,
                            stratified = stratified, stratlev = stratlev, b1_pos = b1_pos, b_constr = b_constr,
                            lambda = lambda, intpos = intpos, wqsint = wqsint, intvars = intvars, bint_cont_pos = bint_cont_pos, bint_cat_pos = bint_cat_pos, intcont = intcont)
  if(!is.null(opt_res)) {
    if(dwqs) par_opt <- c(opt_res$par[(kx + 1):(kx + kw)]^2/sum(opt_res$par[(kx + 1):(kx + kw)]^2),
                        opt_res$par[(kx + kw + 1):(kx + 2*kw)]^2/sum(opt_res$par[(kx + kw + 1):(kx + 2*kw)]^2))
    else{
        if(!is.null(stratified)){
          par_opt <- lapply(1:stratlev, function(i){
            ncolQ <- ncol(bQ)/stratlev
            tmp <- (opt_res$par[(kx + 1 + (i-1)*ncolQ):(kx + kw/stratlev*i)]^2)/sum(opt_res$par[(kx + 1 + (i-1)*ncolQ):(kx + kw/stratlev*i)]^2)
            return(tmp)
          })
          par_opt <- do.call("c", par_opt)
        }
        else par_opt <- (opt_res$par[(kx + 1):(kx + kw)]^2)/sum(opt_res$par[(kx + 1):(kx + kw)]^2)
    }
    conv <- opt_res$convergence
    counts <- opt_res$counts
    val <- opt_res$val
    mex <- opt_res$message
  }
  else{
    if(dwqs) par_opt <- c(initp[(kx + 1):(kx + kw)]^2/sum(initp[(kx + 1):(kx + kw)]^2),
                        initp[(kx + kw + 1):(kx + 2*kw)]^2/sum(initp[(kx + kw + 1):(kx + 2*kw)]^2))
    else par_opt <- (initp[(kx + 1):(kx + kw)]^2)/sum(initp[(kx + 1):(kx + kw)]^2)
    conv <- 1
    counts <- val <- mex <- NULL
  }
  if(any(is.infinite(par_opt))){
    if(dwqs){
        par_opt[is.infinite(par_opt[1:kw])] <- 1/sum(is.infinite(par_opt[1:kw]))
        par_opt[!(is.infinite(par_opt[1:kw]))] <- 0
        par_opt[(kw+1):2*kw][is.infinite(par_opt[(kw+1):2*kw])] <- 1/sum(is.infinite(par_opt[(kw+1):2*kw]))
        par_opt[(kw+1):2*kw][!(is.infinite(par_opt[(kw+1):2*kw]))] <- 0
    }
    else{
        par_opt[is.infinite(par_opt)] <- 1/sum(is.infinite(par_opt))
        par_opt[!(is.infinite(par_opt))] <- 0
    }
  }
  mfit <- model.fit(w = par_opt, dt = bXm, bQ = bQ, Y = bY, family = family, dwqs = dwqs, zilink = zilink,
                    formula = formula, ff = ff, wghts = bwghts, offset = boffset, initp = initp, Xnames = Xnames,
                    stratified = stratified, b1_pos = b1_pos, zero_infl = zero_infl, kx = kx,
                    kw = kw, intpos = intpos, wqsint = wqsint, intvars = intvars)

  out <- list(par_opt = par_opt, conv = conv, counts = counts, val = val, mex = mex, mfit = mfit, bindex = bindex, slctd_vars = slctd_vars)
  return(out)
}


optim.f_multinom <- function(i, objfn, Y, Xm, Q, offset, wghts, initp, n_levels, level_names, wqsvars, pwqsvars,
                             nwqsvars, b1_pos, b_constr, n_vars, dwqs, family, rs, formula, kx, kw, Xnames,
                             stratified, b, optim.method, control, lambda, stratlev, intpos, wqsint, intvars){

  if(rs){
    bindex <- 1:nrow(Xm)
    slctd_vars <- sample(colnames(Q), n_vars, replace=FALSE)
    if(dwqs) initp <- initp[c(rep(T, ncol(Xm)), rep(rep(colnames(Q) %in% slctd_vars, (n_levels-1)), 2))]
    else initp <- initp[c(rep(T, ncol(Xm)), rep(colnames(Q) %in% slctd_vars, (n_levels-1)))]
    kw <- length(slctd_vars)
  }
  else{
    if(b == 1) bindex <- 1:nrow(Xm)
    else bindex <- sample(1:nrow(Xm), nrow(Xm), replace=TRUE)
    slctd_vars <- colnames(Q)
  }
  bXm <- Xm[bindex,]
  bQ <- Q[bindex, slctd_vars]
  bY <- Y[bindex,]
  bwghts <- wghts[bindex]
  boffset <- offset[bindex]
  opt_res <- optim(par = initp, fn = objfn, method = optim.method, control = control,
                            kw = kw, bXm = bXm, bY = bY, boffset = boffset, bQ = bQ, kx = kx,
                            Xnames = Xnames, n_levels = n_levels, level_names = level_names,
                            wqsvars = wqsvars, pwqsvars = pwqsvars, nwqsvars = nwqsvars, dwqs = dwqs,
                            formula = formula, bwghts = bwghts,
                            stratified = stratified, stratlev = stratlev, b1_pos = b1_pos, b_constr = b_constr,
                            lambda = lambda, intpos = intpos, wqsint = wqsint, intvars = intvars)
  if(!is.null(opt_res)) {
    if(dwqs) par_opt <- rbind(apply(matrix(opt_res$par[(kx + 1):(kx + kw*(n_levels-1))]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)),
                       apply(matrix(opt_res$par[(kx + kw*(n_levels-1) + 1):(kx + 2*(n_levels-1)*kw)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)))
    else par_opt <- apply(matrix(opt_res$par[(kx + 1):length(initp)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
    conv <- opt_res$convergence
    counts <- opt_res$counts
    val <- opt_res$val
    mex <- opt_res$message
  }
  else{
    if(dwqs) par_opt <- rbind(apply(matrix(initp[(kx + 1):(kx + kw*(n_levels-1))]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)),
                              apply(matrix(initp[(kx + kw*(n_levels-1) + 1):(kx + 2*(n_levels-1)*kw)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i)))
    else par_opt <- apply(matrix(initp[(kx + 1):length(initp)]^2, kw, n_levels-1), MARGIN = 2, FUN = function(i) i/sum(i))
    conv <- 1
    counts <- val <- mex <- NULL
  }
  if(any(is.infinite(par_opt))){
    if(dwqs){
      par_opt[which(is.infinite(par_opt[1:kw,]), arr.ind = T)] <- 1/colSums(is.infinite(par_opt[1:kw,])[,which(is.infinite(par_opt[1:kw,]), arr.ind = T)[,2]])
      par_opt[which(!is.infinite(par_opt[1:kw,]), arr.ind = T)] <- 0
      par_opt[(kw+1):2*kw,][which(is.infinite(par_opt[(kw+1):2*kw,]), arr.ind = T)] <- 1/colSums(is.infinite(par_opt[(kw+1):2*kw,])[,which(is.infinite(par_opt[(kw+1):2*kw,]), arr.ind = T)[,2]])
      par_opt[(kw+1):2*kw,][which(!is.infinite(par_opt[(kw+1):2*kw,]), arr.ind = T)] <- 0
    }
    else{
      par_opt[which(is.infinite(par_opt), arr.ind = T)] <- 1/colSums(is.infinite(par_opt)[,which(is.infinite(par_opt), arr.ind = T)[,2]])
      par_opt[which(!is.infinite(par_opt), arr.ind = T)] <- 0
    }
  }
  mfit <- model.fit_multinom(w = par_opt, dt = bXm, bQ = bQ, Y = bY, dwqs = dwqs,
                    formula = formula, wghts = bwghts, offset = boffset, initp = initp, Xnames = Xnames,
                    n_levels = n_levels, level_names = level_names, wqsvars = wqsvars, pwqsvars = pwqsvars,
                    nwqsvars = nwqsvars, stratified = stratified, b1_pos = b1_pos, kx = kx,
                    kw = kw, intpos = intpos, wqsint = wqsint, intvars = intvars)

  out <- list(par_opt = par_opt, conv = conv, counts = counts, val = val, mex = mex, mfit = mfit, bindex = bindex, slctd_vars = slctd_vars)
  return(out)
}


# function that fit the wqs model
model.fit <- function(w, dt, bQ, Y, family, dwqs, zilink, formula, ff, wghts, offset, initp, Xnames,
                      stratified, b1_pos, zero_infl, kx, kw, intpos, wqsint, intvars){

  if(is.matrix(dt)){
    dtf <- as.data.frame(dt[,which(colnames(dt)!="(Intercept)")])
    formula <- ff <- Y ~ .
    zero_infl <- FALSE
    intpos <- grepl("wqs", colnames(dtf)) & grepl(":", colnames(dtf))
  }
  else dtf <- dt

    if(dwqs){
      pwqs <- dtf[,"pwqs"] <- as.numeric(bQ%*%w[1:kw])
      nwqs <- dtf[,"nwqs"] <- as.numeric(bQ%*%w[(kw+1):(2*kw)])
      wqs <- NULL
    }
    else{
      wqs <- dtf[,"wqs"] <- as.numeric(bQ%*%w)
      if(wqsint & is.matrix(dt)) dtf[,intpos] <- wqs*dtf[,intvars]
      pwqs <- nwqs <- NULL
    }

  if(zero_infl) m_f <- zeroinfl(ff, dtf, dist = family$family, link = zilink$name)
  else{
    if(family$family == "negbin") m_f = glm.nb(formula, data = dtf, weights = wghts)
    else m_f = glm(formula, data = dtf, family = family, weights = wghts)
  }

  mf_out = list(wqs = wqs, pwqs = pwqs, nwqs = nwqs, m_f = m_f)
  return(mf_out)
}


model.fit_multinom <- function(w, dt, bQ, Y, dwqs, formula, wghts, offset, initp, Xnames, n_levels, level_names,
                      wqsvars, pwqsvars, nwqsvars, stratified, b1_pos, kx, kw, intpos, wqsint, intvars){

  if(is.matrix(dt)){
    dtf <- as.data.frame(dt)
    formula <- Y ~ .
  }
  else{
    if(dwqs){
      tmp <- sapply(1:length(pwqsvars), function(i) gsub("pwqs", pwqsvars[i], format(formula)))
      fm_l <- sapply(1:length(nwqsvars), function(i) as.formula(gsub("nwqs", nwqsvars[i], tmp[i])))
    }
    else fm_l <- sapply(wqsvars, function(i) as.formula(gsub("wqs", i, format(formula))))
    dtfl <- lapply(fm_l, function(i) model.matrix(i, data = dt))
    dtf <- do.call("cbind", dtfl)
  }

  if(dwqs){
    wp <- w[1:kw, 1:(n_levels-1)]
    wn <- w[(kw+1):(2*kw), 1:(n_levels-1)]
    pwqs <- dtf[, pwqsvars] <- bQ%*%wp
    nwqs <- dtf[, nwqsvars] <- bQ%*%wn
    colnames(pwqs) <- paste0("pwqs_", level_names[-1], "_vs_", level_names[1])
    colnames(nwqs) <- paste0("nwqs_", level_names[-1], "_vs_", level_names[1])
    wqs <- NULL
  }
  else{
    w <- matrix(w, kw, n_levels-1)
    wqs <- dtf[, wqsvars] <- bQ%*%w
    colnames(wqs) <- paste0("wqs_", level_names[-1], "_vs_", level_names[1])
    pwqs <- nwqs <- NULL
  }
  initp <- initp[1:kx]

    LL <- function(p, dtfm){
      kx <- length(p)
      tmp <- lapply(1:(n_levels-1), function(i) as.matrix(p[((i-1)*kx/(n_levels-1)+1):(i*kx/(n_levels-1))]))
      b_covs <- as.matrix(bdiag(tmp))
      term = dtfm%*%b_covs + offset
      -sum((diag(Y%*%t(term)) - log(1 + rowSums(exp(term))))*wghts)
    }
    nlm_out = nlm(LL, initp, hessian = TRUE, dtfm = as.matrix(dtf))

    Estimate = nlm_out$estimate
    Standard_Error = sqrt(diag(solve(nlm_out$hessian)))
    stat = Estimate/Standard_Error
    p_value = 2*pnorm(-abs(stat))
    coefficients <- data.frame(Estimate = Estimate, Standard_Error = Standard_Error, stat = stat, p_value = p_value)
    rownames(coefficients) = Xnames
    m_f = list(nlm_out, coefficients)
    names(m_f) = c("nlm_out", "coefficients")

  mf_out = list(wqs = wqs, pwqs = pwqs, nwqs = nwqs, m_f = m_f)
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

set_par_names <- function(i, slctd_vars, param, q_name, family, dwqs){

  temp <- param[[i]]$par_opt
  if(family$family == "multinomial"){
    param[[i]]$par_opt <- matrix(NA, length(q_name), dim(temp)[2])
    param[[i]]$par_opt[which(q_name %in% slctd_vars[[i]]),] <- temp
    rownames(param[[i]]$par_opt) <- q_name
  }
  else{
    param[[i]]$par_opt <- rep(NA, ifelse(dwqs, 2, 1)*length(q_name))
    names(param[[i]]$par_opt) <- rep(q_name, ifelse(dwqs, 2, 1))
    if(dwqs){
      param[[i]]$par_opt[1:length(q_name)][slctd_vars[[i]]] <- temp[1:(length(temp)/2)]
      param[[i]]$par_opt[(length(q_name)+1):(2*length(q_name))][slctd_vars[[i]]] <- temp[(length(temp)/2+1):length(temp)]
    }
    else param[[i]]$par_opt[slctd_vars[[i]]] <- temp
  }

  return(param[[i]])
}

