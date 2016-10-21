
# function to check datatset
.check.function = function(formula, mix_name, data, q, validation, valid_var, b, b1_pos, family, seed,
                           wqs2, plots, tables){

  form = ifelse(class(formula) != "formula", TRUE, FALSE)
  if(form) stop("formula must be of class formula")

  chars = ifelse(!is.character(mix_name), TRUE, FALSE)
  if(chars) stop("mix_name must be a character")

  df = ifelse(!is.data.frame(data), TRUE, FALSE)
  if(df) stop("data must be a data.frame")

  quant = ifelse((!is.numeric(q) | length(q) != 1) & !is.null(q), TRUE, FALSE)
  if(quant) stop("q must be numeric of length 1, or a NULL")

  val = ifelse(!is.numeric(validation) | validation < 0 | validation >= 1, TRUE, FALSE)
  if(val) stop("validation must be numeric between 0 and 1, or 0")

  v_var = ifelse(!is.character(valid_var) & !is.null(valid_var), TRUE, FALSE)
  if(v_var) stop("valid_var must be a character or NULL")

  boot = ifelse(!is.numeric(b) | length(b) != 1, TRUE, FALSE)
  if(boot) stop("b must be numeric of length 1")

  b1p = ifelse(!is.logical(b1_pos), TRUE, FALSE)
  if(b1p) stop("b1_pos must be logical")

  fam = ifelse(family != "gaussian" & family != "binomial", TRUE, FALSE)
  if(fam) stop("family must be gaussian or binomial")

  sd = ifelse(!is.null(seed) & !is.numeric(seed), TRUE, FALSE)
  if(sd) stop("seed must be numeric or NULL")

  wqs_quad = ifelse(!is.logical(wqs2), TRUE, FALSE)
  if(wqs_quad) stop("wqs2 must be logical")

  plts = ifelse(!is.logical(plots), TRUE, FALSE)
  if(plts) stop("plots must be logical")

  tbls = ifelse(!is.logical(tables), TRUE, FALSE)
  if(tbls) stop("tables must be logical")
}


# function to create variables with quantile of the components
quantile_f <- function(data, v_n, q){

  for (i in 1:length(v_n)){
    dat_num = as.numeric(unlist(data[, v_n[i]]))
    data[[paste(v_n[i], "q", sep = "_")]] = cut(dat_num, breaks = unique(quantile(dat_num, probs = seq(0, 1, by = 1/q), na.rm = TRUE)),
                                                labels = FALSE, include.lowest = TRUE) - 1
  }

  return(data)
}


# function to split the dataset
split_f <- function(data, validation, seed){

  if (validation == 0) {
    data_t = data
    data_v = data
  }
  else{
    data$group = 0
    data$group[rownames(data) %in% sample(rownames(data), round(nrow(data)*validation))] = 1
    data_t = data[data$group == 0,]
    data_v = data[data$group == 1,]
  }

  split_out =list(data_t, data_v)
  names(split_out) = c("data_t", "data_v")

  return(split_out)
}


# functtion to define the objective function
objfn <- function(initp, Q, y, family, covrts){

  b0 = initp[1]
  b1 = initp[2]
  w = initp[3:(dim(Q)[2]+2)]
  cvrt_coeff = initp[(dim(Q)[2]+3):length(initp)]

  if(dim(covrts)[2] == 0) term = b0 + b1*Q%*%w
  else term = b0 + b1*Q%*%w + covrts%*%cvrt_coeff

  if (family == "gaussian") f = sum((y - term)^2)
  else if (family == "binomial") {
    p = 1/(1 + exp(-term))
    f = -sum(y * log(p) + (1 - y) * log(1 - p))
  }

  return(f)
}


# function to define the equality constraint
linconst = function(initp, Q, y, family, covrts){

  wsum=sum(initp[3:(dim(Q)[2] + 2)])

  return(wsum)
}


# function to determine bounded prameters
bounded_param = function(initp, Q, y, family, covrts){

  b_p = initp[3:(dim(Q)[2] + 2)]

  return(b_p)
}


# function to define the lower bounds
LBound = function(Q, b1_pos){

  LB = rep(0, dim(Q)[2])

  return(LB)
}


# function to define the upper bounds
UBound = function(Q, b1_pos){

  UB = rep(1, dim(Q)[2])

  return(UB)
}


# function to define the parameters initial values
values.0 = function(lenp, y, covrts, b1_pos, family){

  y = as.matrix(y)
  b1_0 = ifelse(b1_pos, 0.001, -0.001)

  if (dim(covrts)[2] == 0){
    val.0 = c(0, b1_0, rep(1/lenp, lenp))
    names(val.0) = c("b0", "b1", paste0("w", 1:lenp))
  }
  else {
    if (family == "gaussian") fit = glm(y ~ covrts, family = gaussian(link = "identity"))
    else if (family == "binomial") fit = glm(y ~ covrts, family = binomial(link = "logit"))
    bj = coef(fit)[-1]
    val.0 = c(0, b1_0, rep(1/lenp, lenp), bj)
    names(val.0) = c("b0", "b1", paste0("w", 1:lenp), paste0("b", 2:(dim(covrts)[2]+1)))
  }

  return(val.0)
}


# optimization function to estimate the weights
optim.f <- function(Q, y, b1_pos, family, covrts){
  Q = as.matrix(Q)
  if(dim(covrts)[2] > 0) covrts = as.matrix(covrts)
  lenp = dim(Q)[2]
  initp = values.0(lenp, y, covrts, b1_pos, family)
  LowB = LBound(Q, b1_pos)
  UpB = UBound(Q, b1_pos)

  opt_res = suppressWarnings(solnp(pars = initp, fun = objfn, eqfun = linconst, eqB = 1,
                                   ineqfun = bounded_param, ineqLB = LowB, ineqUB = UpB,
                                   control = list(trace = 0), Q = Q, y = y, family = family,
                                   covrts = covrts))

  par_opt = opt_res$pars
  conv = opt_res$convergence
  nfuneval = opt_res$nfuneval

  out = list(par_opt, conv, nfuneval)
  names(out) = c("par_opt", "conv", "nfuneval")

  return(out)
}


# function that fit the wqs model
model.fit <- function(Q, y, wghts, family, covrts, wqs2){

  Q = as.matrix(Q)
  wqs = Q%*%wghts
  wqs = as.data.frame(wqs)
  names(wqs) = "wqs"

  new_data = cbind(as.data.frame(y), wqs)
  names(new_data)[1] = "y"

  if (dim(covrts)[2] > 0) {
    new_data = cbind(new_data, covrts)
    names(new_data)[c((dim(wqs)[2] + 2):dim(new_data)[2])] = names(covrts)
  }

  if (family == "gaussian") m_f = glm(y ~ ., data = new_data, family = gaussian(link = "identity"))
  else if (family == "binomial") m_f = glm(y ~ ., data = new_data, family = binomial(link = "logit"))

  mf_out = list(wqs, m_f)
  names(mf_out) = c("wqs", "m_f")

  if (wqs2 == TRUE) {
    wqs_2 = wqs^2
    wqs_2 = as.data.frame(wqs_2)
    names(wqs_2) = "wqs_2"
    wqs = cbind(wqs, wqs_2)
    new_data = cbind(new_data, wqs_2)

    if (family == "gaussian") m_f2 = glm(y ~ ., data = new_data, family = gaussian(link = "identity"))
    else if (family == "binomial") m_f2 = glm(y ~ ., data = new_data, family = binomial(link = "logit"))

    aov = anova(m_f, m_f2, test = "Chisq")

    mf_out = list(wqs, m_f, m_f2, aov)
    names(mf_out) = c("wqs", "m_f", "m_f2", "aov")
  }

  return(mf_out)
}


# function that calls the optimization function and the function to fit the model for each bootstrap sample
par.modl.est <- function(data, y_name, q_name, cov_name, b, b1_pos, family, seed){

  index_b = vector("list", b)
  param = vector("list", b)
  b1 = rep(NA, b)
  wght_matrix = matrix(NA, b, length(q_name))
  colnames(wght_matrix) = paste("w", 1:length(q_name), sep = "_")
  conv = rep(NA, b)
  p_val = rep(NA, b)
  nfuneval = rep(NA, b)

  if (family == "gaussian") ts = "t"
  else if (family == "binomial") ts = "z"

  for (i in 1:b) {

    # create bootstrap sample
    index_b[[i]] = sample(1:nrow(data), nrow(data), replace=TRUE)
    data_b = data[index_b[[i]],]

    # calculate parameters
    param[[i]] = optim.f(data_b[, q_name, drop = FALSE], data_b[, y_name, drop = FALSE], b1_pos,
                         family, data_b[, cov_name, drop = FALSE])

    # fit the wqs model for each bootstrap sample
    b_fit = model.fit(data_b[, q_name, drop = FALSE], data_b[, y_name, drop = FALSE],
                      param[[i]]$par_opt[3:(2 + length(q_name))], family,
                      data_b[, cov_name, drop = FALSE], wqs2 = FALSE)

    wght_matrix[i,] = param[[i]]$par_opt[3:(2 + length(q_name))]
    b1[i] = b_fit$m_f$coefficients[row.names = "wqs"]
    conv[i] = param[[i]]$conv
    p_val[i] = summary(b_fit$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")]
    nfuneval[i] = param[[i]]$nfuneval
  }

  par_model_out = list(wght_matrix, b1, conv, p_val, index_b, nfuneval)
  names(par_model_out) = c("wght_matrix", "b1", "conv", "p_val", "index_b", "nfuneval")

  return(par_model_out)
}


# Function for plotting the graphs
plots <- function(data_plot, y_adj_wqs_df, q, mix_name, mean_weight){

  bar_plot_h = ggplot(data_plot,aes(x = mix_name, y = mean_weight, fill = mix_name)) +
    geom_bar(stat = "identity", color = "black") + theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(color='black'),
          legend.position = "none") + coord_flip()

  print(bar_plot_h)

  yadj_vs_wqs = ggplot(y_adj_wqs_df, aes_string("wqs", "y_adj")) +
    geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()

  if (!is.null(q)) yadj_vs_wqs = yadj_vs_wqs + coord_cartesian(xlim = c(0, q-1))

  print(yadj_vs_wqs)
}

# Function for creating the file containing the tables
tables <- function(final_weight, mf, mf2, aov){

  final_weight = data.frame(Mix_name = final_weight$mix_name,
                            Final_weight = signif(final_weight$mean_weight, 3))
  write_tableHTML(tableHTML(final_weight, rownames = FALSE), file = "Weights_table.html")

  mf_df = data.frame(Coefficients = names(mf$coefficients), Estimate = signif(mf$coefficients, 3),
                     Standard_Error = signif(coef(summary(mf))[, 2], 3),
                     stat = signif(coef(summary(mf))[, 3], 3),
                     p_value = signif(coef(summary(mf))[, 4], 3))
  names(mf_df)[4] = colnames(coef(summary(mf)))[3]
  write_tableHTML(tableHTML(mf_df, rownames = FALSE), file = "Summary_results.html")
  print(ztable(mf_df, include.rownames = FALSE))

  if (!is.null(mf2) & !is.null(aov)){
    mf2_df = data.frame(Coefficients = names(mf2$coefficients), Estimate = signif(mf2$coefficients, 3),
                       Standard_Error = signif(coef(summary(mf2))[, 2], 3),
                       stat = signif(coef(summary(mf2))[, 3], 3),
                       p_val = signif(coef(summary(mf2))[, 4], 3))
    names(mf2_df)[4] = colnames(coef(summary(mf2)))[3]
    write_tableHTML(tableHTML(mf2_df, rownames = FALSE), file = "Summary_results_quadratic.html")
    print(ztable(mf2_df, include.rownames = FALSE))

    aov_df = data.frame(Model = c("Linear", "Quadratic"), Resid_df = aov[, 1],
                        Resid_Dev = signif(aov[, 2], 5), Df = aov[, 3],
                        Deviance = signif(aov[, 4], 3), p_value = signif(aov[, 5], 3))
    aov_format = format(aov_df)
    aov_format[is.na(aov_df)] = ""
    write_tableHTML(tableHTML(aov_format, rownames = FALSE), file = "Aov_results.html")
  }
}

