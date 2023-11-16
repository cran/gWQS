## ----setup, include=FALSE-----------------------------------------------------
# knitr::opts_chunk$set(cache = TRUE)

library(cowplot)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

library(gWQS)
library(ggplot2)
library(knitr)
library(kableExtra)
library(reshape2)

# we save the names of the mixture variables in the variable "PCBs"
PCBs <- names(wqs_data)[1:34]
# we run the model and save the results in the variable "results2i"
results2i <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
                  q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                  family = "gaussian", seed = 2016)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

summary(results2i)


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  # bar plot
#  gwqs_barplot(results2i)
#  # scatter plot y vs wqs
#  gwqs_scatterplot(results2i)
#  # scatter plot residuals vs fitted values
#  gwqs_fitted_vs_resid(results2i)
#  # boxplot of the weights estimated at each repeated holdout step
#  gwqs_boxplot(results2i)
#  

## ----model1, results='asis', fig.show='hold', fig.height=8, fig.width=8, echo=FALSE, fig.cap="Plots displayed for linear outcomes.", message=FALSE----

PCBs <- names(wqs_data)[1:34]
results2i <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
                q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                family = "gaussian", seed = 2016)

w_ord <- order(results2i$final_weights$`Estimate pos`)
mean_weight_pos <- results2i$final_weights$`Estimate pos`[w_ord]
mean_weight_neg <- results2i$final_weights$`Estimate neg`[w_ord]
mix_name <- factor(results2i$final_weights$mix_name[w_ord], 
                   levels = results2i$final_weights$mix_name[w_ord])
data_plot <- data.frame(mean_weight = c(mean_weight_pos, mean_weight_neg), 
                        mix_name = rep(mix_name, 2),
                        index = factor(rep(c("pwqs", "nwqs"), each = length(w_ord)), 
                                       levels = c("pwqs", "nwqs")))
bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight)) + 
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip() + ggtitle("A") +
  geom_hline(yintercept = 1/length(PCBs), linetype="dashed", color = "red") +
  facet_wrap(~ index)

yadj_vs_wqs <- ggplot(melt(results2i$y_wqs_df, measure.vars = c("pwqs", "nwqs")), aes(value, y_adj)) + 
  geom_point() + facet_wrap(~ variable) + xlab("wqs") + 
  stat_smooth(method = "loess", se = FALSE, linewidth = 1.5) + theme_bw() + ggtitle("B")

fit_df <- data.frame(fitted = fitted(results2i), 
                     resid = residuals(results2i, type = "response"))
res_vs_fitted <- ggplot(fit_df, aes(x = fitted, y = resid)) + geom_point() + 
  theme_bw() + xlab("Fitted values") + ylab("Residuals") + ggtitle("C")

wboxplot <- melt(results2i$wmat, varnames = c("rh", "mix_name"))
wboxplot$mix_name <- factor(wboxplot$mix_name, levels = results2i$final_weights$mix_name)
wboxplot$L1 <- factor(wboxplot$L1, levels = c("wmatpos", "wmatneg"), labels = c("pwqs", "nwqs"))
box_plot <- ggplot(wboxplot, aes(x = mix_name, y = value)) +
  geom_boxplot(outlier.shape = " ") + theme_bw() + facet_wrap(~ L1, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Weight") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) + 
  geom_hline(yintercept = 1/length(PCBs), linetype="dashed", color = "red") +
  geom_jitter(alpha = 0.3) + ggtitle("D")

plot_grid(bar_plot_h, yadj_vs_wqs, res_vs_fitted, box_plot, nrow = 2, ncol = 2)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

head(results2i$final_weights)


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  gwqs_summary_tab(results2i)
#  

## ----sum1, results='asis', echo=FALSE-----------------------------------------

gwqs_summary_tab(results2i, caption = "Summary results of the WQS regression for linear outcomes.")


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  mf_df <- as.data.frame(signif(coef(summary(results2i)), 3))
#  kable_styling(kable(mf_df, row.names = TRUE))
#  

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  gwqs_weights_tab(results2i)
#  

## ----w1, echo=FALSE, eval=TRUE------------------------------------------------

final_weight <- results2i$final_weights
final_weight[, -1] <- signif(final_weight[, -1], 3)
scroll_box(kable_styling(kable(final_weight, row.names = FALSE, caption = "Weights table of the WQS regression for linear outcomes.")), height = "400px")


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  final_weight <- results2i$final_weights
#  final_weight[, -1] <- signif(final_weight[, -1], 3)
#  kable_styling(kable(final_weight, row.names = FALSE))
#  

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # bar plot
#  w_ord <- order(results2i$final_weights$`Estimate pos`)
#  mean_weight_pos <- results2i$final_weights$`Estimate pos`[w_ord]
#  mean_weight_neg <- results2i$final_weights$`Estimate neg`[w_ord]
#  mix_name <- factor(results2i$final_weights$mix_name[w_ord],
#                     levels = results2i$final_weights$mix_name[w_ord])
#  data_plot <- data.frame(mean_weight = c(mean_weight_pos, mean_weight_neg),
#                          mix_name = rep(mix_name, 2),
#                          index = factor(rep(c("pwqs", "nwqs"), each = length(w_ord)),
#                                         levels = c("pwqs", "nwqs")))
#  ggplot(data_plot, aes(x = mix_name, y = mean_weight)) +
#    geom_bar(stat = "identity", color = "black") + theme_bw() +
#    theme(axis.ticks = element_blank(),
#          axis.title = element_blank(),
#          axis.text.x = element_text(color='black'),
#          legend.position = "none") + coord_flip() +
#    geom_hline(yintercept = 1/length(PCBs), linetype="dashed", color = "red") +
#    facet_wrap(~ index)
#  #
#  # scatter plot y vs wqs
#  ggplot(melt(results2i$y_wqs_df, measure.vars = c("pwqs", "nwqs")), aes(value, y_adj)) +
#    geom_point() + facet_wrap(~ variable) + xlab("wqs") +
#    stat_smooth(method = "loess", se = FALSE, linewidth = 1.5) + theme_bw()
#  #
#  # scatter plot residuals vs fitted values
#  fit_df <- data.frame(fitted = fitted(results2i),
#                       resid = residuals(results2i, type = "response"))
#  res_vs_fitted <- ggplot(fit_df, aes(x = fitted, y = resid)) + geom_point() +
#    theme_bw() + xlab("Fitted values") + ylab("Residuals")

## ---- echo=TRUE, eval=TRUE----------------------------------------------------

# we run the model setting the penalization term equal to 90
results2i_l90 <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
                    q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                    family = "gaussian", seed = 2016, lambda = 90)

# we run the model setting the penalization term equal to 900
results2i_l900 <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
                     q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                     family = "gaussian", seed = 2016, lambda = 900)

# we run the model setting the penalization term equal to 900
results2i_l9000 <- gwqs(yLBX ~ pwqs + nwqs, mix_name = PCBs, data = wqs_data, 
                      q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                      family = "gaussian", seed = 2016, lambda = 9000)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

lambda_AIC_2i <- data.frame(lambda = c(0, 90, 900, 9000),
                            AIC = c(results2i$fit$aic, results2i_l90$fit$aic, 
                                    results2i_l900$fit$aic, results2i_l9000$fit$aic))
kable(lambda_AIC_2i) %>% kable_styling()


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

summary(results2i_l90)


## ---- results='asis', fig.show='hold', fig.height=8, fig.width=8, echo=TRUE, message=FALSE----

gwqs_barplot(results2i_l90)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

results1i <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
                  q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                  family = "gaussian", seed = 2016)

# we run the model setting the penalization term equal to 90
results1i_l90 <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
                    q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                    family = "gaussian", seed = 2016, lambda = 90)

# we run the model setting the penalization term equal to 900
results1i_l900 <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
                     q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                     family = "gaussian", seed = 2016, lambda = 900)

# we run the model setting the penalization term equal to 900
results1i_l9000 <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
                      q = 10, validation = 0.6, b = 5, b1_pos = TRUE, rh = 5,
                      family = "gaussian", seed = 2016, lambda = 9000)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

lambda_AIC_1i <- data.frame(lambda = c(0, 90, 900, 9000),
                            AIC = c(results1i$fit$aic, results1i_l90$fit$aic, 
                                    results1i_l900$fit$aic, results1i_l9000$fit$aic))
kable(lambda_AIC_1i) %>% kable_styling()


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

summary(results1i_l90)


## ---- results='asis', fig.show='hold', fig.height=8, fig.width=8, echo=TRUE----

gwqs_barplot(results1i_l90)


