## ----setup, include=FALSE-----------------------------------------------------
# knitr::opts_chunk$set(cache = TRUE)

library(gWQS)
library(ggplot2)
library(cowplot)
library(knitr)
library(kableExtra)

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE---------------------
#  
#  # we save the names of the mixture variables in the variable "PCBs"
#  PCBs <- names(wqs_data)[1:34]
#  # we run the model and save the results in the variable "results"
#  results <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data,
#                  q = 10, validation = 0.6, b = 2, b1_pos = TRUE,
#                  b1_constr = FALSE, family = "gaussian", seed = 2016)
#  # bar plot
#  gwqs_barplot(results)
#  # scatter plot y vs wqs
#  gwqs_scatterplot(results)
#  # scatter plot residuals vs fitted values
#  gwqs_fitted_vs_resid(results)
#  

## ----model1, results='asis', message=FALSE, warning=FALSE, fig.show='hold', fig.height=8, fig.width=8, echo=FALSE, fig.cap="Plots displayed for linear outcomes."----

PCBs <- names(wqs_data)[1:34]
results <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, 
                q = 10, validation = 0.6, b = 2, b1_pos = TRUE, 
                b1_constr = FALSE, family = "gaussian", seed = 2016)

w_ord <- order(results$final_weights$mean_weight)
mean_weight <- results$final_weights$mean_weight[w_ord]
mix_name <- factor(results$final_weights$mix_name[w_ord], 
                   levels = results$final_weights$mix_name[w_ord])
data_plot <- data.frame(mean_weight, mix_name)
bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight)) + 
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip() + ggtitle("A") +
  geom_hline(yintercept = 1/length(PCBs), linetype="dashed", color = "red")

yadj_vs_wqs <- ggplot(results$y_wqs_df, aes(wqs, y_adj)) + geom_point() + 
  stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw() + ggtitle("B")

fit_df <- data.frame(fitted = fitted(results), 
                     resid = residuals(results, type = "response"))
res_vs_fitted <- ggplot(fit_df, aes(x = fitted, y = resid)) + geom_point() + 
  theme_bw() + xlab("Fitted values") + ylab("Residuals") + ggtitle("C")

plot_grid(bar_plot_h, yadj_vs_wqs, res_vs_fitted, nrow = 2, ncol = 2)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

summary(results)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

head(results$final_weights)


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  gwqs_summary_tab(results)
#  

## ----sum1, results='asis', echo=FALSE-----------------------------------------

gwqs_summary_tab(results, caption = "Summary results of the WQS regression for linear outcomes.")


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  mf_df <- as.data.frame(signif(coef(summary(results$fit)), 3))
#  kable_styling(kable(mf_df, row.names = TRUE))
#  

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  gwqs_weights_tab(results)
#  

## ----w1, echo=FALSE, eval=TRUE------------------------------------------------

final_weight <- results$final_weights
final_weight[, -1] <- signif(final_weight[, -1], 3)
scroll_box(kable_styling(kable(final_weight, row.names = FALSE, caption = "Weights table of the WQS regression for linear outcomes.")), height = "400px")


## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  
#  final_weight <- results$final_weights
#  final_weight[, -1] <- signif(final_weight[, -1], 3)
#  kable_styling(kable(final_weight, row.names = FALSE))
#  

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  # bar plot
#  w_ord <- order(results$final_weights$mean_weight)
#  mean_weight <- results$final_weights$mean_weight[w_ord]
#  mix_name <- factor(results$final_weights$mix_name[w_ord],
#                     levels = results$final_weights$mix_name[w_ord])
#  data_plot <- data.frame(mean_weight, mix_name)
#  ggplot(data_plot, aes(x = mix_name, y = mean_weight)) +
#    geom_bar(stat = "identity", color = "black") + theme_bw() +
#    theme(axis.ticks = element_blank(),
#          axis.title = element_blank(),
#          axis.text.x = element_text(color='black'),
#          legend.position = "none") + coord_flip() +
#    geom_hline(yintercept = 1/length(PCBs), linetype="dashed", color = "red")
#  #
#  # scatter plot y vs wqs
#  ggplot(results$y_wqs_df, aes(wqs, y_adj)) + geom_point() +
#    stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()
#  #
#  # scatter plot residuals vs fitted values
#  fit_df <- data.frame(fitted = fitted(results),
#                       resid = residuals(results, type = "response"))
#  ggplot(fit_df, aes(x = fitted, y = resid)) + geom_point() +
#    theme_bw() + xlab("Fitted values") + ylab("Residuals")

