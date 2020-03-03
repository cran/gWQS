### R code from vignette source 'gwqs-vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: gwqs-vignette.Rnw:91-95
###################################################
library(knitr)
opts_chunk$set(
engine='R', tidy=FALSE, concordance=FALSE
)


###################################################
### code chunk number 2: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library(gWQS)
library(ggplot2)
library(broom)
library(gridExtra)
library(plotROC)
library(reshape2)
library(ggrepel)
library(VGAM)


###################################################
### code chunk number 3: linear_reg
###################################################
# we save the names of the mixture variables in the variable "toxic_chems"
toxic_chems <- names(wqs_data)[1:34]
# we run the model and save the results in the variable "results"
results <- gwqs(y ~ wqs, mix_name = toxic_chems,
               data = wqs_data, q = 4, validation = 0.6, b = 2,
               b1_pos = TRUE, b1_constr = FALSE, family = "gaussian",
               seed = 2016, plots = TRUE, tables = TRUE)
#
# bar plot
w_ord <- order(results$final_weights$mean_weight)
mean_weight <- results$final_weights$mean_weight[w_ord]
mix_name <- factor(results$final_weights$mix_name[w_ord],
                   levels = results$final_weights$mix_name[w_ord])
data_plot <- data.frame(mean_weight, mix_name)
ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip()
#
# scatter plot y vs wqs
ggplot(results$y_wqs_df, aes(wqs, y_adj)) + geom_point() +
  stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw()
#
# scatter plot residuals vs fitted values
fit_df <- broom::augment(results$fit)
ggplot(fit_df, aes(x = .fitted, y = .resid)) + geom_point() +
  theme_bw() + xlab("Fitted values") + ylab("Residuals")



###################################################
### code chunk number 4: fig1
###################################################
bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip() + ggtitle("A")


yadj_vs_wqs <- ggplot(results$y_wqs_df, aes(wqs, y_adj)) + geom_point() +
  stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw() + ggtitle("B")


res_vs_fitted <- ggplot(fit_df, aes(x = .fitted, y = .resid)) + geom_point() +
  theme_bw() + xlab("Fitted values") + ylab("Residuals") + ggtitle("C")


grid.arrange(bar_plot_h, yadj_vs_wqs, res_vs_fitted, ncol=2)


###################################################
### code chunk number 5: linear_reg_sum
###################################################
summary(results$fit)


###################################################
### code chunk number 6: linear_reg_weights
###################################################
head(results$final_weights)


###################################################
### code chunk number 7: logistic_reg
###################################################
# we run the logistic model and save the results in the variable
# "results2"
results2 <- gwqs(y_bin ~ wqs + sex, mix_name = toxic_chems,
                data = wqs_data, q = NULL, validation = 0.4, b = 2,
                b1_pos = TRUE, b1_constr = FALSE, family = binomial,
                seed = 2018, plots = TRUE, tables = FALSE, pred = 0.3)
#
# plot ROC curve
gg_roc <- ggplot(results2$df_pred, aes(d=y, m=p_y)) + geom_roc(n.cuts = 0) +
  style_roc(xlab = "1 - Specificity", ylab = "Sensitivity")
auc_est <- plotROC::calc_auc(gg_roc)
gg_roc + annotate("text", x=0.75, y=0.25,
                  label=paste0("AUC = ", round(auc_est[, "AUC"], 3)))



###################################################
### code chunk number 8: fig2
###################################################
w_ord <- order(results2$final_weights$mean_weight)
mean_weight <- results2$final_weights$mean_weight[w_ord]
mix_name <- factor(results2$final_weights$mix_name[w_ord], levels = results2$final_weights$mix_name[w_ord])
data_plot <- data.frame(mean_weight, mix_name)
bar_plot_h <- ggplot(results2$final_weights, aes(x = mix_name, y = mean_weight, fill = mix_name))
bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip() + ggtitle("A")


yadj_vs_wqs <- ggplot(results2$y_wqs_df, aes(wqs, y)) + geom_point() +
  stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw() + ggtitle("B")


gg_roc <- gg_roc + annotate("text", x=0.75, y=0.25, label=paste0("AUC = ", round(auc_est[, "AUC"], 3))) +
  ggtitle("C")


grid.arrange(bar_plot_h, yadj_vs_wqs, gg_roc, ncol=2)


###################################################
### code chunk number 9: logistic_reg_sum
###################################################
summary(results2$fit)


###################################################
### code chunk number 10: logistic_reg_pred
###################################################
# create a dataset exluding the data where we want to apply the prediction
# and define the group variable to identify the test and validation dataset
wqs_data$group <- 0
wqs_data$group[results2$vindex] <- 1
wqs_data_train <- wqs_data[-results2$pindex,]
# fit the model on the training dataset
results2_pred <- gwqs(y_bin ~ wqs + sex, mix_name = toxic_chems,
                data = wqs_data_train, q = NULL, validation = NULL,
                b = 2, valid_var = "group", b1_pos = TRUE,
                b1_constr = FALSE, family = binomial, seed = 2018)
# creat the dataset on which we apply the prediction
wqs_data_pred <- wqs_data[results2$pindex,]
# create wqs variable for the prediction dataset
mix_matrx <- as.matrix(wqs_data_pred[, rownames(results2$final_weights)])
wqs_data_pred$wqs <- as.numeric(mix_matrx%*%results2$final_weights$mean_weight)
# apply the predict() function
pred <- predict(results2$fit, newdata = wqs_data_pred, type = "response")
df_pred <- data.frame(y = wqs_data_pred$y_bin, p_y = pred)
# plot the roc curve
gg_roc <- ggplot(df_pred, aes(d=y, m=p_y)) + geom_roc(n.cuts = 0) +
  style_roc(xlab = "1 - Specificity", ylab = "Sensitivity")
auc_est <- plotROC::calc_auc(gg_roc)
gg_roc + annotate("text", x=0.75, y=0.25,
                  label=paste0("AUC = ", round(auc_est[, "AUC"], 3)))


