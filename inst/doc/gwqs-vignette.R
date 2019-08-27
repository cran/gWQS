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


###################################################
### code chunk number 11: multinom_reg
###################################################
# we create the variable "group" in the dataset to identify the training
# and validation dataset: we choose 300 observations for the validation
# dataset and the remaining 200 for the training dataset
set.seed(123)
wqs_data$group <- 0
wqs_data$group[rownames(wqs_data) %in%
                 sample(rownames(wqs_data), 300)] <- 1
#
# we run the logistic model and save the results in the variable
# "results3"
results3 <- gwqs(y_multinom ~ wqs, mix_name = toxic_chems,
                data = wqs_data, q = NULL, validation = 0.6,
                valid_var = "group", b = 2, b1_pos = c(TRUE, TRUE),
                b1_constr = FALSE, family = "multinomial", seed = 123,
                plots = TRUE, tables = TRUE,
                plan_strategy = "multiprocess")
#
# bar plot
data_plot <- results3$final_weights[order(results3$final_weights[, 2]),]
pos <- match(data_plot$mix_name, sort(data_plot$mix_name))
data_plot$mix_name <- factor(data_plot$mix_name,
                            levels(data_plot$mix_name)[pos])
data_plot_l <- melt(data_plot, id.vars = "mix_name")
ggplot(data_plot_l, aes(x = mix_name, y = value, fill = mix_name)) +
  facet_wrap(~ variable) + geom_bar(stat = "identity", color = "black") +
  theme_bw() + theme(axis.ticks = element_blank(),
                     axis.title = element_blank(),
                     axis.text.x = element_text(color='black'),
                     legend.position = "none") + coord_flip()
#
# scatter plot y vs wqs
ggplot(results3$y_wqs_df, aes(wqs, y)) +
  geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) +
  theme_bw() + facet_wrap(~ level)
#
# scatter plot of weights for the two levels of the dependent variable
ggplot(data_plot, aes_string(names(data_plot)[2], names(data_plot)[3])) +
  geom_point() + theme_bw() + xlab(names(data_plot)[2]) +
  ylab(names(data_plot)[3]) + geom_abline(linetype = 2) +
  ggrepel::geom_text_repel(aes(label=mix_name))



###################################################
### code chunk number 12: fig3
###################################################
bar_plot_h <- ggplot(data_plot_l, aes(x = mix_name, y = value, fill = mix_name)) +
  facet_wrap(~ variable) + geom_bar(stat = "identity", color = "black") +
  theme_bw() + theme(axis.ticks = element_blank(),
                     axis.title = element_blank(),
                     axis.text.x = element_text(color='black'),
                     legend.position = "none") + coord_flip() + ggtitle("A")


yadj_vs_wqs <- ggplot(results3$y_wqs_df, aes(wqs, y)) +
  geom_point() + stat_smooth(method = "loess", se = FALSE, size = 1.5) +
  theme_bw() + facet_wrap(~ level) + ggtitle("B")

grid.arrange(bar_plot_h, yadj_vs_wqs, ncol=2)


###################################################
### code chunk number 13: fig3b
###################################################
w1_vs_w2 <- ggplot(data_plot, aes_string(names(data_plot)[2], names(data_plot)[3])) +
  geom_point() + theme_bw() + xlab(names(data_plot)[2]) +
  ylab(names(data_plot)[3]) + geom_abline(linetype = 2) +
  ggrepel::geom_text_repel(aes(label=mix_name)) + ggtitle("C")

w1_vs_w2


###################################################
### code chunk number 14: multinom_reg_sum
###################################################
results3$fit$sum_stat


###################################################
### code chunk number 15: poisson_reg
###################################################
# we create the sex factor variable sex_factor
wqs_data$sex_factor <- factor(wqs_data$sex, labels = c("F", "M"))
#
# we run the poisson model and save the results in the variable
# "results4"
results4 <- gwqs(y_count ~ wqs, mix_name = toxic_chems,
                stratified = "sex_factor", data = wqs_data, q = 10,
                validation = 0.6, b = 2, b1_pos = TRUE,
                b1_constr = FALSE, family = poisson, seed = 123,
                plots = TRUE, tables = TRUE)



###################################################
### code chunk number 16: fig4
###################################################

w_ord <- order(results4$final_weights$mean_weight)
mean_weight <- results4$final_weights$mean_weight[w_ord]
mix_name <- factor(results4$final_weights$mix_name[w_ord], levels = results4$final_weights$mix_name[w_ord])
data_plot <- data.frame(mean_weight, mix_name)
bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip() + ggtitle("A")


yadj_vs_wqs <- ggplot(results4$y_wqs_df, aes(wqs, y_adj)) + geom_point() +
  stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw() + ggtitle("B")

fit_df <- augment(results4$fit)
res_vs_fitted <- ggplot(fit_df, aes(x = .fitted, y = .resid)) + geom_point() +
  theme_bw() + xlab("Fitted values") + ylab("Residuals") + ggtitle("C")

grid.arrange(bar_plot_h, yadj_vs_wqs, res_vs_fitted, ncol=2, nrow = 2, layout_matrix = cbind(c(1,1), c(2,3)))


###################################################
### code chunk number 17: poisson_reg_sum
###################################################
summary(results4$fit)


###################################################
### code chunk number 18: poisson_disp_test
###################################################
library(AER)
mean(wqs_data$y_count)
var(wqs_data$y_count)
AER::dispersiontest(results4$fit)


###################################################
### code chunk number 19: quasi_poisson_reg
###################################################
# we run the quasi-poisson model and save the results in the variable
# "results5"
results5 <- gwqs(y_count ~ wqs, mix_name = toxic_chems,
                data = wqs_data, q = 10, validation = 0.6, b = 2,
                b1_pos = TRUE, b1_constr = FALSE, family = quasipoisson,
                seed = 123)


###################################################
### code chunk number 20: quasi_poisson_reg_sum
###################################################
summary(results5$fit)


###################################################
### code chunk number 21: nb_reg
###################################################
# generate new variable from normal distribution
set.seed(123)
wqs_data$new_var <- rnorm(500)
wqs_data$y_zinb <- rzinegbin(500, pstr0 = 0.3, mu = 3, size = 10)
#
# we run the zero-inflated negative binomial model and save the results in the variable
# "results6"
results6 <- gwqs(y_zinb ~ wqs + sex | new_var, mix_name = toxic_chems,
                data = wqs_data, q = 10, validation = 0.6, b = 2,
                zero_infl = TRUE, zilink = "logit", b1_pos = FALSE,
                b1_constr = FALSE, family = "negbin", seed = 1234,
                plots = TRUE, tables = TRUE)
#
# scatter plot residuals vs fitted values
fit_df <- data.frame(.fitted = results6$fit$fitted.values,
                    .resid = results6$fit$residuals)
ggplot(fit_df, aes(x = .fitted, y = .resid)) + geom_point() +
  theme_bw() + xlab("Fitted values") + ylab("Residuals")



###################################################
### code chunk number 22: fig5
###################################################
w_ord <- order(results6$final_weights$mean_weight)
mean_weight <- results6$final_weights$mean_weight[w_ord]
mix_name <- factor(results6$final_weights$mix_name[w_ord], levels = results6$final_weights$mix_name[w_ord])
data_plot <- data.frame(mean_weight, mix_name)
bar_plot_h <- ggplot(data_plot, aes(x = mix_name, y = mean_weight, fill = mix_name)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black'),
        legend.position = "none") + coord_flip() + ggtitle("A")

yadj_vs_wqs <- ggplot(results6$y_wqs_df, aes(wqs, y_adj)) + geom_point() +
  stat_smooth(method = "loess", se = FALSE, size = 1.5) + theme_bw() + ggtitle("B")

fit_df <- data.frame(.fitted = results6$fit$fitted.values, .resid = results6$fit$residuals)
res_vs_fitted <- ggplot(fit_df, aes(x = .fitted, y = .resid)) + geom_point() +
  theme_bw() + xlab("Fitted values") + ylab("Residuals") + ggtitle("C")

grid.arrange(bar_plot_h, yadj_vs_wqs, res_vs_fitted, ncol=2)


###################################################
### code chunk number 23: nb_reg_sum
###################################################
summary(results6$fit)


