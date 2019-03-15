## ---- echo=FALSE, results='asis', message=FALSE--------------------------
library(gWQS)
library(Rsolnp)
library(ztable)
library(ggplot2)
library(tableHTML)
library(pander)
knitr::kable(head(wqs_data[, c(37, 36, 35, 1:34)], 10))

## ---- echo=FALSE, results='asis', message=FALSE--------------------------
knitr::kable(results$final_weights, digits = 3, row.names = F)

## ---- results='asis', message=FALSE, eval=F------------------------------
#  summary(results$fit)

## ---- echo=F, results='asis', message=FALSE------------------------------
library(pander)
pander(results$fit)

## ---- results='asis', message=FALSE, eval=F------------------------------
#  summary(results$fit_2)

## ---- echo=F, results='asis', message=FALSE------------------------------
pander(results$fit_2)

## ---- results='asis', message=FALSE, eval=FALSE--------------------------
#  results$aov

## ---- echo=F, results='asis', message=FALSE------------------------------
pander(results$aov)

## ---- results='asis', message=FALSE, eval=F------------------------------
#  summary(results2$fit)

## ---- echo=F, results='asis', message=FALSE------------------------------
pander(results2$fit)

