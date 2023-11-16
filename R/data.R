#' Exposure concentrations of 34 PCB (simulated dataset)
#'
#' We created the `wqs_data` dataset to show how to use this function. These data reflect
#' 59 exposure concentrations simulated from a distribution of 34 PCB exposures and
#' 25 phthalate biomarkers measured in subjects participating in the NHANES study (2001-2002).
#' Additionally, 8 outcome measures were simulated applying different distributions and fixed
#' beta coefficients to the predictors. In particular `y` and `yLBX` were simulated from
#' a normal distribution, `ybin` and `ybinLBX` from a binomial distribution, `ymultinom` and
#' `ymultinomLBX` from a multinomial distribution and `ycount` and `ycountLBX` from a Poisson
#' distribution.
#' The regression coefficients used to generate the outcomes `yLBX`, `ybinLBX` and
#' `ycountLBX` were set to:\cr
#' LBX105LA = 0.3\cr
#' LBX138LA = 0.6\cr
#' LBX157LA = 0.2\cr
#' LBXD02LA = 0.45\cr
#' LBXD04LA = 0.15\cr
#' LBXF06LA = 0.3\cr
#' LBXF07LA = 0.45\cr
#' then the following terms were added to generate the variables `y`, `ybin` and
#' `ycount`:\cr
#' URXMC1 = 0.15\cr
#' URXMOH = 0.45\cr
#' URXP02 = 0.2\cr
#' URXP10 = 0.3\cr
#' URXUCR = 0.2\cr
#' All the remaining coefficients were set to 0.\cr
#' The coefficients to generate `ymultinomLBX` were set as below:\cr
#' level B:\cr
#' LBX138LA = 0.8\cr
#' LBXD04LA = 0.2\cr
#' level C:\cr
#' LBX105LA = 0.4\cr
#' LBX157LA = 0.3\cr
#' LBXD02LA = 0.6\cr
#' LBXF06LA = 0.4\cr
#' LBXF07LA = 0.6\cr
#' and the following terms were added for `ymultinom`:\cr
#' level B:\cr
#' URXMC1 = 0.2\cr
#' URXP02 = 0.3\cr
#' URXP10 = 0.4\cr
#' URXUCR = 0.3\cr
#' level C:\cr
#' URXMOH = 0.6\cr
#' The `sex` variable was also simulated to allow to adjust for a covariate in the model.
#' This dataset can thus be used to test the `gWQS` package by analyzing the mixed effect
#' of the 59 simulated PCBs on the continuous, binary or count outcomes, with adjustments
#' for covariates.
#'
#' \describe{
#' \item{y}{continuous outcome generated considerig all the predictors}
#' \item{yLBX}{continuous outcome generated considerig only PCBs}
#' \item{ybin}{binary outcome generated considerig all the predictors}
#' \item{ybinLBX}{binary outcome generated considerig only PCBs}
#' \item{ymultinom}{multinomial outcome generated considerig all the predictors}
#' \item{ymultinomLBX}{multinomial outcome generated considerig only PCBs}
#' \item{ycount}{count outcome generated considerig all the predictors}
#' \item{ycountLBX}{count outcome generated considerig only PCBs}
#' \item{sex}{covariate, gender of the subject}
#' \item{LBX}{34 exposure concentrations of PCB}
#' \item{URX}{25 exposure concentrations of phthalates}
#' ...
#' }
#'
#' @format A data frame with 500 rows and 68 variables
"wqs_data"



#' Measurement of 38 nutrients (NHANES dataset)
#'
#' We created the `tiwqs_data` dataset to show how apply the two-indices WQS (2iWQS).
#' These data reflect 38 nutrients measured in 5960 subjects participating in the NHANES study
#' (2011-2012, 2013-2014, 2015-2016 cycles). Nutrients were estimated from the dietary intake
#' data that considered the types and amounts of foods and beverages consumed during the 24-hour
#' period prior to the interview. Two interviews were performed: the first one was collected
#' in-person while the second interview was collected by telephone 3 to 10 days later. In this
#' study we averaged the two nutrients when both evaluations were considered as usual food
#' consumption (only one measurement was included in the analysis if the other one was not usual
#' while the observation was dropped if both evaluations were not usual) and we added the dietary
#' supplement intake when applicable.
#' Additionally, BMI as both a continuous variable (`bmxbmi`) and categorical variable (`bmi_cat`)
#' was included as outcome variable.
#' A total of 10 covariates can also be found in the dataset
#'
#' @docType data
#'
#' @usage data(tiwqs_data)
#'
#' @format A data frame with 5960 rows and 50 variables
"tiwqs_data"


