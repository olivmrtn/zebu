## ----opts = TRUE, setup = TRUE, include = FALSE-------------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "")
library(ggplot2)

## ----include=FALSE, paged.print=FALSE-----------------------------------------
set.seed(63) # Set seed for reproducibility

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("zebu")

## ----eval=FALSE---------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("oliviermfmartin/zebu")

## -----------------------------------------------------------------------------
library(zebu) # Load zebu

## -----------------------------------------------------------------------------
data(trial) # Load trial dataset
head(trial) # Show head of trial dataset

## ----biomarker-histograms, echo=FALSE-----------------------------------------
ggplot(trial, aes(prebiom, fill = interaction(drug, resistance))) + 
  geom_histogram(alpha=0.5, position="identity", bins = 20) +
  xlab("Biomarker levels before treatment") +
  ylab("Number of Patients") +
  xlim(c(0, 1)) +
  scale_fill_discrete(name = "Patient", 
                      labels = c("Resistant and Drug", 
                                 "Resistant and Placebo",
                                 "Sensitive and Drug", 
                                 "Sensitive and Placebo"))

ggplot(trial, aes(postbiom, fill = interaction(drug, resistance))) + 
  geom_histogram(alpha=0.5, position="identity", bins = 20) +
  xlab("Biomarker levels after treatment") +
  ylab("Number of Patients") +
  xlim(c(0, 1)) +
  scale_fill_discrete(name = "Patient", 
                      labels = c("Resistant and Drug", 
                                 "Resistant and Placebo",
                                 "Sensitive and Drug", 
                                 "Sensitive and Placebo")) +
  geom_vline(xintercept = 0.7)

## -----------------------------------------------------------------------------
las <- lassie(trial, 
              select = c("drug", "postbiom"), 
              continuous = "postbiom", 
              breaks = c(0, 0.7, 1), 
              measure = "z")

## -----------------------------------------------------------------------------
las <- permtest(las, 
                nb = 1000, 
                p_adjust = "BH", 
                progress_bar = FALSE)

## ----plot-local-association---------------------------------------------------
print(las)
plot(las)

## -----------------------------------------------------------------------------
sub <- subgroups(las = las, 
                 x = trial, 
                 select = "resistance", 
                 thresholds = c(-0.01, 0.01),
                 significance = TRUE,
                 alpha = 0.01)

## -----------------------------------------------------------------------------
sub <- permtest(sub, nb = 1000)

## ----plot-subgroups-----------------------------------------------------------
print(sub)
plot(sub)

## -----------------------------------------------------------------------------
las2 <- lassie(trial, 
               select = c("drug", "postbiom", "resistance"), 
               continuous = "postbiom", 
               breaks = c(0, 0.7, 1))
las2 <- permtest(las2, 
                 group = list("drug", c("postbiom", "resistance")), progress_bar = FALSE)
print(las2)

