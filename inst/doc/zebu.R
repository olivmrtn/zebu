## ----opts = TRUE, setup = TRUE, include = FALSE--------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "")

## ------------------------------------------------------------------------
set.seed(63) # Set seed for reproducibility

library(zebu) # Load zebu
data(trial) # Load trial dataset
head(trial) # Show head of trial dataset

## ------------------------------------------------------------------------
las <- lassie(trial, 
              select = c("drug", "postbiom"), 
              continuous = "postbiom", 
              breaks = c(0, 0.7, 1), 
              measure = "z")

## ------------------------------------------------------------------------
las <- permtest(las, 
                nb = 1000, 
                p_adjust = "BH", 
                parallel = FALSE, 
                progress_bar = FALSE)

## ------------------------------------------------------------------------
print(las)
plot(las)

## ------------------------------------------------------------------------
sub <- subgroups(las = las, 
                 x = trial, 
                 select = "resistance", 
                 thresholds = c(-0.05, 0.05),
                 significance = TRUE,
                 alpha = 0.01)

## ------------------------------------------------------------------------
sub <- permtest(sub)

## ------------------------------------------------------------------------
print(sub)
plot(sub)

## ------------------------------------------------------------------------
las2 <- lassie(trial, 
               select = c("drug", "postbiom", "resistance"), 
               continuous = "postbiom", 
               breaks = c(0, 0.7, 1))
las2 <- permtest(las2, 
                 group = list("drug", c("postbiom", "resistance")), progress_bar = FALSE)
print(las2)

