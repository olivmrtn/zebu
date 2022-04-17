## ----opts = TRUE, setup = TRUE, include = FALSE-------------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("zebu")

## -----------------------------------------------------------------------------
library(zebu) 

## -----------------------------------------------------------------------------
starter_prob <- c("Tomato Mozzarella Salad" = 1/3, "Rice Tuna Salad" = 1/3, "Lentil Salad" = 1/3)
starter_prob

## -----------------------------------------------------------------------------
main_given_starter_prob <- matrix(c(5/11, 1/11, 5/11,
                                    5/11, 5/11, 1/10,
                                    1/11, 5/11, 5/11),
                                  3, 3, byrow = TRUE)
rownames(main_given_starter_prob) <- names(starter_prob)
colnames(main_given_starter_prob) <- c("Sausage and Lentil Stew", "Pizza Margherita", "Pilaf Rice")
main_given_starter_prob

## -----------------------------------------------------------------------------
dessert_given_main <- matrix(c(2/6, 2/6, 2/6,
                               7/12, 1/12, 2/6, 
                               1/12, 7/12, 2/6),
                             3, 3, byrow = TRUE)
rownames(dessert_given_main) <- colnames(main_given_starter_prob)
colnames(dessert_given_main) <- c("Rice Pudding", "Apple Pie", "Fruit Salad")
dessert_given_main

## -----------------------------------------------------------------------------
set.seed(0)
sample_size <- 1000
df <- t(sapply(seq_len(sample_size), function(i) {
  
  starter <- sample(names(starter_prob), size = 1, prob = starter_prob)
  main <- sample(colnames(main_given_starter_prob), size = 1, prob = main_given_starter_prob[starter, ])
  dessert <- sample(colnames(dessert_given_main), size = 1, prob = dessert_given_main[main, ])
  
  c(Starter = starter, Main = main, Dessert = dessert)
}))
df <- as.data.frame(df)

## -----------------------------------------------------------------------------
head(df)

## -----------------------------------------------------------------------------
table(df)

## -----------------------------------------------------------------------------
las <- lassie(df, select = c("Main", "Dessert"), measure = "z")

## -----------------------------------------------------------------------------
las <- permtest(las, 
                nb = 5000, 
                p_adjust = "BH")

## ----plot-local-association---------------------------------------------------
print(las)
plot(las)

## -----------------------------------------------------------------------------
las2 <- lassie(df, measure = "z")
las2 <- permtest(las2, nb = 5000)
print(las2, what_sort = "local_p", decreasing = FALSE)

