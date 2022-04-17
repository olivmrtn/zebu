## Test association measures

test_that(desc = "Chi-squared", {

  # compute chi-squared using zebu
  l <- zebu::lassie(trees, 1:2, 1:2, measure = "chisq")

  # compute mutual information using svs
  tab <- table(l$data$pp)
  r_chisq_test <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  r_chisq_res <- residuals(r_chisq_test)
  class(r_chisq_res) <- "matrix"

  # global association test
  expect_equal(l$global, as.numeric(r_chisq_test$statistic))

  # local association test
  expect_equal(l$local, r_chisq_res)
})

test_that(desc = "Mutual Information", {

  # compute chi-squared using zebu
  l <- zebu::lassie(trees, 1:2, 1:2, measure = "pmi")

  # compute chi-squared using stats::chisq.text
  tab <- table(l$data$pp)
  r_pmi <- svs::pmi(tab, normalize = FALSE, base = 2)
  class(r_pmi) <- "matrix"
  r_pmi[l$prob$observed == 0] <- -Inf
  r_mi <- svs::MI(tab)

  # global association test
  expect_equal(l$global, r_mi)

  # local association test
  expect_equal(l$local, r_pmi)
})

test_that(desc = "Normalized Mutual Information", {

  # compute chi-squared using zebu
  l <- zebu::lassie(trees, 1:2, 1:2, measure = "npmi")

  # compute normalized mutual information using svs
  tab <- table(l$data$pp)
  r_npmi <- svs::pmi(tab, normalize = TRUE, base = 2)
  class(r_npmi) <- "matrix"
  r_npmi[l$prob$observed == 0] <- -1

  prob <- prop.table(table(l$data$pp))
  r_nmi <- sum(prob * r_npmi)

  # global association test
  expect_equal(l$global, r_nmi)

  # local association test
  expect_equal(l$local, r_npmi)
})
