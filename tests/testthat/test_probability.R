test_that(desc = "Probabilities", {

  # compute probabilities using zebu
  l <- zebu::lassie(trees, 1:3, 1:3)

  # compute probabilities using R functions
  observed <- prop.table(table(l$data$pp))
  class(observed) <- class(l$prob$observed)

  margins <- lapply(seq_along(dim(observed)), function(i)
    apply(observed, i, sum))
  names(margins) <- colnames(l$data$pp)

  expected <- Reduce(outer, margins)
  dimnames(expected) <- dimnames(observed)

  # tests
  expect_equal(l$prob$observed, observed)
  expect_equal(l$prob$margins, margins)
  expect_equal(l$prob$expected, expected)

})
