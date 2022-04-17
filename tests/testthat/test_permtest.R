skip_on_cran()

test_that(desc = "Permutation test", {

  res <- pbapply::pblapply(seq_len(10), function(trial) {

    # simulate date
    data <- sapply(seq_len(2), function(i) sample(0:3, 200, TRUE))

    # chisq using lassie
    l <- lassie(data, measure = "chisq")
    l <- permtest(l, 5000, p_adjust = "none")

    # chisq using R
    chisq <- suppressWarnings(chisq.test(table(data[, 1], data[, 2])))

    list(
      lassie_global_p = l$global_p,
      r_global_p = as.numeric(chisq$p.value),
      lassie_local_p = l$local_p,
      r_local_p = as.numeric(2*(1- pnorm(abs(chisq$stdres)))))
  })

  # retrieve data
  lassie_global_p <- sapply(res, function(x) as.numeric(x$lassie_global_p))
  r_global_p <- sapply(res, function(x) as.numeric(x$r_global_p))

  lassie_local_p <- sapply(res, function(x) as.numeric(x$lassie_local_p))
  r_local_p <- sapply(res, function(x) as.numeric(x$r_local_p))

  # compute correlation between R and lassie
  rho_global_p <- cor(lassie_global_p, r_global_p)

  rho_local_p <- sapply(seq_len(ncol(lassie_local_p)), function(i)
    cor(lassie_local_p[, i], r_local_p[, i]))

  # test
  expect_true(all(rho_global_p > 0.95))
  expect_true(all(rho_local_p > 0.90))
})
