context("test-summary")

test_that("Summary is identical", {
  dyn.out <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    constant = FALSE, trend = FALSE, simulate = FALSE) 
  expect_identical(summary(dyn.out), summary(dyn.out$model))
})
