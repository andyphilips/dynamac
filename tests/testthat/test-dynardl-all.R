context("test-dynardl-all")


test_that("Warnings are issued correctly", {
  model.notdyn <- lm(concern ~ urate, data = ineq)
  expect_error(dynardl.all.plots(model.notdyn), "provide a dynardl model object")
  
  model.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1), ec = FALSE,
                     simulate = FALSE)
  expect_error(dynardl.all.plots(model.1), "does not include simulation")  # no simulation
  
  model.4 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1, "concern" = 1), ec = FALSE,
                     simulate = TRUE, fullsims = FALSE, shockvar = "urate")  
  expect_error(dynardl.all.plots(model.4), "object must have fullsims")  # no fullsims provided
  
})