context("test-dynardl-simulation")


test_that("Warnings are issued correctly", {
  model.notdyn <- lm(concern ~ urate, data = ineq)
  expect_error(dynardl.simulation.plot(model.notdyn), "provide a dynardl model object")
  
  model.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1), ec = FALSE,
                     simulate = FALSE)
  expect_error(dynardl.simulation.plot(model.1), "does not include simulation")  # no simulation
  
  model.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1, "concern" = 1), ec = FALSE,
                     simulate = TRUE, shockvar = "urate")
  expect_error(dynardl.simulation.plot(model.2, response = "leveels"), "Response must be one of")  # not valid response 
  expect_error(dynardl.simulation.plot(model.2, type = "spikee"), "type must be either an area plot")  # not valid type 
  
  model.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1), ec = FALSE,
                     simulate = TRUE, noLDV = TRUE, shockvar = "urate")
  expect_warning(dynardl.simulation.plot(model.3), "from a model with no lagged dependent variable")  # no simulation
  
  expect_warning(dynardl.simulation.plot(model.2, last.period = 400), "requested exceeds simulation range")  # too long of response
  expect_warning(dynardl.simulation.plot(model.2, last.period = 400), "regardless of if Y might")  # too long of response
  
  model.4 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1, "concern" = 1), ec = FALSE,
                     simulate = TRUE, fullsims = FALSE, shockvar = "urate")  
  expect_error(dynardl.simulation.plot(model.4, response = "cumulative.diffs"), "object must have fullsims")  # no fullsims provided
  expect_error(dynardl.simulation.plot(model.4, response = "cumulative.abs.diffs"), "object must have fullsims")  # no fullsims provided
  
  model.5 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
                     lags = list("incshare10" = 1, "urate" = 1, "concern" = 1), ec = FALSE,
                     simulate = TRUE, fullsims = TRUE, shockvar = "urate")    
  expect_warning(dynardl.simulation.plot(model.5, response = "cumulative.abs.diffs", last.period = 12), "last.period")  # finish calc. on 12 
  expect_warning(dynardl.simulation.plot(model.5, response = "cumulative.abs.diffs"), "by tolerance")  # finish whenever

})
