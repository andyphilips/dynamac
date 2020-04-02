context("test-dynardl")

test_that("Correct number of parameters estimated: intercept vs. no intercept", {
  dyn.out.1.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)
  dyn.out.1.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, constant = FALSE)
  expect_equal(length(coef(dyn.out.1.2$model)), 3)
  expect_equal(length(coef(dyn.out.1.1$model)), 4)
  expect_equal(("(Intercept)" %in% names(coefficients(dyn.out.1.2))), FALSE)
})

test_that("Correct number of parameters estimated: trend vs. no trend", {
  dyn.out.2.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, trend = TRUE)
  dyn.out.2.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, trend = FALSE)
  expect_equal(length(coef(dyn.out.2.1$model)), 5)
  expect_equal(length(coef(dyn.out.2.2$model)), 4)
  expect_equal(("trendvar" %in% names(coefficients(dyn.out.2.1$model))), TRUE)
})

test_that("Correct number of parameters estimated: multiple lags", {
  dyn.out.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = c(1, 2, 3), "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)
  expect_equal(length(coef(dyn.out.3$model)), 6)
})

test_that("Correct length of simulation: default burnin", {
  dyn.out.4 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = c(1, 2, 3), "incshare10" = 1, "urate" = 1), ec = TRUE,
    shockvar = "urate", simulate = TRUE, range = 30)
  expect_equal(length(dyn.out.4$simulation$central), 30)
})

test_that("Correct length of simulation: custom burnin", {
  dyn.out.5 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = c(1, 2, 3), "incshare10" = 1, "urate" = 1), ec = TRUE,
    shockvar = "urate", simulate = TRUE, burnin = 20, range = 30)
  expect_equal(length(dyn.out.5$simulation$central), 30)
})

test_that("Correct number of bounds for custom significance levels", {
  dyn.out.6 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = c(1, 2, 3), "incshare10" = 1, "urate" = 1), ec = TRUE,
    shockvar = "urate", simulate = TRUE, burnin = 20, range = 30, sig = 80)
  expect_equal(ncol(dyn.out.6$simulation), 20)
})

test_that("Correctly simulating and including rawsims", {
  dyn.out.8.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    shockvar = "urate", simulate = TRUE, burnin = 20, range = 30, sig = 80, fullsims = FALSE)
  dyn.out.8.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    shockvar = "urate", simulate = TRUE, burnin = 20, fullsims = TRUE)
  dyn.out.8.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    shockvar = "urate", simulate = FALSE, burnin = 20, range = 30, sig = 80, fullsims = FALSE)
  expect_equal(dyn.out.8.1$rawsims, NULL)
  expect_equal(ncol(dyn.out.8.2$rawsims), 21) # 20 time periods + central tendency
  expect_equal(dyn.out.8.3$simulation, NULL)
})

test_that("Not including an LDV, if requested", {
  dyn.out.9.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, noLDV = TRUE)
  dyn.out.9.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, noLDV = FALSE)
  expect_equal(("l.1.concern" %in% names(coefficients(dyn.out.9.1$model))), FALSE)
  expect_equal(("l.1.concern" %in% names(coefficients(dyn.out.9.2$model))), TRUE)
})

test_that("Error correction is differencing Y", {
  dyn.out.10.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)  
  dyn.out.10.2 <- lm(dshift(ineq$concern) ~ lshift(ineq$concern, 1) + lshift(ineq$incshare10, 1) + lshift(ineq$urate, 1))
  expect_equal((coefficients(dyn.out.10.1$model) %in% coefficients(dyn.out.10.2)), rep(TRUE, length(coefficients(dyn.out.10.1$model))))
})

test_that("DV named correctly (for pssbounds)", {
	dyn.out.10.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)  
  expect_equal(dyn.out.10.1$model$y.name, "concern")	
})

test_that("Warnings are issued correctly", {
  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    diffs = list("incshare" = 1),
    simulate = FALSE), "vector of variables to be") # error: diffs can't be a list
  
  expect_warning(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    diffs = c("urate"),
    simulate = FALSE), "the diffs list not found in model") # warning: diffs not in formula

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = c("incshare10" = 1, "urate" = 1), ec = TRUE,
    diffs = c("urate"),
    simulate = FALSE), "must be a list") # error: lags must be a list
 
  expect_warning(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    diffs = c("urate"), noLDV = TRUE,
    simulate = FALSE), "Are you sure you want this") # warning: noLDV bad idea
  
  expect_warning(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE), "lags list not found in model formula") # warning: lags not in formula

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    lagdiffs = c("urate"),
    simulate = FALSE), "must be a list") # error: lagdiffs must be a list
  
  expect_warning(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    lagdiffs = list("urate" = 1),
    simulate = FALSE), "lagdiffs list not found in model") # warning: lagdiffs not in formula

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1),
    levels = list("urate"),
    simulate = FALSE), "levels as a vector of variables") # error: levels can't be a list

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    levels = c("concern"),
    simulate = FALSE), "variable cannot appear in contemporaneous") # error: DV can't be in levels

  expect_warning(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1),
    levels = c("urate"),
    simulate = FALSE), "levels list not found in model") # warning: levels not in formula
  
  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = TRUE, forceset = list("concern" = 1)), "LDV cannot be forceset") # error: can't force DV
  
  expect_warning(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE, 
    shockvar = "incshare10", sims = 10, range = 10, time = 5, burnin = 1,
    simulate = TRUE, forceset = list("urate" = 1)), "the forceset list not found") # warning: forceset not in formula 
  
  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, shockvar = incshare10), "must be specified in quotation") # error: shockvar not in quotes 
    
  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE, shockvar = c("incshare10", "urate")), 
    	"specified more than one shockvar") # warning: more than one shockvar, but not simulating
    
  expect_warning(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    simulate = FALSE, shockvar = c("urate")),
 			"is not in the model formula, but since you are not simulating") # warning: shockvar not found in formula, but not simulating
       
  expect_warning(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    simulate = FALSE, shockvar = c("urate")), 
 			"shockvar is not found in lags") # warning: shockvar in formula but not in model, but not simulating

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = TRUE), "must be specified for a dynamic") # error: no shockvar 

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = TRUE, shockvar = incshare10), "must be specified in quotation") # error: shockvar not in quotes 
    
  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = TRUE, shockvar = c("incshare10", "urate")), "specify one shockvar") # error: more than one shockvar
    
  expect_error(dynardl(concern ~ incshare10, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    simulate = TRUE, shockvar = c("urate")), "not found in model formula") # error: shockvar not found in formula
       
  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    simulate = TRUE, shockvar = c("urate")), "Shock variable required for dynamic") # error: shockvar in formula but not in model

  expect_output(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1, "concern" = 1), ec = TRUE, sims = 10, range = 10, time = 5, burnin = 1,
    simulate = TRUE, shockvar = c("urate")), "shocked by one standard deviation") # message: shockvar shocked by 1 sd

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE, range = 10, time = 30,
    simulate = TRUE, shockvar = c("urate")), "range of simulation must be longer") # error: shocktime greater than range  

  expect_warning(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1), ec = TRUE,
    simulate = FALSE), "in the formula list not found lags") # message: urate isn't anywhere in data
 
  expect_output(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE), "dependent variable to be run in differences") # message: EC message
 
  expect_output(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = FALSE,
    simulate = FALSE), "Dependent variable to be run in levels") # message: levels message
 
  expect_warning(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE, noLDV = FALSE,
    simulate = FALSE), "Lagged dependent variable added to model formula") # message: adding variable (LDV not in list)
 
  expect_warning(dynardl(concern ~ incshare10 + urate, data = ineq, 
    diffs = c("incshare10", "urate"), ec = FALSE, noLDV = FALSE,
    simulate = FALSE), "Lagged dependent variable added to model formula") # message: adding variable (no lags list)

  expect_error(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 2, "incshare10" = 1, "urate" = 1), ec = TRUE, noLDV = FALSE,
    simulate = FALSE), "must include first lag of dependent") # stop: no LDV variable (LDV lag 2, not 1)

  expect_output(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = FALSE, trend = TRUE,
    simulate = FALSE), "Deterministic linear trend added to model formula") # message: constant message

  expect_output(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = FALSE, constant = FALSE,
    simulate = FALSE), "Constant suppressed from model formula") # message: constant message

  expect_warning(dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = c(1, 2, 3), "urate" = 1), ec = TRUE, noLDV = FALSE,
    simulate = FALSE), "lags included in implied") # message: multiple lags in cointegrated

})

