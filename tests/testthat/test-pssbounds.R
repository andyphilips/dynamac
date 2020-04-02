context("test-pssbounds")


test_that("Correct output for pssbounds", {
  dyn.out.7.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    constant = FALSE, trend = FALSE, simulate = FALSE) # case 1
  pss.7.1 <- pssbounds(dyn.out.7.1, object.out = TRUE)

  dyn.out.7.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    constant = TRUE, trend = FALSE, simulate = FALSE) # case 2
  pss.7.2 <- pssbounds(dyn.out.7.2, restriction = TRUE, object.out = TRUE)

  dyn.out.7.3 <- dynardl(concern ~ incshare10, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1), ec = TRUE,
    constant = TRUE, trend = FALSE, simulate = FALSE) # case 3
  pss.7.3 <- pssbounds(dyn.out.7.3, object.out = TRUE)

  dyn.out.7.4 <- dynardl(concern ~ incshare10, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1), ec = TRUE,
    constant = TRUE, trend = TRUE, simulate = FALSE) # case 4
  pss.7.4 <- pssbounds(dyn.out.7.4, restriction = TRUE, object.out = TRUE)

  dyn.out.7.5 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = c(1, 2, 3), "urate" = 1), ec = TRUE,
    constant = TRUE, trend = TRUE, simulate = FALSE) # case 5
  pss.7.5 <- pssbounds(dyn.out.7.5, object.out = TRUE)

  expect_equal(pss.7.1$case, 1)
  expect_equal(pss.7.2$case, 2)
  expect_equal(pss.7.3$case, 3)
  expect_equal(pss.7.4$case, 4)
  expect_equal(pss.7.5$case, 5)
  expect_equal(pss.7.1$k, 2)
  expect_equal(pss.7.2$k, 2)
  expect_equal(pss.7.3$k, 1)
  expect_equal(pss.7.4$k, 1)
  expect_equal(pss.7.5$k, 2)
  expect_equal(round(coef(summary(dyn.out.7.1))["l.1.concern", 3], digits = 3), pss.7.1$tstat)
  expect_equal(round(coef(summary(dyn.out.7.2))["l.1.concern", 3], digits = 3), pss.7.2$tstat)
  expect_equal(round(coef(summary(dyn.out.7.3))["l.1.concern", 3], digits = 3), pss.7.3$tstat)
  expect_equal(round(coef(summary(dyn.out.7.4))["l.1.concern", 3], digits = 3), pss.7.4$tstat)
  expect_equal(round(coef(summary(dyn.out.7.5))["l.1.concern", 3], digits = 3), pss.7.5$tstat)
  expect_equal(pss.7.1$obs, length(dyn.out.7.1$model$residuals))
  expect_equal(pss.7.2$obs, length(dyn.out.7.2$model$residuals))
  expect_equal(pss.7.3$obs, length(dyn.out.7.3$model$residuals))
  expect_equal(pss.7.4$obs, length(dyn.out.7.4$model$residuals))
  expect_equal(pss.7.5$obs, length(dyn.out.7.5$model$residuals))
  # fstat should be equal to fstat from regression with just vars. in levels
  expect_equal(unname(round(summary(lm(dshift(ineq$concern) ~ lshift(ineq$concern, 1) + lshift(ineq$incshare10, 1)))$fstatistic[1], 3)), 
  		round(pss.7.3$fstat, 3))
  expect_output(pssbounds(dyn.out.7.1), "No intercept; no trend")
  expect_output(pssbounds(dyn.out.7.2, restriction = TRUE), "Intercept included in F-stat restriction; no trend")
  expect_output(pssbounds(dyn.out.7.3), "Unrestricted intercept; no trend")
  expect_output(pssbounds(dyn.out.7.4, restriction = TRUE), "Unrestricted intercept; trend included in F-stat restriction")
  expect_output(pssbounds(dyn.out.7.5), "Unrestricted intercept; unrestricted trend")
})


test_that("Warnings are issued correctly", {
	
  model.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = FALSE,
    simulate = FALSE)
  expect_error(pssbounds(model.1), "only for error-correcting relationships")  # not ECM

  model.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lagdiffs = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)
  expect_error(pssbounds(model.2), "no variables in first lags") # k = 0

  model.3 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lagdiffs = list("incshare10" = 1, "urate" = 1), ec = TRUE, noLDV = TRUE,
    simulate = FALSE)
  expect_error(pssbounds(model.3), "only for models that include a lagged dependent variable") # no LDV

  model.4 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)
  expect_warning(pssbounds(model.4, obs = 4), "different from observations calculated") # supplied wrong obs + model
  expect_warning(pssbounds(model.4, fstat = 4), "different from F") # supplied wrong fstat + model
  expect_warning(pssbounds(model.4, tstat = 4), "different from t") # supplied wrong tstat + model
  expect_warning(pssbounds(model.4, case = 4), "different from case calculated") # supplied wrong case + model
  expect_warning(pssbounds(model.4, k = 400), "different from k calculated") # supplied wrong k + model

  expect_error(pssbounds(tstat = 6), "Provide either number of observations")
  expect_error(pssbounds(obs = 6), "Provide either fstat on lagged variables")
  expect_error(pssbounds(obs = 6, fstat = 5), "Provide either case of regression")
  expect_error(pssbounds(obs = 6, fstat = 5, case = 2), "Provide either k number of lagged")
  expect_error(pssbounds(obs = 6, fstat = 5, case = 7), "Provide either k number of lagged")

  model.bad <- lm(rnorm(500) ~ rnorm(500))
  expect_error(pssbounds(model.bad, "supply either a dynardl model"))
})

test_that("Correct values pulled", {
  # Too many combinations from the matrix to test all. Just pick a lot.
  expect_output(pssbounds(obs = 6, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 1, k = 1, tstat = 5), "-1.62")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 2, k = 4, tstat = 5), "4.22")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 3, k = 3, tstat = 5), "-4.37")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 4, k = 4, tstat = 5), "3.10")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 6, fstat = 5, case = 5, k = 9, tstat = 5), "-5.79")

  expect_output(pssbounds(obs = 31, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 1, k = 3, tstat = 5), "4.84")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 2, k = 4, tstat = 5), "2.95")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 3, k = 7, tstat = 5), "-4.23")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 4, k = 3, tstat = 5), "3.94")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 31, fstat = 5, case = 5, k = 1, tstat = 5), "-3.40")

  expect_output(pssbounds(obs = 39, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 1, k = 4, tstat = 5), "-3.26")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 2, k = 4, tstat = 5), "3.40")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 3, k = 9, tstat = 5), "-4.56")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 4, k = 7, tstat = 5), "3.65")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 39, fstat = 5, case = 5, k = 5, tstat = 5), "4.21")

  expect_output(pssbounds(obs = 44, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 1, k = 7, tstat = 5), "2.54")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 2, k = 3, tstat = 5), "3.08")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 3, k = 6, tstat = 5), "2.76")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 4, k = 3, tstat = 5), "3.82")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 44, fstat = 5, case = 5, k = 6, tstat = 5), "-4.69")

  expect_output(pssbounds(obs = 49, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 1, k = 3, tstat = 5), "3.63")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 2, k = 5, tstat = 5), "3.26")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 3, k = 4, tstat = 5), "-2.86")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 4, k = 4, tstat = 5), "3.38")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 49, fstat = 5, case = 5, k = 6, tstat = 5), "-4.37")

  expect_output(pssbounds(obs = 54, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 1, k = 5, tstat = 5), "1.81")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 2, k = 4, tstat = 5), "2.76")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 3, k = 6, tstat = 5), "3.49")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 4, k = 5, tstat = 5), "3.66")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 54, fstat = 5, case = 5, k = 4, tstat = 5), "-4.04")


  expect_output(pssbounds(obs = 59, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 1, k = 4, tstat = 5), "2.26")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 2, k = 4, tstat = 5), "3.79")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 3, k = 6, tstat = 5), "5.08")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 4, k = 5, tstat = 5), "3.09")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 59, fstat = 5, case = 5, k = 5, tstat = 5), "4.63")

  expect_output(pssbounds(obs = 64, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 1, k = 2, tstat = 5), "2.17")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 2, k = 3, tstat = 5), "2.98")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 3, k = 4, tstat = 5), "4.19")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 4, k = 3, tstat = 5), "5.84")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 64, fstat = 5, case = 5, k = 6, tstat = 5), "4.36")

  expect_output(pssbounds(obs = 69, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 1, k = 5, tstat = 5), "2.82")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 2, k = 5, tstat = 5), "4.72")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 3, k = 3, tstat = 5), "-4.37")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 4, k = 4, tstat = 5), "4.27")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 69, fstat = 5, case = 5, k = 7, tstat = 5), "3.77")

  expect_output(pssbounds(obs = 74, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 1, k = 5, tstat = 5), "3.34")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 2, k = 1, tstat = 5), "3.78")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 3, k = 5, tstat = 5), "5.21")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 4, k = 4, tstat = 5), "5.38")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 74, fstat = 5, case = 5, k = 2, tstat = 5), "-4.53")

  expect_output(pssbounds(obs = 79, fstat = 5, case = 1, k = 1, tstat = 5), "Small-sample critical values not provided for Case I")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 1, k = 6, tstat = 5), "2.66")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 2, k = 4, tstat = 5), "3.70")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 3, k = 1, tstat = 5), "critical values not provided for Case III")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 3, k = 2, tstat = 5), "6.78")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 4, k = 3, tstat = 5), "3.84")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 5, k = 1, tstat = 5), "critical values not provided for Case V")
  expect_output(pssbounds(obs = 79, fstat = 5, case = 5, k = 3, tstat = 5), "-4.73")

  expect_output(pssbounds(obs = 100, fstat = 5, case = 1, k = 1, tstat = 5), "Asymptotic critical values used")
  expect_output(pssbounds(obs = 222, fstat = 5, case = 1, k = 1, tstat = 5), "3.15")
  expect_output(pssbounds(obs = 3232, fstat = 5, case = 2, k = 1, tstat = 5), "Critical values do not currently exist for Case II")
  expect_output(pssbounds(obs = 2335, fstat = 5, case = 2, k = 5, tstat = 5), "2.39")
  expect_output(pssbounds(obs = 235352, fstat = 5, case = 3, k = 1, tstat = 5), "Asymptotic critical values used")
  expect_output(pssbounds(obs = 1224, fstat = 5, case = 3, k = 3, tstat = 5), "-3.43")
  expect_output(pssbounds(obs = 42323, fstat = 5, case = 4, k = 1, tstat = 5), "Critical values do not currently exist for Case IV")
  expect_output(pssbounds(obs = 233, fstat = 5, case = 4, k = 6, tstat = 5), "4.39")
  expect_output(pssbounds(obs = 5434, fstat = 5, case = 5, k = 1, tstat = 5), "Asymptotic critical values used")
  expect_output(pssbounds(obs = 453, fstat = 5, case = 5, k = 7, tstat = 5), "2.38")

})
