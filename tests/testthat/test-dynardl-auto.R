context("test-dynardl-auto")

test_that("Stop for non-dynardl objects", {
  y <- rnorm(500)
  expect_error(dynardl.auto.correlated(y), "provide a dynardl model object")
})

test_that("Correct fit statistics reported", {
  dyn.out.1.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)  
  expect_output(dynardl.auto.correlated(dyn.out.1.1), paste(round(AIC(dyn.out.1.1$model), 3)))
  expect_output(dynardl.auto.correlated(dyn.out.1.1), paste(round(BIC(dyn.out.1.1$model), 3)))
  # Test rounding because why not
  expect_output(dynardl.auto.correlated(dyn.out.1.1, digits = 5), paste(round(AIC(dyn.out.1.1$model), 5)))
  expect_output(dynardl.auto.correlated(dyn.out.1.1, digits = 5), paste(round(BIC(dyn.out.1.1$model), 5)))
})

test_that("Correct bgtest reported", {  
  dyn.out.1.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)  
  expect_output(dynardl.auto.correlated(dyn.out.1.1, digits = 3), paste(round(bgtest(dyn.out.1.1$model)$statistic, 3)))
  expect_output(dynardl.auto.correlated(dyn.out.1.1, digits = 3, order = 5), paste(round(bgtest(dyn.out.1.1$model, order = 5)$statistic, 3)))
  expect_output(dynardl.auto.correlated(dyn.out.1.1, digits = 3, bg.type = "F"), paste(round(bgtest(dyn.out.1.1$model, type = "F")$statistic, 3)))
})

test_that("Correct shapiro test reported", { 
  dyn.out.1.1 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)   
  expect_output(dynardl.auto.correlated(dyn.out.1.1, digits = 3), paste(round(shapiro.test(dyn.out.1.1$model$residuals)$statistic, 3)))
})

test_that("Information out correctly", { 
  dyn.out.1.2 <- dynardl(concern ~ incshare10 + urate, data = ineq, 
    lags = list("concern" = 1, "incshare10" = 1, "urate" = 1), ec = TRUE,
    simulate = FALSE)  
  info <- dynardl.auto.correlated(dyn.out.1.2, object.out = TRUE)
  expect_equal(names(info), c("bg", "sw", "AIC", "BIC", "logLik"))
})

test_that("Does clear autocorrelation fail the test?", { 
  dyn.out.1.3 <- dynardl(urate ~ concern, data = ineq, 
    lags = list("urate" = 1), ec = FALSE,
    simulate = FALSE)   
  expect_output(dynardl.auto.correlated(dyn.out.1.3), "Breusch-Godfrey test indicates we reject the null hypothesis")
  expect_output(dynardl.auto.correlated(dyn.out.1.3), "Shapiro-Wilk test indicates we reject the null hypothesis")
})

