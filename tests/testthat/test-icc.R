# Tests for lme4icc package

require(lme4)
require(nlme)

context("Tests for icc() with object of class merMod")

# Example from lme4::VarCorr() helpfile - removing random slope
data(Orthodont, package = "nlme")
fm1 <- lme4::lmer(distance ~ age + (1 | Subject), data = Orthodont)
vc <- lme4::VarCorr(fm1) ## default print method: standard dev and corr
# print(vc, comp = c("Variance", "Std.Dev."), digits = 2)
vcdf <- as.data.frame(vc, order = "lower.tri")
iccfit <- vcdf[1, 4] / (vcdf[1, 4] + vcdf[2, 4])

test_that("Check no error from running icc function", {
  expect_error(icc(fm1), NA)
})

iccpkg <- icc(fm1)

test_that("ICC for lme4::VarCorr() helpfile example", {
  expect_equivalent(iccpkg$icc, iccfit, tol =  1e-3)
})

test_that("Pass model of incorrect class to icc()", {
  # Example from lm helpfile
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  lm.D9 <- lm(weight ~ group)
  expect_error(icc(lm.D9))
})

test_that("Test print method for class iccmlm", {
  expect_output(print(iccpkg))
})

test_that("Check exactly what is printed", {
  expect_output(print(iccpkg, digits = 2), 
                "ICC: 0.69")  
})

test_that("Check percent option to print gives percent sign at end", {
  expect_output(print(iccpkg, percent = TRUE, digits = 2),
                "ICC: 68.57%")  
})

test_that("Check bootstrap default CI", {
  iccbci <- icc(fm1, bci = TRUE, seed = 20200510)
  expect_equal(iccbci$bci[2], .5496791, tol = 1e-4)
  expect_equal(iccbci$bci[3], .8217991, tol = 1e-4)
})

test_that("Check bootstrap 90% CI", {
  iccbci <- icc(fm1, bci = TRUE, seed = 20200510, conf = .9)
  expect_equal(iccbci$bci[2], .5715540, tol = 1e-4)
  expect_equal(iccbci$bci[3], .7999, tol = 1e-4)
})
