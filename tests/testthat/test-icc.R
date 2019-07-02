# Tests for lme4icc package
# 2019-01-22

require(lme4)
require(nlme)

context("Tests for icc() with object of class merMod")

# Example from lme4::VarCorr() helpfile - removing random slope
data(Orthodont, package = "nlme")
fm1 <- lme4::lmer(distance ~ age + (1 | Subject), data = Orthodont)
vc <- VarCorr(fm1) ## default print method: standard dev and corr
# print(vc, comp = c("Variance", "Std.Dev."), digits = 2)
vcdf <- as.data.frame(vc, order = "lower.tri")

test_that("ICC for lme4::VarCorr() helpfile example",
          {
            iccfit <- vcdf[1,4] / (vcdf[1,4] + vcdf[2,4])
            iccpkg <- iccmlm::icc(fm1)
            expect_equivalent(iccpkg, iccfit, tol =  1e-3)
          })
