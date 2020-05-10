#' Calculate ICC from an object of class lmerMod or glmerMod.
#'
#' @param model Object of class either lmerMod or glmerMod (fitted with family = "gaussian").
#' 
#' @export
#' @return Object of class iccmlm containing the ICC.
icc <- function(model) {

  clmod <- class(model)[1]

  if (!(clmod %in% c("lmerMod", "glmerMod"))) {
    stop("Currently, this function only works with models of class lmerMod or glmerMod (with family gaussian) fit by the lme4 package.")
  }
  
  if (stats::family(model)$family != "gaussian") {
    stop("The family must be gaussian.")
  }

  revars <- names(model@flist)
  reintvar <- revars[1]
    
  vc <- lme4::VarCorr(model)
  resvar <- attr(vc, "sc")^2
  
  intvarexpr <- paste0("vc$", reintvar, "[1,1]")
  
  intvar <- eval(parse(text = intvarexpr))
  icc <- intvar/(resvar + intvar)  
  
  class(icc) <- append(class(icc), "iccmlm")
  
  return(icc)
}

#' print method for an object of class iccmlm
#'
#' @param x Object of class iccmlm
#' @param digits minimal number of \emph{significant} digits, see \code{\link{print.default}}.
#' @param percent Logical indicating whether the ICC should be reported as a percentage.
#' 
#' @export
print.iccmlm <- function(x, digits = getOption("digits"), percent = FALSE, ...) {
  if (!percent) {
    cat("\nIntra-class correlation coefficient: ")
    cat(round(x[1], digits = digits))
  } else {
    cat("\nIntra-class correlation coefficient: ")
    cat(round(x[1] * 100, digits = digits))
    cat("%")
  }
  invisible(x)
}

#' Calculate ICC with indices (i.e. subset of the data)
iccmodel <- function(data, formula, indices){
  d <- data[indices, ]
  fit <- lme4::lmer(mathach ~ (1 | school), data = d)
  vc <- lme4::VarCorr(fit)
  resvar <- attr(vc, "sc")^2
  intvar <- vc$school[1,1]
  icc <- intvar/(resvar + intvar)
  return(icc)
}

#' Bootstrap standard error for the ICC
#' 
#' @export
bootse.iccmlm <- function(...) {
  if (!exists("R")) R <- 50
  bsrun <- boot::boot(data = dat, iccmodel, R = R, ...)
  return(bsrun)
}
