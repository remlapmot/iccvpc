#' Calculate ICC from an object of class lmerMod or glmerMod
#'
#' @param model Object of class either `lmerMod` or `glmerMod` (fitted with `family = "gaussian"`).
#' 
#' @export
#' @return Object of class `iccmlm` containing the ICC.
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
#' @param ... Additional arguments passed to [`print`].
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

#' Calculate ICC with indices (i.e. on a subset of the data)
#' 
#' @param data The data the model is fitted to.
#' @param indices The indices that boot::boot() uses to select the replicate samples.
#' @param fit The fitted model object which the ICC is for.
#' @importFrom stats update
#' @keywords internal
iccmodel <- function(data, indices, fit){
  d <- data[indices, ]
  fit <- update(fit, data = d)
  icc(fit)
}

#' Bootstrap standard error for the ICC
#' 
#' @param x An object of class `iccmlm`.
#' @param ... Further arguments passed to [boot::boot()].
#' @export
bootci <- function(x, ...) {
  UseMethod("bootci")  
}
 
#' @export
bootci.iccmlm <- function(x, fit, data, seed, ...) {
  if (!missing(seed)) set.seed(seed)
  if (!exists("R")) R <- 50
  bsrun <- boot::boot(data = data, 
                      statistic = iccmodel, 
                      R = R, 
                      fit = fit, 
                      ...)
  bias <- mean(bsrun$t) - bsrun$t0
  bci <- boot::boot.ci(bsrun, type = 'norm')
  return(list(bci = c(bci$normal[1], bci$normal[2:3] + bias)))
}
