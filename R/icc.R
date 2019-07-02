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
#' 
#' @export
print.iccmlm <- function(x, digits = getOption("digits"), percent = FALSE, ...) {
  if (!percent) {
    cat("\n Intra-class correlation coefficient:", print(x, digits = digits))  
  } else {
    cat("\n Intra-class correlation coefficient:", print(x * 100, digits = digits), "%")
  }
}
