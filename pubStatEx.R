logLik.lm <- function(object, REML = FALSE, ...)
{
  # log likelihood function for mlm objects, i.e. lm with multiple responses
  # modified from logLik.lm to handle multiple response lm (mlm) object 
  # reference: https://github.com/lgautier/R-3-0-branch-alt/blob/master/src/library/stats/R/logLik.R
  # version 1.0, 20150103, YF Li

  res <- object$residuals # not resid(object) because of NA methods
  p <- object$rank
  if(inherits(object, "mlm")){
    cat("'logLik.lm' has been upgrades to support multiple responses")
    N.s = ncol(res)
    N <- nrow(res)
  }else{
    N.s = 1
    N <- length(res)    
  }
  
  if(is.null(w <- object$weights)) {
    w <- matrix(1, nrow = N, ncol=N.s)
  } else {
    excl <- w == 0
    if (any(excl)) {
      N <- colSums(!excl)
      N.s = rowSums(!excl)
    }
  }
  N0 <- N
  if(REML) N <- N - p
  val <- .5* (colSums(log(w)) - N * (log(2 * pi) + 1 - log(N) +
                                       log(colSums(w*res^2))))
  if(REML) val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
  attr(val, "nall") <- N0 # NB, still omits zero weights
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}
