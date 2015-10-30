#' Predict with a supervised quantiled linear model.
#' 
#' Make prediction on new samples with a quantiled linear model.
#' 
#' @param m The model, a list that contains three fields: \code{f}, the quantile
#'   function, \code{b}, the linear model, \code{b0}, the intercept
#' @param x Samples (one per row) for which the prediction is made.
#'   
#' @return A vector of predictions, for all samples.
#'   
#' @export
#' 
#' @examples
#' n <- 100
#' p <- 40
#' x <- matrix(rnorm(n*p),n,p)
#' y <- rbinom(n,1,0.5)*2-1
#' lambda <- 0.2*max(abs(crossprod(y,x)))/n
#' m <- suquan(x, y, family="gaussian", penalty="elasticnet", opts=list(alpha=0))
#' predict.suquan(m,x)
#' 
predict.suquan <- function(m, x) {
    
    rankx <- t(apply(x,1,rank))
#    matrix(m$b[rankx],nrow=nrow(x))%*%m$f + m$a0
    matrix(m$f[rankx],nrow=nrow(x)) %*% m$b + m$a0    
}
    