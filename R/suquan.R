#' Fit a supervised quantiled linear model.
#' 
#' Fit a generalized linear model via penalized maximum likelihood, with joint
#' optimization of full quantile normalization.
#' 
#' @param x The input matrix, each row is a sample, each column a feature.
#' @param y The response variable. Quantitative for \code{family="gaussian"},
#'   binary with values \code{+1} and \code{-1} for \code{family="binomial"}
#' @param family The response type. For \code{family="gaussian"}, ...
#' @param penalty The penalty type.
#' @param lambda The scaling of the penalty (default \code{1})
#' @param intercept Should intercept(s) be fitted (default=\code{TRUE}) or set
#'   to zero (\code{FALSE})
#' @param opts List of parameters, which must include: \itemize{ \item }
#'   
#' @return toto
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
#' 
suquan <- function(x, y, family=c("gaussian", "binomial"), penalty="elasticnet", lambda=1, intercept=TRUE, f_init = NULL, opts=list()) {
    
    MAXITER <- 100
    EPS <- 1e-6
    p <- ncol(x)
    penalty <- match.arg(penalty)

    # Initiate the quantile function
    if (is.null(f_init)) {
        # By default, itiniate with Gaussian distribution
        f <- qnorm(seq(from=1/(p+1), by=1/(p+1), length.out=p))
    } else {
        if (length(f_init) != p) {
            stop('The length of the initial quantile function should be the same as the dimension of samples in x.')
        }
        f <- f_init
    }
    
    # Sort samples
    rankx <- t(apply(x,1,rank))
    orderx <- t(apply(rankx, 1, function(v) {match(seq(p),v)}))
    
    # Main loop
    for (iter in seq(MAXITER)) {
        
        cat('iter ',iter,'\n')
        f_old <- f

        # Optimize model (b,a0) for a fixed quantile function f
        newx <- matrix(f[orderx],nrow=n)
        m <- glm.apg(newx, y, family=family, penalty=penalty, intercept=intercept, lambda=lambda, opts=opts)

        # Optimize quantile function f for fixed model b
        newx <- matrix(m[["b"]][rankx],nrow=n)
        m2 <- glm.apg(newx, y, family=family, penalty="boundednondecreasing", intercept=intercept, opts=opts)
        
        f <- m2[["b"]]
        
        if (sqrt(sum((f-f_old)^2)) < EPS) break
    }
    
    return(list(f=f_old, b=m[["b"]], a0=m[["a0"]]))
    
}