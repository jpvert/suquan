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
#' @param f_init The initial distribution for the quantile transformation 
#'   (default is a Gaussian cdf)
#' @param maxiter Maximum number of times the optimization loop (once in the 
#'   quantile distribution and once in the model) is performed (default 
#'   \code{10})
#' @param eps Stopping criterion: the loop will stop when the norm of the 
#'   difference between the quantile distribution after two successive 
#'   optimization is smaller than \code{eps} (default \code{1e-6})
#' @param use.glmnet Whether the \code{glmnet} package should be used to fit the
#'   model for a given quantile distribution. If \code{FALSE}, then the less
#'   optimized \code{apg} optimizer is used (default \code{TRUE}).
#'   
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

suquan <- function(x, y, family=c("gaussian", "binomial"), penalty="elasticnet", lambda=1, intercept=TRUE, f_init = NULL, maxiter = 10, eps = 1e-6, use.glmnet=TRUE, opts=list()) {
    
    p <- ncol(x)
    n <- nrow(x)
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
    rankx <- t(apply(x,1, function(v){rank(v, na.last=NA, ties.method='first')}))
    orderx <- t(apply(rankx, 1, function(v) {match(seq(p),v)}))
    
    # Main loop
    for (iter in seq(maxiter)) {
        
        cat('iter ',iter,'\n')
        f_old <- f

        # Optimize model (b,a0) for a fixed quantile function f
        newx <- matrix(f[rankx],nrow=n)
        if (use.glmnet) {
            mm <- glmnet(newx, y, family=family, intercept=intercept, standardize=FALSE, lambda=lambda, alpha=opts[["alpha"]])
            m <- list(b=coef(mm)[-1,], a0=coef(mm)[1,])
            
        } else {
            if (iter >= 2){
                opts$X_INIT <- m[['b']] 
            }
            m <- glm.apg(newx, y, family=family, penalty=penalty, intercept=intercept, lambda=lambda, opts=opts)
        }

        # Optimize quantile function f for fixed model b
        if (iter < maxiter) {
            newx <- matrix(m[["b"]][orderx],nrow=n)
            if (iter >= 2){
                opts$X_INIT <- c(m2[['b']], m2[['a0']])
            }
            m2 <- glm.apg(newx, y, family=family, penalty="boundednondecreasing", intercept=intercept, opts=opts)
            f <- m2[["b"]]
        }
        
        if (sqrt(sum((f-f_old)^2)) < eps) break
        if (sum(f == mean(f)) == length(f)) break
    }
    
    return(list(f=f_old, b=m[["b"]], a0=m[["a0"]]))
    
}
