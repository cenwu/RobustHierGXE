#=========================================== gALAL =====================================#

gALAL <- function (xx, y, tau, beta=0.9995, eps, wpp) 
{
    n <- length(y)
    p <- ncol(xx)
    if (n != nrow(xx)) 
        stop("xx and y don't match n")
    lambda =   wpp * n
    if (length(lambda) != p) 
        stop(paste("lambda must be either of length ", p, " or length one"))
    if (any(lambda < 0)) 
        stop("negative lambdas disallowed")
    if (tau < eps || tau > 1 - eps) 
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    R <- matrix(0,length(lambda),length(lambda))
    diag(R) <- lambda
    index <- which(lambda != 0)
    len = length(index)
  
    if (len>0)
       {
    R <- R[index, ]
    if (len==1){R=as.matrix(t(R))}
    r <- rep(0, len)
    X <- rbind(xx, R)
    Y <- c(y, r)
    N <- length(Y)
    rhs <- (1 - tau) * apply(xx, 2, sum) + 0.5 * apply(R, 2, sum)
    d <- rep(1, N)
    u <- rep(1, N)
    wn <- rep(0, 10 * N)
    wn[1:N] <- c(rep(1 - tau, n), rep(0.5, nrow(R)))
    z <- .Fortran("rqfnb", as.integer(N), as.integer(p), a = as.double(t(as.matrix(X))), 
        c = as.double(-Y), rhs = as.double(rhs), d = as.double(d), 
        as.double(u), beta = as.double(beta), eps = as.double(eps), 
        wn = as.double(wn), wp = double((p + 3) * p), aa = double(p * p), 
        it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
    if (z$info != 0) 
        stop(paste("Error info = ", z$info, "in stepy2: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(xx)[[2]]
    residuals <- y - xx %*% coefficients
      }
    list(coefficients = coefficients, residuals = residuals)
}