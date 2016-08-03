#' Take n draws from G ~ DP(alpha, H)
#'
#' @examples
#'     X1 <- rdp(1000, mean=4, sd=20)
#'     X2 <- rdp(1000, H=rbeta, shape1=2, shape2=2)
#' @param n number of draws from Dirichlet Process
#' @param alpha precision parameter
#' @param H baseline distribution
#' @param k number of stick breaks and draws from baseline distribution
#' @param ... named parameters passed to baseline distribution
#' @return A vector of draws from DP
rdp <- function(n, alpha=1, H=rnorm, k=1000, ...) {

    # beta_k ~ Beta(1, alpha)
    beta.k <- rbeta(k, 1, alpha)
    
    # pi_k = beta_k prod_1^{k-1} (1-beta_t)
    pi.k <- numeric(k)
    pi.k[1] <- beta.k[1]
    for (j in 2:k) {
        pi.k[j] <- beta.k[j] * prod(1-beta.k[1:(j-1)])
    }

    # theta_k^* ~ H
    theta.star.k <- H(k, ...)

    # G = sum_1^{infty} pi_k delta_{theta_k^*}
    G <- sample(theta.star.k, prob=pi.k, replace=TRUE)

    return(G)
}

#' Take n draws from G ~ DP(alpha, H) using CRP
#'
#' @examples
#'     X1 <- rcrp(k=1000, mean=4, sd=20)
#'     X2 <- rcrp(k=1000, H=rbeta, shape1=2, shape2=2)
#' @param alpha precision parameter
#' @param H baseline distribution
#' @param k number of customers to seat
#' @param ... named parameters passed to baseline distribution
#' @return A vector of draws from DP
rcrp <- function(alpha=1, H=rnorm, k=1000, ...) {

    # number of customers at each table
    m <- numeric(k)
    
    # assign first customer to first table w/ prob 1
    m[1] <- 1
    
    # now conditionally assign (k-1) customers to tables
    for (i in 2:k) {
        # choose already seated table w/ prob prop to customers
        prob <- m
        
        # choose new table with prob proportional to alpha
        prob[sum(m > 0) + 1] <- alpha
        
        m <- m + rmultinom(1, 1, prob)
    }
    
    m <- m[m > 0]
    m <- m / sum(m)
    
    # phi_k ~ H
    phi.k <- H(length(m), ...)
    
    return(cbind(m, phi.k))
}
