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

# compare alpha parameters
set.seed(90210)
n <- 10000
G.0 <- rbeta(n, 2, 2)
Y1 <- rdp(n, alpha=0.01, H=rbeta, shape1=2, shape2=2)
Y2 <- rdp(n, alpha=0.05, H=rbeta, shape1=2, shape2=2)
Y3 <- rdp(n, alpha=0.5, H=rbeta, shape1=2, shape2=2)
Y4 <- rdp(n, alpha=1, H=rbeta, shape1=2, shape2=2)
Y5 <- rdp(n, alpha=10, H=rbeta, shape1=2, shape2=2)
Y6 <- rdp(n, alpha=100, H=rbeta, shape1=2, shape2=2)
Y7 <- rdp(n, alpha=1000, H=rbeta, shape1=2, shape2=2)
Y8 <- rdp(n, alpha=10000, H=rbeta, shape1=2, shape2=2)

library(ggplot2)
theme_set(theme_classic())

data <- rbind(data.frame(Y=G.0, what="baseline"),
              data.frame(Y=Y1, what="alpha=0.01"),
              data.frame(Y=Y2, what="alpha=0.05"),
              data.frame(Y=Y3, what="alpha=0.5"),
              data.frame(Y=Y4, what="alpha=1"),
              data.frame(Y=Y5, what="alpha=10"),
              data.frame(Y=Y6, what="alpha=100"),
              data.frame(Y=Y7, what="alpha=1000"),
              data.frame(Y=Y8, what="alpha=10000"))

ggplot(data, aes(Y)) +
    geom_density() +
    facet_wrap(~what)
