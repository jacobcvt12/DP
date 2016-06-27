# make n draws from G ~ DP(alpha, H)
# where alpha is [0.05, 0.5, 1] and H ~ Normal(0, 1)
n <- 500

rdp <- function(n, alpha=1, H="normal", k=1000, ...) {

    .dens <- list("normal"=list("r"=rnorm,
                                "d"=dnorm),
                  "t"=list("r"=rt,
                           "d"=dt))

    # beta_k ~ Beta(1, alpha)
    beta.k <- rbeta(k, 1, alpha)
    
    # pi_k = beta_k prod_1^{k-1} (1-beta_t)
    pi.k <- numeric(k)
    pi.k[1] <- beta.k[1]
    for (j in 2:k) {
        pi.k[j] <- beta.k[j] * prod(1-beta.k[1:(j-1)])
    }

    # theta_k^* ~ H
    theta.star.k <- .dens[[H]][["r"]](k, ...)

    # G = sum_1^{infty} pi_k delta_{theta_k^*}
    G <- samply(theta.star.k, prob=pi.k, replace=TRUE)

    return(G)
}
