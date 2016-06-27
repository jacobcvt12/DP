# make n draws from G ~ DP(alpha, H)
# where alpha is [0.05, 0.5, 1] and H ~ Normal(0, 1)
n <- 500

dp <- function(n, alpha=1, H="normal", k=100, ...) {
    G <- numeric(n)

    .dens <- list("normal"=list("r"=rnorm,
                                "d"=dnorm),
                  "t"=list("r"=rt,
                           "d"=dt))

    for (i in 1:n) {
        # beta_k ~ Beta(1, alpha)
        beta.k <- rbeta(k, 1, alpha)
        
        # pi_k = beta_k prod_1^{k-1} (1-beta_t)
        pi.k <- numeric(k)
        pi.k[1] <- beta.k[1]
        for (j in 2:k) {
            pi.k[j] <- beta.k[j] * prod(1-beta.k[1:(j-1)])
        }

        theta.star.k <- H(k, ...)

        G[i] <- sum(pi.k * theta.star.k)
    }

    return(G)
}
