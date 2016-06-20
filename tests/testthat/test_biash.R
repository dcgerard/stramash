library(stramash)
context("biash")

test_that("biash works", {
    set.seed(223)
    p <- 100
    beta <- rnorm(p)
    ctl <- rep(FALSE, length = p)
    ctl[1:50] <- TRUE
    beta[ctl] <- 0
    sebetahat <- sqrt(rchisq(p, df = 5) / 5)
    betahat <- beta + rnorm(n = p, sd = sebetahat)

    bout <- biash(betahat = betahat, ctl = ctl, sebetahat = sebetahat - min(sebetahat) + 0.05)
    mixsd(bout$inflation_dist)
    min(sebetahat)
}
)
