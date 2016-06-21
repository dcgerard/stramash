library(stramash)
context("Mixture Approximations")

test_that("Laplace distribution is well approximated", {
    set.seed(747)
    mu <- 1
    sig <- 5
    n <- 10000

    lmix <- laplace_to_mix(mu = mu, sig = sig, gridsize = 10000)

    expect_equal(pnormalmix(q = 5, lmix),
                 VGAM::plaplace(q = 5, location = mu, scale = sig),
                 tol = 10 ^ -5)

    expect_equal(pnormalmix(q = -1, lmix),
                 VGAM::plaplace(q = -1, location = mu, scale = sig),
                 tol = 10 ^ -4)

    expect_equal(pnormalmix(q = 0, lmix),
                 VGAM::plaplace(q = 0, location = mu, scale = sig),
                 tol = 10 ^ -4)

    pnorm(q = 0, mean = mu, sd = sig)
    VGAM::plaplace(q = 0, location = mu, scale = sig)
    pnormalmix(q = 0, lmix)

    lvar <- sum((lmix$mean ^ 2 + lmix$sd ^ 2) * lmix$pi) - sum(lmix$mean * lmix$pi) ^ 2
    theo_var <- 2 * sig ^ 2
    lvar
    theo_var

    ## lsamp <- rnormalmix(n = 1000, mixdense = lmix)
    ## vsamp <- VGAM::rlaplace(n = 1000, location = mu, scale = sig)
    ## plot(sort(lsamp), sort(vsamp))
    ## abline(0, 1)
}
)


test_that("t approximation is accurate", {
    set.seed(600)
    mu <- 0
    sig <- 1
    df <- 6
    n <- 10000

    tmix <- t_to_mix(mu = mu, sig = sig, df = df, gridsize = 10000)

    expect_equal(pnormalmix(q = 3, tmix),
                 stats::pt(q = 3, df = df),
                 tol = 10 ^ -5)

    expect_equal(pnormalmix(q = -1, tmix),
                 stats::pt(q = -1, df = df),
                 tol = 10 ^ -6)

    tvar <- sum((tmix$mean ^ 2 + tmix$sd ^ 2) * tmix$pi) - sum(tmix$mean * tmix$pi) ^ 2
    theo_var <- sig * df / (df - 2)
    expect_equal(theo_var, tvar, tol =  10 ^ -3)
}
)


test_that("general quantile to mix works", {
    set.seed(101)
    mixdense <- quantile_to_mix(q_func = qnorm, gridsize = 1000)
    x <- seq(-8, 8, length = 100)
    y <- dunimix(x = x, mixdense = mixdense)
    z <- dnorm(x)
    expect_equal(y, z, tol = 10 ^ -2)

    plot(x, y, type = "l")
    lines(x, z, col = 2, lty = 2)

    set.seed(51)
    p <- 100
    beta <- rnorm(p)
    betahat <- beta + rnorm(p)

    mixdense <- quantile_to_mix(q_func = qnorm, gridsize = 20)
    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- mixdense
    }
    sout <- stramash.workhorse(betahat = betahat, errordist = errordist, mixcompdist = "normal")
    ash_out <- ashr::ash.workhorse(betahat = betahat, sebetahat = 1, mixcompdist = "normal")

    plot(sout$lfdr, ash_out$lfdr)
    abline(0, 1)

    plot(ash_out$PosteriorMean, sout$PosteriorMean)
    abline(0, 1)
    plot(beta, sout$PosteriorMean)
    plot(beta, ash_out$PosteriorMean)
}
)
