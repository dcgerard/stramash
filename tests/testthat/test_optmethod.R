library(stramash)
context("optmethod")

test_that("optmethod switches to EM when given negative probabilities", {
	set.seed(777)
        betahat <- c(0.7775, 0.6091, 0.6880, 0.6897, 0.6565, 0.7505,
                     0.7125, 0.7201, 0.7498, 0.7553)
        sebetahat <- rep(0.01, length = length(betahat))
        stramash_out <- stramash.workhorse(betahat = betahat, sebetahat = sebetahat,
                             mixcompdist = "uniform", errordist = NULL)
        expect_true(min(stramash_out$fitted.g$pi) > -10 ^ -12)
})
