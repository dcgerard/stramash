context("VB")

test_that("stramash emits error if try to use prior=unit without VB", {
    expect_error(stramash.workhorse(betahat = rnorm(10), sebetahat = 1,
                                    optmethod="mixEM", prior = "unit"))
})

test_that("stramash emits error if try to use prior<1 without VB", {
    expect_error(stramash.workhorse(betahat = rnorm(10), sebetahat = 1,
                                    optmethod="mixEM", nullweight = 0.5))
})
