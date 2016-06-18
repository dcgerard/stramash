context("+Uniform")

test_that("mixcompdist=+uniform gives all non-negative values for b and zero for a", {
    set.seed(1); z=rnorm(10); z.stramash=stramash.workhorse(betahat = z,sebetahat = rep(1, length(z)),mixcompdist="+uniform")
  k = length(z.stramash$fitted.g$pi)
  expect_true(all(z.stramash$fitted.g$b >= rep(0,k)))
  expect_equal(z.stramash$fitted.g$a,rep(0,k))
})

test_that("mixcompdist=-uniform gives all non-positive values for a and zero for b", {
  set.seed(1); z=rnorm(10); z.stramash=stramash.workhorse(betahat = z,sebetahat = rep(1, length(z)),mixcompdist="-uniform")
  k = length(z.stramash$fitted.g$pi)
  expect_equal(z.stramash$fitted.g$b,rep(0,k))
  expect_true(all(z.stramash$fitted.g$a <= 0))
})
