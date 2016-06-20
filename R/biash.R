


#' Use control genes to estimate the unimodal inflation. Follow this
#' by running ash on the non-control genes.
#'
#' @inheritParams stramash
#' @param ctl A vector of logicals of length \code{length(betahat)}. A
#'     \code{TRUE} at position \eqn{i} indicates that \eqn{\beta_i =
#'     0}, i.e. that it is a control gene. A \code{FALSE} indicates
#'     that it is not a control gene.
#' @param stramash_args A list containing additional arguments to pass to
#'     \code{\link{stramash.workhorse}}
#'
biash <- function(betahat, ctl, sebetahat, df = NULL, stramash_args = list()) {
    p <- length(betahat)
    assertthat::are_equal(length(ctl), p)
    assertthat::are_equal(length(sebetahat), p)
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(is.list(stramash_args))

    if (!is.null(stramash_args$method)) {
        message("Overwriting your \"method\" option in stramash_args")
    }

    if (!is.null(df)) {
        if (length(df) == 1) {
            df <- rep(df, p)
        }
        df_ctl     <- df[ctl]
        df_unknown <- df[!ctl]
    } else {
        df_ctl     <- NULL
        df_unknown <- NULL
    }

    stramash_args$betahat     <- betahat[ctl]
    stramash_args$sebetahat   <- sebetahat[ctl]
    stramash_args$df          <- df_ctl
    stramash_args$method      <- "shrink"
    stramash_args$mixcompdist <- "normal"
    stramash_out <- do.call(what = stramash.workhorse, args = stramash_args)

    inflation_dist <- stramash_out$fitted.g

    betahat_unknown    <- betahat[!ctl]
    sebetahat_unknown <- sebetahat[!ctl]

    which_pi0 <- inflation_dist$pi < 10 ^ -12
    inflation_dist$mean <- inflation_dist$mean[!which_pi0]
    inflation_dist$sd <- inflation_dist$sd[!which_pi0]
    inflation_dist$pi <- inflation_dist$pi[!which_pi0]

    errordist <- list()
    for(index in 1:length(betahat_unknown)) {
        current_dist <- inflation_dist
        current_dist$sd <- sqrt(inflation_dist$sd ^ 2 + sebetahat_unknown[index] ^ 2)
        errordist[[index]] <- current_dist
    }


    stramash_new_args           <- list()
    stramash_new_args$betahat   <- betahat_unknown
    stramash_new_args$errordist <- errordist
    stramash_new_args$method    <- "fdr"
    stramash_new_args$df        <- df_unknown

    stramash_new_out <- do.call(what = stramash.workhorse, args = stramash_new_args)
    stramash_new_out$inflation_dist <- inflation_dist

    return(stramash_new_out)
}
