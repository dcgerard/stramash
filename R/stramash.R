# avoid "no visible binding for global variable" note in CRAN check
# These variables are actually defined in process_args
if (getRversion() >= "2.15.1") utils::globalVariables(c("completeobs",
                                                       "controlinput",
                                                       "sebetahat.orig",
                                                       "excludeindex"))

#' @import ashr
#
#


#' @title Main Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their
#'     standard errors (sebetahat), and applies shrinkage to them,
#'     using Empirical Bayes methods, to compute shrunk estimates for
#'     beta.
#'
#' @details This function is actually just a simple wrapper that
#'     passes its parameters to \code{\link{stramash.workhorse}} which
#'     provides more documented options for advanced use. See readme
#'     for more details.
#'
#' @param betahat a p vector of estimates
#' @param sebetahat a p vector of corresponding standard errors
#' @param errordist A list of objects of either class \code{normalmix}
#'     or \code{unimix}. The length of this list must be the length of
#'     \code{betahat}. \code{errordist[[i]]} is the \eqn{i}th error
#'     distribution of \code{betahat[i]}. Defaults to \code{NULL}, in
#'     which case \code{stramash.workhorse} will assume either a
#'     normal or t likelihood, depending on the value for \code{df}.
#' @param mixcompdist distribution of components in mixture
#'     ("uniform","halfuniform" or "normal"; "+uniform" or
#'     "-uniform"), the default is "uniform". If you believe your
#'     effects may be asymmetric, use "halfuniform". If you want to
#'     allow only positive/negative effects use "+uniform"/"-uniform".
#'     The use of "normal" is permitted only if df=NULL.
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL which is actually treated as
#'     infinity (Gaussian)
#' @param ... Further arguments to be passed to
#'     \code{\link{stramash.workhorse}}.
#'
#' @return stramash returns an object of \code{\link[base]{class}}
#'     "stramash", a list with some or all of the following elements
#'     (determined by outputlevel) \cr
#' \item{fitted.g}{fitted mixture, either a normalmix or unimix}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{logLR}{log[P(D|mle(pi))/P(D|beta==0)]}
#' \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#' \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#' \item{PositiveProb}{A vector of posterior probability that beta is positive}
#' \item{NegativeProb}{A vector of posterior probability that beta is negative}
#' \item{ZeroProb}{A vector of posterior probability that beta is zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{excludeindex}{the vector of index of observations with 0 standard error; if none, then returns NULL}
#' \item{model}{either "EE" or "ET", denoting whether exchangeable effects (EE) or exchangeable T stats (ET) has been used}
#' \item{optmethod}{the optimization method used}
#' \item{data}{a list consisting the input betahat and sebetahat (only included if outputlevel>2)}
#'
#' @seealso \code{\link{stramash.workhorse}} for complete specification of
#'     stramash function
#'
#' @export
#' @examples
#' beta = c(rep(0, 100), stats::rnorm(100))
#' sebetahat = abs(stats::rnorm(200, 0, 1))
#' betahat = stats::rnorm(200, beta, sebetahat)
#' beta.stramash = stramash.workhorse(betahat = betahat, sebetahat = sebetahat)
#' summary(beta.stramash)
#' names(beta.stramash)
#' graphics::plot(betahat, beta.stramash$PosteriorMean, xlim = c(-4, 4), ylim = c(-4, 4))
#'
stramash <- function(betahat, errordist = NULL, sebetahat = NULL,
               mixcompdist = c("uniform", "halfuniform", "normal", "+uniform", "-uniform"),
               df = NULL, ...) {
  return(utils::modifyList(stramash.workhorse(betahat, sebetahat,
                                  mixcompdist = mixcompdist, df = df,
                                  ...), list(call = match.call())))
}


#' @title Detailed Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their
#'     standard errors (sebetahat), and applies shrinkage to them,
#'     using Empirical Bayes methods, to compute shrunk estimates for
#'     beta. This is the more detailed version of stramash for "research"
#'     use.  Most users will be happy with the stramash function, which
#'     provides the same usage, but documents only the main options
#'     for simplicity.
#'
#' @details See readme for more details.
#'
#' @param betahat a p vector of estimates
#' @param sebetahat a p vector of corresponding standard errors
#' @param method specifies how stramash is to be run. Can be "shrinkage"
#'     (if main aim is shrinkage) or "fdr" (if main aim is to assess
#'     fdr or fsr) This is simply a convenient way to specify certain
#'     combinations of parameters: "shrinkage" sets pointmass=FALSE
#'     and prior="uniform"; "fdr" sets pointmass=TRUE and
#'     prior="nullbiased".
#' @param mixcompdist distribution of components in mixture (
#'     "uniform","halfuniform","normal" or "+uniform"), the default
#'     value is "uniform" use "halfuniform" to allow for assymetric g,
#'     and "+uniform"/"-uniform" to constrain g to be
#'     positive/negative.
#' @param optmethod specifies optimization method used. Default is
#'     "mixIP", an interior point method, if REBayes is installed;
#'     otherwise an EM algorithm is used. The interior point method is
#'     faster for large problems (n>2000).
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL(Gaussian)
#' @param nullweight scalar, the weight put on the prior under
#'     "nullbiased" specification, see \code{prior}
#' @param randomstart logical, indicating whether to initialize EM
#'     randomly. If FALSE, then initializes to prior mean (for EM
#'     algorithm) or prior (for VBEM)
#' @param nonzeromode logical, indicating whether to use a non-zero
#'     unimodal mixture(default is "FALSE")
#' @param pointmass logical, indicating whether to use a point mass at
#'     zero as one of components for a mixture distribution
#' @param prior string, or numeric vector indicating Dirichlet prior
#'     on mixture proportions (defaults to "uniform", or (1,1...,1);
#'     also can be "nullbiased" (nullweight,1,...,1) to put more
#'     weight on first component), or "unit" (1/K,...,1/K) [for
#'     optmethod=mixVBEM version only]
#' @param mixsd vector of sds for underlying mixture components
#' @param gridmult the multiplier by which the default grid values for
#'     mixsd differ by one another. (Smaller values produce finer
#'     grids)
#' @param outputlevel determines amount of output [0=just fitted g;
#'     1=also PosteriorMean and PosteriorSD; 2= everything usually
#'     needed; 3=also include results of mixture fitting procedure
#'     (includes matrix of log-likelihoods used to fit mixture); 4=
#'     output additional things required by flash (flash.data)]
#' @param g the prior distribution for beta (usually estimated from
#'     the data; this is used primarily in simulated data to do
#'     computations with the "true" g)
#' @param fixg if TRUE, don't estimate g but use the specified g -
#'     useful for computations under the "true" g in simulations
#' @param VB (deprecated, use optmethod) whether to use Variational
#'     Bayes to estimate mixture proportions (instead of EM to find
#'     MAP estimate), see \code{\link{mixVBEM}} and
#'     \code{\link{mixEM}}
#' @param cxx flag (deprecated, use optmethod) to indicate whether to
#'     use the c++ (Rcpp) version. After application of Squared
#'     extrapolation methods for accelerating fixed-point iterations
#'     (R Package "SQUAREM"), the c++ version is no longer faster than
#'     non-c++ version, thus we do not recommend using this one, and
#'     might be removed at any point.
#' @param model c("EE","ET") specifies whether to assume exchangeable
#'     effects (EE) or exchangeable T stats (ET).
#' @param control A list of control parameters for the optmization
#'     algorithm. Default value is set to be control.default=list(K =
#'     1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4,
#'     kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). User
#'     may supply changes to this list of parameter, say,
#'     control=list(maxiter=10000,trace=TRUE)
#' @param errordist A list of objects of either class \code{normalmix}
#'     or \code{unimix}. The length of this list must be the length of
#'     \code{betahat}. \code{errordist[[i]]} is the \eqn{i}th error
#'     distribution of \code{betahat[i]}. Defaults to \code{NULL}, in
#'     which case \code{stramash.workhorse} will assume either a normal or
#'     t likelihood, depending on the value for \code{df}.
#' @param likelihood One of the pre-specified likelihoods available.
#' @param gridsize The size of the grid if you are using one of the
#'     pre-specified likelihoods.
#'
#' @return stramash returns an object of \code{\link[base]{class}} "stramash", a list
#' with some or all of the following elements (determined by outputlevel) \cr
#' \item{fitted.g}{fitted mixture, either a normalmix or unimix}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{logLR}{log[P(D|mle(pi))/P(D|beta==0)]}
#' \item{PosteriorMean}{A vector consisting the posterior mean of beta
#' from the mixture}
#' \item{PosteriorSD}{A vector consisting the corresponding posterior
#' standard deviation}
#' \item{PositiveProb}{A vector of posterior probability that beta is
#' positive}
#' \item{NegativeProb}{A vector of posterior probability that beta is
#' negative}
#' \item{ZeroProb}{A vector of posterior probability that beta is
#' zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#' \item{svalue}{A vector of s values}
#' \item{call}{a call in which all of the specified arguments are
#' specified by their full names}
#' \item{excludeindex}{the vector of index of observations with 0
#' standard error; if none, then returns NULL}
#' \item{model}{either "EE" or "ET", denoting whether exchangeable
#' effects (EE) or exchangeable T stats (ET) has been used}
#' \item{optmethod}{the optimization method used}
#' \item{data}{a list consisting the input betahat and sebetahat (only
#' included if outputlevel>2)}
#' \item{fit}{a list containing results of mixture optimization, and
#' matrix of component log-likelihoods used in this optimization}
#'
#' @seealso \code{\link{stramash}} for simplified specification of
#'     stramash function
#'
#' @export
#' @examples
#' beta = c(rep(0, 100), stats::rnorm(100))
#' sebetahat = abs(stats::rnorm(200, 0, 1))
#' betahat = stats::rnorm(200, beta, sebetahat)
#' beta.stramash = stramash.workhorse(betahat = betahat, sebetahat = sebetahat)
#' names(beta.stramash)
#' summary(beta.stramash)
#' head(as.data.frame(beta.stramash))
#' graphics::plot(betahat, beta.stramash$PosteriorMean, xlim=c(-4, 4), ylim=c(-4, 4))
#'
#' ## Testing the non-zero mode feature
#' betahat = betahat + 5
#' beta.stramash = stramash.workhorse(betahat = betahat, sebetahat = sebetahat)
#' graphics::plot(betahat, beta.stramash$PosteriorMean)
#' summary(beta.stramash)
#'
#' ## Running stramash with a pre-specified g, rather than estimating it
#' beta = c(rep(0, 100), stats::rnorm(100))
#' sebetahat = abs(stats::rnorm(200, 0, 1))
#' betahat = stats::rnorm(200, beta, sebetahat)
#' true_g = ashr::normalmix(c(0.5, 0.5), c(0, 0), c(0, 1)) # define true g
#' ## Passing this g into stramash causes it to
#' ## i) take the sd and the means for each component from this g, and
#' ## ii) initialize pi to the value from this g.
#' beta.stramash = stramash.workhorse(betahat = betahat, sebetahat = sebetahat,
#'                                    g = true_g, fixg = TRUE)
stramash.workhorse <- function(betahat, errordist = NULL, sebetahat = NULL,
                         likelihood = c("normal", "t"),
                         gridsize = 100,
                         method = c("fdr", "shrink"),
                         mixcompdist = c("uniform", "halfuniform", "normal", "+uniform", "-uniform"),
                         optmethod = c("mixIP", "cxxMixSquarem", "mixEM", "mixVBEM"),
                         df = NULL, randomstart = FALSE,
                         nullweight = 10, nonzeromode = FALSE,
                         pointmass = NULL,
                         prior = c("nullbiased", "uniform", "unit"),
                         mixsd = NULL, gridmult = sqrt(2),
                         outputlevel = 2, g = NULL, fixg = FALSE,
                         cxx = NULL, VB = NULL, model = c("EE", "ET"),
                         control = list()){

    likelihood  <- match.arg(likelihood)


    if (is.null(errordist) & is.null(sebetahat)) {
        stop("Either errordist or sebetahat needs to be specified")
    }
    if (likelihood == "t" & is.null(df)) {
        stop("If likelihood = \"t\" then df needs to be specified")
    }

    ## 1.Handling Input Parameters

    method      <- match.arg(method)
    mixcompdist <- match.arg(mixcompdist)
    optmethod   <- match.arg(optmethod)
    model       <- match.arg(model)


    assertthat::assert_that(is.null(errordist) | is.list(errordist))

    ## Various likelihoods are allowed
    if (likelihood == "t" & is.null(errordist) & !is.null(sebetahat)) {
        assertthat::assert_that(length(df) == 1 | length(df) == length(betahat))
        errordist <- list()
        if (length(df) == 1) {
            df <- rep(df, length = length(betahat))
        }
        for (index in 1:length(betahat)) {
            errordist[[index]] <- t_to_mix(mu = 0, sig = sebetahat[index],
                                           df = df[index], gridsize = 70)
        }
    } else if (likelihood == "normal" & is.null(errordist) & !is.null(sebetahat)) {
        errordist <- list()
        for (index in 1:length(betahat)) {
            errordist[[index]] <- normalmix(pi = 1, mean = 0, sd = sebetahat[index])
        }
    }


    assertthat::are_equal(length(betahat), length(errordist))
    class_vec <- sapply(errordist, class)
    assertthat::assert_that(all(class_vec == "normalmix" | class_vec == "unimix"))
    class_e <- unique(class_vec)

    ## set grid based on variance of errordist
    if (!is.null(sebetahat)) {
        message("Overwriting sebetahat because errordist is specified.")
    }

    sebetahat <- rep(NA, length = length(betahat))
    if (class_e == "normalmix") {
        for (seindex in 1:length(sebetahat)) {
            cdist <- errordist[[seindex]]
            second_moment <- sum(cdist$pi * (cdist$sd ^ 2 + cdist$mean ^ 2))
            first_moment2 <- sum(cdist$pi * cdist$mean) ^ 2
            sebetahat[seindex] <- sqrt(second_moment - first_moment2)
        }
    } else if (class_e == "unimix") {
        for (seindex in 1:length(sebetahat)) {
            cdist <- errordist[[seindex]]
            second_moment <- sum(cdist$pi *
                                 (cdist$a ^ 2 + cdist$a * cdist$b + cdist$b ^ 2)) / 3
            first_moment2 <- (sum(cdist$pi * (cdist$a + cdist$b)) / 2) ^ 2
            sebetahat[seindex] <- sqrt(second_moment - first_moment2)
        }
    }

    ## Make sure errordist and g arguments are correct
    if (nonzeromode) {
        stop("nonzeromode = TRUE not supported when errordist is not NULL")
    } else if (!is.null(g)) { # check where pointmass is
        if (class(g) == "normalmix") {
            which_pointmass <- abs(g$sd) < 10 ^ -14
            if (sum(pointmass) > 1) {
                stop("g may only have a pointmass at only one location")
            } else if (sum(pointmass) == 1) {
                if (abs(g$mean[which_pointmass]) > 10 ^ -14) {
                    stop("g may only have a pointmass at 0.")
                }
            }
        } else if (class(g) == "unimix") {
            which_pointmass <- abs(g$a - g$b) < 10 ^ -14
            if (sum(pointmass) > 1) {
                stop("g may only have a pointmass at only one location")
            } else if (sum(pointmass) == 1) {
                if (abs(g$a[which_pointmass]) > 10 ^ -14) {
                    stop("g may only have a pointmass at 0.")
                }
            }
        }
    } else if (outputlevel > 3) {
        stop("outputlevel > 3 not supported when errordist is not NULL")
    } else if (model == "ET") {
        stop("model = \"ET\" not supported when errordist is not NULL")
    }
    assertthat::are_equal(length(betahat), length(sebetahat))


    ## Capture all arguments into a list
    oldargs <- mget(names(formals()), sys.frame(sys.nframe()))
    newargs <- process_args(oldargs)

    ## Assign each argument in returned list to a variable used by the code next
    for (i in 1:length(newargs)) {
        assign(names(newargs)[i], newargs[[i]])
    }

    ## 2. Generating mixture distribution

    if (fixg & missing(g)) {
        stop("if fixg = TRUE then you must specify g!")
    }

    if (!is.null(g)) {
        k <- ncomp(g)
        null.comp <- 1 #null.comp not actually used unless randomstart true
        prior <- setprior(prior, k, nullweight, null.comp)
        if (randomstart) {
            pi <- initpi(k, n, null.comp, randomstart)
            g$pi <- pi
        } #if g specified, only initialize pi if randomstart is TRUE
    } else {
        if (is.null(mixsd)) {
            if (nonzeromode) {
                mixsd <- autoselect.mixsd(betahat[completeobs] - mean(betahat[completeobs]),
                                         sebetahat[completeobs],
                                         gridmult)
                if (pointmass) {
                    mixsd <- c(0, mixsd)
                }
                nonzeromode.fit <- nonzeromodeEM(betahat[completeobs],
                                                sebetahat[completeobs],
                                                mixsd = mixsd,
                                                mixcompdist = mixcompdist,
                                                df = df,
                                                control = controlinput)
                betahat[completeobs] <- betahat[completeobs] - nonzeromode.fit$nonzeromode
            }
            else if (nonzeromode & !is.null(df)) {
                ## stop("Error: Nonzero mean only implemented for df=NULL")
            }
            mixsd <- autoselect.mixsd(betahat[completeobs], sebetahat[completeobs], gridmult)
        }
        if (pointmass){
            mixsd <- c(0, mixsd)
        }


        null.comp <- which.min(mixsd) #which component is the "null"

        k <- length(mixsd)
        prior <- setprior(prior, k, nullweight, null.comp)
        pi <- initpi(k, n, null.comp, randomstart)


        if (!is.element(mixcompdist, c("normal", "uniform", "halfuniform", "+uniform", "-uniform"))) {
            stop("Error: invalid type of mixcompdist")
        }
        if (mixcompdist == "normal") {
            g <- normalmix(pi, rep(0, k), mixsd)
        } else if (mixcompdist == "uniform") {
            g <- unimix(pi, -mixsd, mixsd)
        } else if (mixcompdist == "+uniform") {
            g <- unimix(pi, rep(0, k), mixsd)
        } else if (mixcompdist == "-uniform") {
            g <- unimix(pi, -mixsd, rep(0, k))
        } else if (mixcompdist == "halfuniform") {
            if (min(mixsd) > 0) { #simply reflect the components
                g <- unimix(c(pi, pi) / 2, c(-mixsd, rep(0, k)), c(rep(0, k), mixsd))
                prior <- rep(prior, 2)
                pi <- rep(pi, 2)
            } else { #define two sets of components, but don't duplicate null component
                null.comp <- which.min(mixsd)
                g <- unimix(c(pi, pi[-null.comp]) / 2,
                            c(-mixsd, rep(0, k - 1)),
                            c(rep(0, k), mixsd[ - null.comp]))
                prior <- c(prior, prior[-null.comp])
                pi <- c(pi, pi[-null.comp])
            }
        }
    }

    ## check that all prior are >=1 (as otherwise have problems with infinite penalty)
    if (!all(prior >= 1) & optmethod != "mixVBEM") {
        stop("Error: prior must all be >= 1 (unless using optmethod mixVBEM)")
    }

    ## 3. Fitting the mixture
    if (!fixg){
        pi.fit <- estimate_mixprop(betahat[completeobs], g = g,
                                   prior = prior,
                                   null.comp = null.comp,
                                   optmethod = optmethod, df = df,
                                   control = controlinput,
                                   errordist = errordist)
    } else {
        pi.fit <- list(g = g)
    }

    ##4. Computing the posterior

    n <- length(betahat)

    ## specialized functions when likelihood is mixture
    postmixout    <- post_mix_dist(g = pi.fit$g, betahat = betahat,
                                   errordist = errordist)
    PosteriorMean <- mix_mean_array(mixdist = postmixout)
    PosteriorSD   <- mix_sd_array(mixdist = postmixout)
    ZeroProb      <- mix_probzero_array(mixdist = postmixout)
    NegativeProb  <- mix_cdf_array(mixdist = postmixout, q = 0) - ZeroProb
    lfsr          <- compute_lfsr(NegativeProb, ZeroProb)
    PositiveProb  <- 1 - NegativeProb - ZeroProb
    PositiveProb  <- ifelse(PositiveProb < 0, 0, PositiveProb)
    lfdr          <- ZeroProb
    qvalue        <- qval.from.lfdr(lfdr)
    svalue        <- qval.from.lfdr(lfsr)
    loglik        <- calc_loglik_array(g = pi.fit$g, betahat = betahat,
                                       errordist = errordist)
    logLR         <- loglik - calc_nulllik_array(betahat = betahat,
                                                 errordist = errordist)

    ## 5. Returning the result

    result <- list(fitted.g = pi.fit$g, call = match.call())
    if (outputlevel > 0) {
        result <- c(result, list(PosteriorMean = PosteriorMean,
                                 PosteriorSD = PosteriorSD,
                                 loglik = loglik,
                                 logLR = logLR))
    }
    if (outputlevel > 1) {
        result <- c(result, list(PositiveProb = PositiveProb,
                              NegativeProb = NegativeProb,
                              ZeroProb = ZeroProb,
                              lfsr = lfsr, lfdr = lfdr,
                              qvalue = qvalue,
                              svalue = svalue,
                              excludeindex = excludeindex,
                              model = model,
                              optmethod = optmethod))
    }
    if (outputlevel > 1.5){
        result <- c(result, list(data = list(betahat = betahat,
                                             sebetahat = sebetahat,
                                             df = df)))
    }
    if (outputlevel > 2) {
        result <- c(result, list(fit = pi.fit))
    }
    class(result) <- "ash"
    return(result)
}

#' Estimate mixture proportions of sigmaa by EM algorithm
#'
#' @param betahat (n vector of observations)
#' @param g the prior distribution for beta (usually estimated from
#'     the data
#' @param prior string, or numeric vector indicating Dirichlet prior
#'     on mixture proportions (defaults to "uniform", or (1,1...,1);
#'     also can be "nullbiased" (nullweight,1,...,1) to put more
#'     weight on first component)
#' @param null.comp the position of the null component
#' @param optmethod name of function to use to do optimization
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL(Gaussian)
#' @param control A list of control parameters for the SQUAREM
#'     algorithm, default value is set to be control.default=list(K =
#'     1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4,
#'     kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE).
#' @inheritParams stramash.workhorse
#' @return A list, including the final loglikelihood, the null
#'     loglikelihood, a n by k likelihoodmatrix with (j,k)th element
#'     equal to \eqn{f_k(x_j)},and a flag to indicate convergence.
#'
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)
#if cxx TRUE use cpp version of R function mixEM
estimate_mixprop <- function(betahat, g, prior,
                             optmethod = c("mixEM", "mixVBEM", "cxxMixSquarem", "mixIP"),
                             null.comp = 1, df = NULL,
                             control = list(), errordist = NULL){
    control.default <- list(K = 1, method = 3, square = TRUE, step.min0 = 1,
                            step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                            tol = 1.e-07, maxiter = 5000, trace = FALSE)
    optmethod <- match.arg(optmethod)
    namc <- names(control)
    if (!all(namc %in% names(control.default))) {
        stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
    }
    controlinput <- utils::modifyList(control.default, control)

    pi_init <- g$pi
    if (optmethod == "mixVBEM") {
        pi_init <- NULL
    }  #for some reason pi_init doesn't work with mixVBEM

    k <- ncomp(g)
    n <- length(betahat)
    controlinput$tol <- min(0.1 / n, 1.e-7) # set convergence criteria to be more stringent for larger samples

    if (controlinput$trace == TRUE) {
        tic()
    }

    matrix_llik <- t(log_compdens_conv_mix(g, betahat, errordist))

    assertthat::are_equal(nrow(matrix_llik), n)
    assertthat::are_equal(ncol(matrix_llik), k)
    matrix_llik <- matrix_llik - apply(matrix_llik, 1, max) #avoid numerical issues by subtracting max of each row
    matrix_lik <- exp(matrix_llik)

    ## the last of these conditions checks whether the gradient at the null is negative wrt pi0
    ## to avoid running the optimization when the global null (pi0=1) is the optimal.
    if (optmethod == "mixVBEM" || max(prior[-1]) > 1 || min(gradient(matrix_lik) + prior[1] - 1, na.rm = TRUE) < 0) {
        fit <- do.call(optmethod, args = list(matrix_lik = matrix_lik, prior = prior, pi_init = pi_init, control = controlinput))
    } else {
        fit <- list(converged = TRUE, pihat = c(1, rep(0, k - 1)))
    }

    ## check if IP method returns negative mixing proportions. If so, run EM.
    if (optmethod == "mixIP" & (min(fit$pihat) < -10 ^ -12)) {
        message("Interior point method returned negative mixing proportions.\n Switching to EM optimization.")
        optmethod <- "mixEM"
        fit <- do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                              prior = prior, pi_init = pi_init,
                                              control = controlinput))
    }

    if (!fit$converged) {
        warning("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
    }

    pi <- fit$pihat
    converged <- fit$converged
    niter <- fit$niter

    loglik.final <-  penloglik(pi, matrix_lik, 1) #compute penloglik without penalty
    null.loglik <- sum(log(matrix_lik[, null.comp]))
    g$pi <- pi
    if (controlinput$trace == TRUE) {
        toc()
    }

    return(list(loglik = loglik.final, null.loglik = null.loglik,
                matrix_lik = matrix_lik, converged = converged, g = g,
                niter = niter))
}

#' This is a copy from the ashr package because it is not exported and
#' I need it.
#'
#' @inheritParams stramash
#' @param mult The amount to multiply by.
autoselect.mixsd <- function (betahat, sebetahat, mult) {
    sebetahat <- sebetahat[sebetahat != 0]
    sigmaamin <- min(sebetahat) / 10
    if (all(betahat ^ 2 <= sebetahat ^ 2)) {
        sigmaamax <- 8 * sigmaamin
    } else {
        sigmaamax <- 2 * sqrt(max(betahat ^ 2 - sebetahat ^ 2))
    }
    if (mult == 0) {
        return(c(0, sigmaamax / 2))
    } else {
        npoint <- ceiling(log2(sigmaamax / sigmaamin) / log2(mult))
        return(mult ^ ( (-npoint):0) * sigmaamax)
    }
}

#' This is a copy from the ashr package because it is not exported and I need it.
#'
#' @param prior Either \code{"nullbiased"}, \code{"uniform"}, or
#'     \code{"unit"}.
#' @param k The size of the grid.
#' @param nullweight The amount to penalize the null.
#' @param null.comp The position of the null component.
setprior <- function (prior, k, nullweight, null.comp) {
    if (!is.numeric(prior)) {
        if (prior == "nullbiased") {
            prior <- rep(1, k)
            prior[null.comp] <- nullweight
        } else if (prior == "uniform") {
            prior <- rep(1, k)
        } else if (prior == "unit") {
            prior <- rep(1 / k, k)
        }
    }
    if (length(prior) != k | !is.numeric(prior)) {
        stop("invalid prior specification")
    }
    return(prior)
}

#' This is copied from the ashr package because I need it and it is
#' not exported.
#'
#' @param k The size of the grid.
#' @param n The sample size.
#' @param null.comp The position of the null component.
#' @param randomstart A logical. Should we start at a random location?
initpi <- function (k, n, null.comp, randomstart) {
    if (randomstart) {
        pi <- stats::rgamma(k, 1, 1)
    } else {
        if (k < n) {
            pi <- rep(1, k) / n
            pi[null.comp] <- (n - k + 1) / n
        } else {
            pi <- rep(1, k) / k
        }
    }
    pi <- normalize(pi)
    return(pi)
}

#' Copied from ashr package because I need it and it is not exported.
#'
#' @param x A vector of numerics.
normalize <- function (x) {
    return(x / sum(x))
}

#' Copied from ashr package because I need it and it is not exported.
#'
#' @param matrix_lik The \eqn{(i, j)}th element of \code{matrix_like}
#'     is the contribution to the likelihood of the \eqn{i}th
#'     observation and the \eqn{j}th component.
#'
#' @return The kth element of this vector is the derivative
#'     of the loglik for \eqn{\pi=(\pi_0,...,1-\pi_0,...)} with respect to
#'     \eqn{\pi_0} at \eqn{\pi_0 = 1}.
gradient <- function (matrix_lik) {
    n <- nrow(matrix_lik)
    grad <- n - colSums(matrix_lik / matrix_lik[, 1])
    return(grad)
}

#' Copied from ashr package because I need it and it is not exported.
#'
#' @param pi A vector of porportions that sum to 1.
#' @param prior A vector of weights of the same length as \code{pi}.
#' @param matrix_lik The \eqn{(i, j)}th element of \code{matrix_like}
#'     is the contribution to the likelihood of the \eqn{i}th
#'     observation and the \eqn{j}th component.
penloglik <- function (pi, matrix_lik, prior) {
    pi <- normalize(pmax(0, pi))
    m <- t(pi * t(matrix_lik))
    m.rowsum <- rowSums(m)
    loglik <- sum(log(m.rowsum))
    subset <- (prior != 1)
    priordens <- sum( (prior - 1)[subset] * log(pi[subset]))
    return(loglik + priordens)
}

#' Function to compute the local false sign rate.
#'
#' Copied from the ashr package.
#'
#' @param NegativeProb A vector of posterior probability that beta is
#'     negative.
#' @param ZeroProb A vector of posterior probability that beta is
#'     zero.
#' @return The local false sign rate.
compute_lfsr <- function (NegativeProb, ZeroProb) {
    ifelse(NegativeProb > 0.5 * (1 - ZeroProb), 1 - NegativeProb,
        NegativeProb + ZeroProb)
}


#' From the ashr package because not exported but I need it.
#'
#' @param gcFirst A loogical.
#' @param type One of \code{"elapsed"}, \code{"user.self"}, or
#'     \code{"sys.self"}.
tic <- function (gcFirst = TRUE, type = c("elapsed", "user.self", "sys.self")) {
    type <- match.arg(type)
    assign(".type", type, envir = baseenv())
    if (gcFirst)
        gc(FALSE)
    tic <- proc.time()[type]
    assign(".tic", tic, envir = baseenv())
    invisible(tic)
}

#' From teh ashr package because not exported but I need it.
#'
#'
toc <- function() {
  type <- get(".type", envir = baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir = baseenv())
  print(toc - tic)
  invisible(toc)
}
