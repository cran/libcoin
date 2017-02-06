
.LinStatExpCov1d <- function(X, Y, weights = integer(0), subset = integer(0), block = integer(0),
                             varonly = FALSE, B = 0L, standardise = FALSE,
                             tol = sqrt(.Machine$double.eps))
{

    if (NROW(X) != NROW(Y))
        stop("dimensions of X and Y don't match")

    if (is.integer(X)) {
        if (is.null(attr(X, "levels")))
            attr(X, "levels") <- 1:max(X)
    }

    if (length(weights) > 0) {
        if (!((NROW(X) == length(weights)) &&
              is.integer(weights) &&
              all(weights >= 0)))
            stop("incorrect weights")
    }

    if (length(subset) > 0) {
        rs <- range(subset)
        if (!((rs[2] <= NROW(X)) &&
              (rs[1] >= 1L) &&
              is.integer(subset)))
            stop("incorrect subset")
        if (rs[1] == 0) stop("subset has start 1 index")
        subset <- subset - 1L
    }

    if (length(block) > 0) {
        if (!((NROW(X) == length(block)) &&
              is.factor(block)))
            stop("incorrect block")
    }

    ms <- !(complete.cases(X) & complete.cases(Y))
    if (all(ms))
        stop("all observations are missing")
    if (any(ms)) {
        if (length(subset) > 0) {
            if (all(subset %in% (which(ms) - 1L)))
                stop("all observations are missing")
            subset <- subset[!(subset %in% (which(ms) - 1L))]
        } else {
            subset <- (0:(NROW(X) - 1))[-which(ms)]
        }
    }

    ret <- .Call(R_ExpectationCovarianceStatistic, X, Y, weights, subset,
                 block, as.integer(varonly), as.double(tol))
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (B > 0)
        ret$PermutedLinearStatistic <-
            .Call(R_PermutedLinearStatistic, ret, X, Y, weights, subset,
                  block, as.integer(B), as.integer(standardise))
    ret
}

.LinStatExpCov2d <- function(X = numeric(0), Y, ix, iy, weights = integer(0), subset = integer(0),
                             block = integer(0), varonly = FALSE, B = 0,
                             standardise = FALSE,
                             tol = sqrt(.Machine$double.eps))
{
    if (!((length(ix) == length(iy)) &&
          is.integer(ix) && is.integer(iy)))
        stop("incorrect ix and/or iy")

    if (is.null(attr(ix, "levels")))
        attr(ix, "levels") <- 1:max(ix)
    if (is.null(attr(iy, "levels")))
        attr(iy, "levels") <- 1:max(iy)

    if (length(X) > 0) {
        if (!((min(ix) >= 0 && nrow(X) == (length(attr(ix, "levels")) + 1)) &&
              all(complete.cases(X)) &&
              (nrow(X) == (length(attr(ix, "levels")) + 1))))
            stop("incorrect X")
    }

    if (!(all(complete.cases(Y))) &&
          (nrow(Y) == (length(attr(iy, "levels")) + 1)) &&
          (min(iy) >= 0L && nrow(Y) == (length(attr(iy, "levels")) + 1)))
        stop("incorrect Y")

    if (length(weights) > 0) {
        if (!((length(ix) == length(weights)) &&
              is.integer(weights) &&
              all(weights >= 0)))
            stop("incorrect weights")
    }

    if (length(subset) > 0) {
        rs <- range(subset)
        if (!((rs[2] <= length(ix)) &&
              (rs[1] >= 1L) &&
              is.integer(subset)))
            stop("incorrect subset")
        if (rs[1] == 0) stop("subset has start 1 index")
        subset <- subset - 1L
    }

    if (!missing(block) && length(block) > 0) {
        if (!((length(ix) == length(block)) &&
              is.factor(block)))
            stop("incorrect block")
    }

    ret <- .Call(R_ExpectationCovarianceStatistic_2d, X, ix, Y, iy,
                 weights, subset, block, as.integer(varonly), as.double(tol))
    ret$varonly <- as.logical(ret$varonly)
    ret$Xfactor <- as.logical(ret$Xfactor)
    if (B > 0)
        ret$PermutedLinearStatistic <-
            .Call(R_PermutedLinearStatistic_2d, ret, X, ix, Y, iy,
                  block, as.integer(B), as.integer(standardise))
    ret
}

LinStatExpCov <- function(X, Y, ix = NULL, iy = NULL, weights = integer(0),
                          subset = integer(0), block = integer(0),
                          varonly = FALSE, B = 0, standardise = FALSE,
                          tol = sqrt(.Machine$double.eps))
{

    if (is.null(ix) & is.null(iy))
        return(.LinStatExpCov1d(X = X, Y = Y, weights = weights,
                                subset = subset, block = block,
                                varonly = varonly, B = B,
                                standardise = standardise, tol = tol))

    if (!is.null(ix) & !is.null(iy))
        return(.LinStatExpCov2d(X = X, Y = Y, ix = ix, iy = iy,
                                weights = weights, subset = subset,
                                block = block, varonly = varonly, B = B,
                                standardise = standardise, tol = tol))

    if (missing(X) & !is.null(ix))
        return(.LinStatExpCov1d(X = ix, Y = Y, weights = weights,
                                subset = subset, block = block,
                                varonly = varonly, B = B,
                                standardise = standardise, tol = tol))

    stop("incorrect call to LinStatExpCov")
}


### note: lower = FALSE => p-value; lower = TRUE => 1 - p-value
doTest <- function(object, teststat = c("maximum", "quadratic", "scalar"),
                   alternative = c("two.sided", "less", "greater"),
                   pvalue = TRUE, lower = FALSE, log = FALSE,
                   minbucket = 10L, ordered = TRUE, pargs = GenzBretz())
{

    ### avoid match.arg for performance reasons
    teststat <- teststat[1]
    if (!any(teststat == c("maximum", "quadratic", "scalar")))
        stop("incorrect teststat")
    alternative <- alternative[1]
    if (!any(alternative == c("two.sided", "less", "greater")))
        stop("incorrect alternative")

    if (teststat == "quadratic") {
        if (alternative != "two.sided")
            stop("incorrect alternative")
    }

    test <- which(c("maximum", "quadratic", "scalar") == teststat)
    if (test == 3) {
        if (length(object$LinearStatistic) != 1)
            stop("scalar test statistic not applicable")
        test <- 1L ### scalar is maximum internally
    }
    alt <- which(c("two.sided", "less", "greater") == alternative)

    if (!pvalue & (NCOL(object$PermutedLinearStatistic) > 0))
        object$PermutedLinearStatistic <- matrix(NA_real_, nrow = 0, ncol = 0)

    if (!object$Xfactor) {
        if (teststat == "quadratic") {
            ret <- .Call(R_QuadraticTest, object,
                         as.integer(pvalue), as.integer(lower), as.integer(log))
        } else {
            ret <- .Call(R_MaximumTest, object,
                         as.integer(alt), as.integer(pvalue), as.integer(lower),
                         as.integer(log), as.integer(pargs$maxpts),
                         as.double(pargs$releps), as.double(pargs$abseps))
            if (teststat == "scalar") {
                var <- if (object$varonly) object$Variance else object$Covariance
                ret$TestStatistic <- object$LinearStatistic - object$Expectation
                ret$TestStatistic <-
                    if (var > object$tol) ret$TestStatistic / sqrt(var) else NaN
            }
        }
    } else {
        ret <- .Call(R_MaximallySelectedTest, object, as.integer(ordered),
                     as.integer(test), as.integer(minbucket),
                     as.integer(lower), as.integer(log))
    }
    ret
}
