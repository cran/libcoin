
ctabs <- function(ix, iy = integer(0), weights = integer(0),
                  subset = integer(0), block = integer(0))
{

    if (is.null(attr(ix, "levels")))
            attr(ix, "levels") <- 1:max(ix)

    if (length(iy) > 0) {
        if (is.null(attr(iy, "levels")))
            attr(iy, "levels") <- 1:max(iy)
    }

    if (length(subset) > 0) subset <- subset - 1L

    ret <- .Call(R_tables, ix, iy, weights, subset, block)

    if (length(block) > 0) {
        if (length(iy) == 0)
            ret <- ret[,,-1, drop = FALSE]
        else
            ret <- ret[,,-dim(ret)[3], drop = FALSE]
    } else {
        ret <- ret[,,,drop = TRUE]
    }
    ret
}
