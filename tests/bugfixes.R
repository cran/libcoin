
library("libcoin")

### by Henric Winell
X <- runif(10)
Y <- runif(10)
o <- LinStatExpCov(X, Y)
ov <- LinStatExpCov(X, Y, varonly = TRUE)
stopifnot(all.equal(doTest(o, teststat = "scalar"),
                    doTest(ov, teststat = "scalar")))

### all weights = 0 and no weights at all was treated the same
X <- as.double(1:10)
Y <- as.double(10:1)
sum(X*Y)
cl <- gl(2, 5)

### linstat = 220
w <- as.integer(rep(1, 10))
LinStatExpCov(X = X, Y = Y)
LinStatExpCov(X = X, Y = Y, weights = w)
LinStatExpCov(X = X, Y = Y, weights = w, block = cl)

### linstat = 0
w <- as.integer(rep(0, 10))
LinStatExpCov(X = X, Y = Y, weights = w)
LinStatExpCov(X = X, Y = Y, weights = w, block = cl)

### linstat = 110
w <- as.integer(rep(0, 10))
w[1:5] <- 1L
LinStatExpCov(X = X, Y = Y, subset = 1:5)
LinStatExpCov(X = X, Y = Y, weights = w)
LinStatExpCov(X = X, Y = Y, weights = w, block = cl)

### linstat = 0
LinStatExpCov(X = X, Y = Y, weights = w, subset = 6:10)
LinStatExpCov(X = X, Y = Y, weights = w, block = cl, subset = 6:10)
