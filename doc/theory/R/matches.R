# probability of getting a match
p <- 1/4

# number of balls
N <- 20

# vector of balls
x <- c(0:N)

# probability distribution for matches
y <- dbinom(x, N, p)

# function to calculate probability of seeing a particular match in L tries.
PL <- function(L, binp) {
    1 - exp(-L * binp)
    }

# Generates probability that in time L we see M matches == best when picking
# from N with probability p of a match.  M is a vector.
PLmax <- function(L, M, N, p) {
  total <- 0;
  result <- c(1:length(M)) * 0;
  for (i in 1:length(M)) {
     m <- M[i];
     s <- 0
     if ((m + 1) <= N) {
        x <- c((m + 1):N);
        y <- dbinom(x, N, p);
        s <- PL(L, sum(y))
     }
     t <- PL(L, dbinom(m, N, p))
     result[i] <- t - t * s
  }
  result
}

# generates probability distribution against sequence length L (vector) of
# best matches, picking N bp with probability p of match.
LEmax <- function(L, N, p) {
   x <- c(0:N)
   result <- c(1:length(L)) * 0
   
   for (i in 1:length(L)) {
      l <- L[i]

      result[i] <- sum(x * PLmax(l, x, N, p))
   }
   result
}

#
# Output:
#

postscript(file='bestmatches-20.ps')
plot(x, PLmax(1e12, c(0:20), 20, 1/4), type="b")
dev.off()

postscript(file='bestmatches-50.ps')
plot(c(0:50), PLmax(1e12, c(0:50), 50, 1/4), type="b")
dev.off()

# t: 0 - 10billion, spaced 1000 apart.
t <- c(0:1e3) * 1e7

dist20 <- LEmax(t, 20, 1/4)
dist50 <- LEmax(t, 50, 1/4)

percent20 <- dist20 * 5
percent50 <- dist50 * 2

postscript(file='window=20.ps')
plot(t, percent20, ylim=range(0, 100))
dev.off()

postscript(file='window=50.ps')
plot(t, percent50, ylim=range(0, 100))
dev.off()
