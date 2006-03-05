P <- function(p) 1 - 3/4 * p

# pv vector, L, N constants.  pv elements are probabilities of a *match*.
DIVEmax <- function(L, N, pv) {
   x <- c(0:N)
   result <- c(1:length(pv)) * 0
   
   for (i in 1:length(pv)) {
      p <- pv[i]

      result[i] <- sum(x * PLmax(L, x, N, P(p)))
   }
   result
}

p <- c(0:100) / 100;

postscript(file='divergedist-50.ps', onefile=FALSE)
percentdiv <- DIVEmax(1e12, 50, p) * 2
plot(p, percentdiv, type="l")
dev.off()

postscript(file='divergedist-20.ps', onefile=FALSE)
percentdiv <- DIVEmax(1e12, 20, p) * 5
plot(p, percentdiv, type="l")
dev.off()
