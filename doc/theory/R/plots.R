postscript(file='trials.ps')
L <- c(0:1e3) * 1e3

r <- range(0, 1)
par(new=FALSE)
plot(L, PL(L, y[16]), "l", ylim=r)
par(new=TRUE)
plot(L, PL(L, y[17]), "l", ylim=r)
par(new=TRUE)
plot(L, PL(L, y[18]), "l", ylim=r)
par(new=TRUE)
plot(L, PL(L, y[19]), "l", ylim=r)
par(new=TRUE)
plot(L, PL(L, y[20]), "l", ylim=r)
par(new=TRUE)
plot(L, PL(L, y[21]), "l", ylim=r)

dev.off()
