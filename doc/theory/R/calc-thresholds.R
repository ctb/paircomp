threshold <- function(seqsize, windowsize) {
   v <- PLmax(seqsize, c(1:windowsize), windowsize, 1/4)

   for (i in windowsize:1) {
      if (v[i] > 0.05) break
   }
   i
}

seqsizes <- c(1e3, 1e4, 1e5, 1e6)
ws <- c(2:10) * 5

for (i in seqsizes) {
   for (j in ws) {
      th <- threshold(i, j)
      print(c(i, j, th))
   }
}
