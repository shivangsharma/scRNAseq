calc_CI <- function(vector, alpha) {

  se <- sd(vector)/sqrt(length(vector))
  t <- qt((1 - alpha)/2 + 0.5, length(vector) - 1)
  CI = t*se
  CI

}

plot_density <- function(x, y, w) {

  dens <- ks::kde(x = matrix(data = c(x, y), ncol = 2), w = w / sum(w) * length(w), bgridsize = rep(1000, 2))
  ix <- findInterval(x, dens$eval.points[[1]])
  iy <- findInterval(y, dens$eval.points[[2]])
  ii <- cbind(ix, iy)
  z <- dens$estimate[ii]
  names(z) <- colnames(obj)
  z

}

