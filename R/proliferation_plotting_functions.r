

#-------------------------------------------------------------------------------
#' Plot optimization of parameters
#'
opt.plot.final <- function(x.mids, y, y.pred.optim) {
  plot(x.mids, y,
       main="Optimal fit vs data",
       xlab="Bin average log10(CTV)",
       ylab="Count",
       pch=20,
       col="grey", bty="n")
  lines(y.pred.optim ~ x.mids, lwd=2, col="red")

  legend("topleft",
         legend=c("Data", "Optimized fit"),
         fill=c("grey", "red"),
         bty="n")
}

#-------------------------------------------------------------------------------
#' Plot optimization of parameters
#'
#'
opt.plot.params <- function(x, opt.env) {

  n     <- ncol(opt.env[["means"]])
  cols  <- palette.colors(n = n, recycle = T)

  #--------------------------------------------------------------------------
  # Score
  plot(1:length(opt.env[["score"]]),
       opt.env[["score"]],
       xlab="Iteration",
       ylab="NLL / RSS",
       pch=20,
       col="grey",
       bty="n",
       main="Neg. Log Lik. / Res. Sum of Sqr.",
       type="l")

  #--------------------------------------------------------------------------
  n.iter <- nrow(opt.env[["means"]])
  # Mean
  plot(1:n.iter,
       opt.env[["means"]][,1],
       xlab="Iteration",
       ylab="Mean",
       pch=20,
       col=cols[1],
       bty="n",
       main="Optim - mean",
       type="l",
       ylim=c(min(x), max(x)))

  for (i in 2:ncol(opt.env[["means"]])) {
    lines(1:n.iter,
          opt.env[["means"]][,i],
          col=cols[i])
  }

  legend("right",
         legend=paste0("Gen", 0:(n-1)),
         fill=cols,
         bty="n")
  #--------------------------------------------------------------------------
  # SD
  plot(1:n.iter,
       opt.env[["sd"]][,1],
       xlab="Iteration",
       ylab="Peak SD",
       pch=20,
       col=cols[i],
       bty="n",
       main="Optim - SD",
       type="l",
       ylim=c(min(opt.env[["sd"]]), max(opt.env[["sd"]])))

  for (i in 2:ncol(opt.env[["sd"]])) {
    lines(1:n.iter,
          opt.env[["sd"]][,i],
          col=cols[i])
  }

  legend("right",
         legend=paste0("Gen", 0:(n-1)),
         fill=cols,
         bty="n")
  #--------------------------------------------------------------------------
  # Summit
  plot(1:n.iter,
       opt.env[["summits"]][,1],
       xlab="Iteration",
       ylab="Summit",
       pch=20,
       col=cols[1],
       bty="n",
       main="Optim - summit",
       type="l",
       ylim=c(min(opt.env[["summits"]]), max(opt.env[["summits"]])))

  for (i in 2:ncol(opt.env[["summits"]])) {
    lines(1:n.iter,
          opt.env[["summits"]][,i],
          col=cols[i])
  }

  legend("right",
         legend=paste0("Gen", 0:(n-1)),
         fill=cols,
         bty="n")
}
