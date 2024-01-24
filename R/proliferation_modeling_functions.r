#-------------------------------------------------------------------------------
#' Find Local Max
#'
#' Find the maximum value of y within a range of x values
#'
#' @param x vector of numbers matching y in length.
#' @param y vector of numbers matching x in length.
#' @param lower.x lower limit of x to search.
#' @param upper.x upper limit of x to search.
#' @returns the value of x where y is maximised.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' find.local.max(x, y, -1, 1)
find.local.max <- function(x, y, lower.x, upper.x){
  y.f <- y[x >= lower.x & x <= upper.x]
  x.f <- x[x >= lower.x & x <= upper.x]
  ymax <- max(y.f)
  return(mean(x.f[y.f==ymax]))
}

#-------------------------------------------------------------------------------
#' Find Local Min
#'
#' Find the minimum value y of x within a range of x values
#'
#' @param x vector of numbers matching y in length.
#' @param y vector of numbers matching x in length.
#' @param lower.x lower limit of x to search.
#' @param upper.x upper limit of x to search.
#' @returns the value of x where y is minimized
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' find.local.min(x, y, -1, 1)
find.local.min <- function(x, y, lower.x, upper.x){
  y.f <- y[x >= lower.x & x <= upper.x]
  x.f <- x[x >= lower.x & x <= upper.x]
  ymin <- min(y.f)
  return(mean(x.f[y.f==ymin]))
}

#-------------------------------------------------------------------------------
#' Nearest Neighbour smoother
#'
#' Smooth a trace using a kernel of size window.
#'
#' @param y a vector y
#' @param window how many values up and down are averaged
#' @returns the mean of each y
#' value and the surrounding values specified by window
#'
#' @examples
#' y <- rnorm(100)
#' nn.smoother(y, 2)
#' @export
nn.smoother <- function(y, window=2) {
  # y <- y[order(y, decreasing = T)]
  y.out <- c()

  for(i in 1:length(y)) {
    i.min <- i-window
    i.min <- ifelse(i.min < 0, 0, i.min)
    i.max <- i+window
    i.max <- ifelse(i.max > length(y), length(y), i.max)

    y.out[i] <- mean(y[i.min:i.max])
  }
  return(y.out)
}

#-------------------------------------------------------------------------------
#' Find nearest index
#'
#' Find the index of the value in x that is closest to value
#'
#' @param x a vector
#' @param the value to search the nearest value in x
#' @returns the index in x closest to value
#'
#' @examples
#' x <- rnorm(100)
#' nearest.index(x, 0)
nearest.index <- function(x, value) {
  which.min(abs(x-value))
}

#-------------------------------------------------------------------------------
#' Find the mode of the next peak using the previous peak and the average peak
#' distance
find.next.peak <- function(x, y, prev.peak, peak.dist) {

  est.window      <- peak.dist/2
  est.curpeak     <- log10((10^prev.peak)/2)
  curpeak.mean    <- find.local.max(x, y, est.curpeak-est.window, est.curpeak+est.window)
  curpeak.summit  <- y[nearest.index(x, curpeak.mean)]

  # calculate enrichment
  fold.change  <- find.enrichment(x, y, curpeak.mean, est.window)

  return(c(curpeak.mean, curpeak.summit, fold.change))
}

#-------------------------------------------------------------------------------
#' Find the relative enrichment over a valley
find.enrichment <- function(x, y, curpeak.mean, est.window) {

  curpeak.summit  <- y[nearest.index(x, curpeak.mean)]

  # calculate enrichment
  a <- nearest.index(x, curpeak.mean + est.window)
  b <- nearest.index(x, curpeak.mean - est.window)
  valley.count <- mean(y[c(a, b)])
  fold.change  <- curpeak.summit/valley.count

  return(fold.change)
}

#-------------------------------------------------------------------------------
#' Scaled distribution of a single peak in prolif model
#'
#' @param x x values (log scaled intensity values)
#' @param mean the mean of the peak
#' @param sd the sd of the peak
#' @param the peak height (estimated relative to the other peaks)
#' @returns scaled density for peak height
prolif.single.peak <- function(x, mean, sd, summit) {
  # Gaussian density
  cur.density <- dnorm(x, mean=mean, sd=sd)

  # Density so the density on the mode / mean is equal to one
  cur.density <- cur.density / dnorm(mean, mean=mean, sd=sd)

  # Scale by the summit
  return(cur.density*summit)
}

#-------------------------------------------------------------------------------
#' Mixture distribution of Gaussian
#'
#' Gaussian prolif model, it is a mixture distribution scaled by relative
#' peak heights.
prolif.model <- function(means, sd, summits, x) {

  y.pred <- 0

  for (i in 1:length(means)) {
    cur.pred <- prolif.single.peak(x, means[i], sd[i], summits[i])
    y.pred   <- y.pred + cur.pred
  }

  return(y.pred)
}

#-------------------------------------------------------------------------------
#' Density of mixture distribution of Gaussian
#'
prolif.model.density <- function(means, sd, summits, x, verbose=F) {
  density <- 0

  summits <- summits/sum(summits)

  if (verbose) {
    cat("[INFO] means: ",   means, "\n")
    cat("[INFO] sd: ",      sd, "\n")
    cat("[INFO] summits: ", summits, "\n")
    cat("[INFO] --------------------------------------------- \n")

  }

  density <- 0
  for (i in 1:length(means)) {
    cur.d   <- dnorm(x, mean=means[i], sd=sd[i])
    density  <- density + (summits[i] * cur.d)
  }

  # Hack to deal with numeric precision limit on very small densities
  if (sum(density==0) >0) {
    msg <- "Zero densities detected. This is due to a numeric precision limit in dnorm().\nAs a hack setting zeroes to  4.761311e-320. This does invalidate the PDF"
    warning(simpleWarning(msg))
    density[density==0] <- 4.761311e-320
  }

  return(density)
}

#-------------------------------------------------------------------------------
#' Negative log likelihood of mixture distribution of Gaussian
#'
prolif.model.nll <- function(means, sd, summits, x, verbose=F) {
  density <- prolif.model.density(means, sd, summits, x, verbose=verbose)

  nll <- -sum(log(density))

  if (verbose) {
    cat("[INFO] nll: ", nll, "\n")
  }
  return(nll)
}

#-------------------------------------------------------------------------------
#' Wrapper that returns a gaussian proliferation model.
#' Takes a vector of paremeters par, number of peaks, and data X
prolif.model.wrapper <- function(par, n.peaks, x, fixed=NULL, type="prolif_model", verbose=F, opt.env=NULL) {

  names      <- names(par)
  par        <- as.numeric(par)
  names(par) <- names


  mean.names   <- paste0("gen", 0:(n.peaks-1), ".mean")
  summit.names <- paste0("gen", 0:(n.peaks-1), ".summit")

  if (!is.null(fixed)) {
    par <- c(par, fixed)
  }

  mean.vec   <- as.numeric(par[mean.names])
  summit.vec <- as.numeric(par[summit.names])
  sd.vec     <- as.numeric(par["peak.sd"])

  if (length(sd.vec)==1) {
    sd.vec <- rep(sd.vec, length(mean.vec))
  }

  # Log the values to the optimization environment
  if (!is.null(opt.env)) {
    opt.env[["means"]]   <- rbind(opt.env[["means"]], mean.vec)
    opt.env[["sd"]]      <- rbind(opt.env[["sd"]], sd.vec)
    opt.env[["summits"]] <- rbind(opt.env[["summits"]], summit.vec)
  }

  if (type == "prolif_model") {
    return(prolif.model(mean.vec, sd.vec, summit.vec, x))

  } else if (type == "density") {
    return(prolif.model.density(mean.vec, sd.vec, summit.vec, x, verbose=verbose))

  } else if (type == "neg_log_likelihood") {
    score <- prolif.model.nll(mean.vec, sd.vec, summit.vec, x, verbose=verbose)
    if (!is.null(opt.env)) {
      opt.env[["score"]]   <- c(opt.env[["score"]], score)
    }
    return(score)
  } else {

    e <- simpleError("Invalid type argument.")
    stop(e)
  }
}

#-------------------------------------------------------------------------------
#' Residuals between observed CTV trace and Gaussian prolif model
prolif.resid <- function(par, n.peaks, x, y, fixed=NULL, opt.env=NULL) {
  residuals <- y - prolif.model.wrapper(par, n.peaks, x, fixed, type="prolif_model", opt.env=opt.env)
  return(residuals)
}

#-------------------------------------------------------------------------------
#' Estimate initial parameters for proliferation model
#'
#' @param cur.hist histogram object
#' @param y a smoothed or raw count trace
#' @param peak.thresh.enrich fold change over valley to call peak in initial estimation (default 1)
#' @param peak.thresh.summit minimum height of a peak in percentage of total heights (default 0.05)
#' @param peak.max the maximum number of peaks to search for (default 12)
#' @param plot should plot be generated
#'
#' @returns data frame with initial peak estimates

find.initial.peaks <- function(cur.hist, y, peak.0.lower.bound,peak.thresh.enrich=1, peak.thresh.summit=0.05, peak.max=12, plot=F) {
  peak0.mean       <- find.local.max(cur.hist$mids, y, peak.0.lower.bound, max(cur.hist$mids))

  if (is.na(peak0.mean)) {
    cat("[WARN] No values found in peak.0.lower.bound. Returning NA.\n")
    return(NULL)
  }

  peak0.summit     <- y[nearest.index(cur.hist$mids, peak0.mean)]

  # Estimate position of next peak
  est.peak1     <- log10(10^peak0.mean/2)
  est.peak.dist <- peak0.mean - est.peak1
  #est.window    <- est.peak.dist/2
  peak0.enrichment <- find.enrichment(cur.hist$mids, y, peak0.mean, est.peak.dist/2)

  prev.peak.mean <- peak0.mean
  peak.stats     <- data.frame()
  lim <- min(cur.hist$mids)

  for (i in 1:peak.max) {

    if (prev.peak.mean-est.peak.dist < lim) {
      next;
    }
    res            <- find.next.peak(cur.hist$mids, y, prev.peak.mean, est.peak.dist/2)
    peak.stats     <- rbind(peak.stats, res)
    prev.peak.mean <- res[1]

  }

  # Add gen 0 to the output table
  peak.stats            <- rbind(c(peak0.mean,
                                   peak0.summit,
                                   peak0.enrichment),
                                 peak.stats)
  colnames(peak.stats)  <- c("est_mean", "est_summit", "est_enrichment")
  peak.stats$index      <- 1:(nrow(peak.stats))
  peak.stats$generation <- 0:(nrow(peak.stats)-1)

  peak.stats$est_summit_percentage <- peak.stats$est_summit/sum(peak.stats$est_summit)

  # Filter peaks based on summit height, while taking the next peak into account.
  # It can happen the earlier generations are low in height due to high proliferation.
  j <- 1
  for (i in 2:nrow(peak.stats)) {

    if (peak.stats[i, "est_summit_percentage"] > peak.thresh.summit) {
      j <- j+1
      next
    } else {
      if (peak.stats[min(i+1, nrow(peak.stats)), "est_summit_percentage"] > peak.thresh.summit) {
        j <- j+1
        next
      } else {
        break
      }
    }
  }
  peak.stats <- peak.stats[1:j,]

  # Filter peaks on threshold, always keep gen 0. This can be filtered if it is really small
  #peak.stats <- peak.stats[(peak.stats$est_enrichment > peak.thresh.enrich) & (peak.stats$est_summit_percentage  > peak.thresh.summit) ,, drop=F]
  peak.stats <- peak.stats[(peak.stats$est_enrichment > peak.thresh.enrich) | peak.stats$generation == 0,, drop=F]

  # Remove peaks that are not sequential in order
  n.peaks  <- 0
  for (i in 1:(nrow(peak.stats)-1)) {
    cur.idx <- peak.stats[i, "index"]

    if (cur.idx==1) {
      n.peaks  <- 1
    }

    if (i+1 > nrow(peak.stats)) {
      break
    }

    if (cur.idx +1 == peak.stats[i+1, "index"]){
      n.peaks  <- n.peaks +1
    }  else {
      break
    }
  }

  # Always force starting generation to be present
  if (n.peaks==0) {
    n.peaks <- 1
  }

  peak.stats <- peak.stats[1:n.peaks,, drop=F]

  #####
  # Plot
  if (plot) {
    plot(cur.hist$mids,
         cur.hist$counts,
         main="Initial peak estimates",
         ylab="Count",
         xlab="Bin average log10(CTV)",
         bty="n",
         pch=20,
         col="grey")
    lines(cur.hist$mids, y, new=F, lwd=3, col="blue")
    abline(v=peak0.mean, col="red", lwd=2)
    abline(v=peak.0.lower.bound, col="red", lty=2)

    for (i in 2:nrow(peak.stats)) {
      abline(v=peak.stats[i, 1], lwd=2, col="grey")
    }

    legend("topleft",
           legend=c("Smoothed data", "Gen0"),
           fill=c("blue", "red"),
           bty="n")

  }
  #####

  return(peak.stats)
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
       type="l")

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

#-------------------------------------------------------------------------------
#' Create a new env for storing optimization trace
#'
#'
opt.new.env <- function(n.peaks) {
  opt.env              <- new.env()
  opt.env[["score"]]   <- c()
  opt.env[["means"]]   <- matrix(ncol=n.peaks, nrow=0)
  opt.env[["sd"]]      <- matrix(ncol=n.peaks, nrow=0)
  opt.env[["summits"]] <- matrix(ncol=n.peaks, nrow=0)

  return(opt.env)
}

#-------------------------------------------------------------------------------
#' Optimize a proliferation model using maximum likelihood
#'
#'
opt.prolif.model.mle <- function(x, starts, upper, lower, fixed, peak.stats, plot=F, opt.env=NULL, verbose=F) {

  # Find optimal parameters using L-BFGS-B
  res <- optim(
    par = starts,
    fn = prolif.model.wrapper,
    x = x,
    fixed=fixed,
    n.peaks = nrow(peak.stats),
    upper=upper,
    lower=lower,
    opt.env=opt.env,
    type="neg_log_likelihood",
    method="L-BFGS-B",
    verbose=verbose
  )

  final.par <- c(res$par, fixed)

  # Plot
  if (plot) {

    # Plot parameters over optimization
    if (!is.null(opt.env)) {
      opt.plot.params(x, opt.env)
    }

    ###########
    # Model fit
    hd <- hist(x, breaks=250, plot=F)
    plot(hd$mids,
         hd$density,
         xlab="Log10(intensity)",
         ylab="Density",
         pch=20,
         col="grey", bty="n",
         main="Estimated fit vs optimal fit")

    mod.d <- prolif.model.wrapper(starts, nrow(peak.stats), x = hd$mids, type="density")
    lines(hd$mids,
          mod.d,
          col="blue")

    mod.d <- prolif.model.wrapper(final.par, nrow(peak.stats), x = hd$mids, type="density")
    lines(hd$mids,
          mod.d,
          col="red")

    legend("topleft",
           legend=c("Initial fit", "Optimized fit", "Data"),
           fill=c("blue", "red", "grey"),
           bty="n")
  }


  return(final.par)
}

#-------------------------------------------------------------------------------
#' Optimize a proliferation model using least squares
#'
#'
opt.prolif.model.ls <- function(x, y, starts, upper, lower, fixed, peak.stats, plot=F, opt.env=NULL) {

  # Find optimal model parameters using Levenberg–Marquardt algo
  # This algo tries to find the minimal least squares solution
  res <- nls.lm(
    par = starts,
    fn = prolif.resid,
    y = y,
    x = x,
    fixed=fixed,
    n.peaks = nrow(peak.stats),
    upper=upper,
    lower=lower,
    opt.env=opt.env,
    control = nls.lm.control(
      nprint = F,
      maxiter = 1024,
      factor = 0.01,
      maxfev=1000000,
      ptol=.Machine$double.xmin,
      gtol= 0,
      ftol=.Machine$double.xmin)
  )

  final.par <- c(res$par, fixed)

  if (!is.null(opt.env)) {
    opt.env[["score"]]     <- log2(res$rsstrace)
  }

  #####
  # See how residual sum of squares changes
  if (plot) {

    # Plot parameters over optimization
    if (!is.null(opt.env)) {
      opt.plot.params(x, opt.env)
    }

    # Plot estimated parameters to optimal fit
    plot(x, y,
         main="Estimated fit vs optimal fit",
         xlab="Bin average log10(CTV)",
         ylab="Count",
         pch=20,
         col="grey", bty="n")

    y.pred <- sapply(x, function(x){prolif.model(means=peak.stats$est_mean,
                                                 sd=peak.stats$est_peak_sd,
                                                 summits=peak.stats$est_summit,
                                                 x)})

    lines(y.pred ~ x, lwd=2, col="blue")

    y.pred.optim <- sapply(x, function(x){prolif.model.wrapper(final.par, nrow(peak.stats), x=x, fixed=fixed)})

    lines(y.pred.optim ~ x, lwd=2, col="red")
    legend("topleft",
           legend=c("Initial fit", "Optimized fit", "Smoothed data"),
           fill=c("blue", "red", "grey"),
           bty="n")
  }
  #####

  return(final.par)
}


#-------------------------------------------------------------------------------
#' Fit a proliferation model using least squares
#'
#' Fit peaks on a CTV or CFSE or other tracking dye trace using a binned approach.
#' The trace should be raw per event data from the FCS file on the CTV+ population.
#'
#' @details
#' As a proliferation model, a guassian mixture distribution is fit.
#' Initial parameters (number of peaks, peak means, peak summits) are estimated
#' using a simple approach where first the generation zero peak is estimated based
#' on 'peak.0.lower.bound' on smoothed data. The bins are smoothed using a
#' nearest neighbour smoother in the window provided.
#'
#' Peaks are only called if they exceed peak.thresh.enrich, the relative enrichment
#' between the valley's either side of the peak and the peak summit. Peaks
#' must also be adjacent. So if peak 2 does not pass peak.thresh.enrich but peak
#' 3 does, only 2 peaks are fit. This is done to avoid fitting many peaks in the
#' marginal count range.
#'
#' Using the gen0 peak the subsequent peak position is estimated, after which
#' the summit is found and the peak mean updated to the local summit. This is
#' repeated using each subsequent peak.
#'
#' After estimating starting values, a gaussian mixture distribution is fit to
#' the trace (a sum of the PDF of individuals gaussians). The starting values
#' for the peak standard deviations are estimated from the trace in the bound
#' of peak0. These values are input into the Levenberg-Marquardt algo to
#' optimise with respect to the resdiual sum of squares between the model
#' and the smoothed trace.
#'
#' opt.peak.pos.dev controls the range peaks are allowed to vary in position
#' during optimization. Defaults to +- half the estimated peak distance (NULL).
#' To ignore this, set to Inf. This setting avoids the optimizer putting peaks
#' in the left tail that sometimes might be present.
#'
#' @param x raw FACS intensities
#' @param peak.0.lower.bound the value on log10 scale where the first valley is
#' @param peak.thresh.enrich fold change over valley to call peak in initial estimation (default 1)
#' @param peak.thresh.summit minimum height of a peak in percentage of total heights (default 0.05)
#' @param peak.max the maximum number of peaks to search for (default 12)
#' @param bins number of bins for fitting (default 250)
#' @param smoothing.window number of values up and downstream for the NN smoother (default 2)
#' @param plot should results be plotted
#' @param plot.main title for the plot
#' @param opt.peak.pos.dev limits the range where peak positions can be. See @details
#' @param opt.trim.left.tail Should values that fall outside the model range bet set to zero @details

#' @returns A data frame with peak statistics
#'
#' @examples
#'
#' # Simulate proliferation data
#' y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
#' y <- c(y/4, y/2, y)
#'
#' # Fit peaks
#' peaks <- fit.peaks(y, peak.0.lower.bound=9.8)
#'
#' # Simulate proliferation data #2
#' y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
#' y <- c((y/8)[1:500], y/4, y/2, y[1:500])
#'
#' peaks <- fit.peaks(y, peak.0.lower.bound=9.9, mode="MLE")
#' @export
fit.peaks  <- function(trace,
                      peak.0.lower.bound,
                      bins=250,
                      smoothing.window=2,
                      plot=T,
                      plot.main = "Proliferation model",
                      opt.peak.pos.dev=NULL,
                      opt.trim.left.tail=T,
                      mode="LS",
                      ...) {

  # Plotting
  if (plot) {par(mfrow=c(2,3), mar=c(5,5,5,1))}

  # Construct histogram
  cur.hist <- hist(log10(trace), breaks=bins, plot=F)

  # Nearest neighbor smoother
  y.smth     <- nn.smoother(cur.hist$counts, window=smoothing.window)

  peak.stats <- find.initial.peaks(cur.hist, y.smth, peak.0.lower.bound, plot = plot, ...)

  # If just generation zero is left, do not proceed to optimization
  if (nrow(peak.stats)==1) {
    dummy <- data.frame('opt_mean'=NA, 'opt_summit'=NA, 'opt_peak_sd'=NA, 'peak_area'=NA, 'peak_area_error'=NA, 'peak_area_prop'=NA, 'peak_events'=NA, 'peak_ancestors'=NA)
    return(cbind(peak.stats, dummy))
  }

  #-----------------------------------------------------------------------------
  # Estimate starting parameters
  #-----------------------------------------------------------------------------
  x.mids        <- cur.hist$mids
  y.smth        <- y.smth
  est.gen0.mean <- peak.stats[1, "est_mean"] # Starting generation mean
  est.peak.dist <- est.gen0.mean - peak.stats[2, "est_mean"] # Distance between peaks

  #-----------------------------------------------------------------------------
  # Filter data that cannot be captured by a peak
  # Usually this does not give issues, but if there is a strong CTV- bleed
  # the optimization can reach a wrong optimum by increasing the peak sd too
  # much, so we only want to optimize around the peaks
  if (opt.trim.left.tail) {
    lower.lim <- peak.stats[nrow(peak.stats),"est_mean"] -  (est.peak.dist/2)
    y.smth[x.mids < lower.lim] <- 0
  }

  #-----------------------------------------------------------------------------
  # Roughly estimate the SD of generation zero
  peak0.mean  <- peak.stats[1,"est_mean"]
  h.dist      <- est.peak.dist / 2
  yy          <- y.smth[(x.mids > est.gen0.mean-h.dist) & (x.mids < est.gen0.mean + h.dist)]
  xx          <- x.mids[(x.mids > peak0.mean-h.dist) & (x.mids < peak0.mean+h.dist)]
  est.peak.sd <- sd(xx)
  peak.stats[,"est_peak_sd"] <- est.peak.sd

  #-----------------------------------------------------------------------------
  # Starting parameter list
  starts          <- c(as.list(peak.stats$est_mean), as.list(peak.stats$est_summit))
  mean.names      <- paste0("gen", 0:(nrow(peak.stats)-1), ".mean")
  summit.names    <- paste0("gen", 0:(nrow(peak.stats)-1), ".summit")
  sd.names        <- c("peak.sd")
  names(starts)   <- c(mean.names, summit.names)
  starts$peak.sd  <- est.peak.sd

  # Fix the position of generation 0
  #fixed <- starts[c("gen0.mean")]
  #fixed <- starts[grep("mean", names(starts))]
  fixed  <- NULL
  starts <- starts[!names(starts) %in% names(fixed)]

  #-----------------------------------------------------------------------------
  # Determine the limits of optimization of peak positions
  upper        <- rep(Inf, length(starts))
  names(upper) <- names(starts)
  lower        <- rep(0, length(starts))
  names(lower) <- names(starts)

  # Allow initial peak positions to vary by opt.peak.pos.dev
  # Defaults to the estimated sd of the peak
  # Prevents peaks shifting too much
  if (is.null(opt.peak.pos.dev)) {
    opt.peak.pos.dev <- starts$peak.sd
  }

  par.means        <- grep("mean", names(starts))
  upper[par.means] <- as.numeric(starts[par.means]) + opt.peak.pos.dev
  lower[par.means] <- as.numeric(starts[par.means]) - opt.peak.pos.dev


  # Create environment to store iterations
  opt.env <- opt.new.env(nrow(peak.stats))
  #opt.new.env(nrow(peak.stats))

  if (mode == "LS") {
    #-----------------------------------------------------------------------------
    # Optimize parameters using least squares
    final.par <- opt.prolif.model.ls(x.mids,
                                     y.smth,
                                     starts=starts,
                                     upper=upper,
                                     lower=lower,
                                     fixed=fixed,
                                     peak.stats=peak.stats,
                                     plot=plot,
                                     opt.env=opt.env)


    # total events captured by the model
    #total.events <- sum(y[x > peak.stats[nrow(peak.stats),"opt_mean"] - 3*(peak.stats[nrow(peak.stats),"opt_peak_sd"]) ])

    # Get the model prediction
    y.pred.optim <- prolif.model.wrapper(final.par,
                                         nrow(peak.stats),
                                         x=x.mids)
    total.events <- sum(y.pred.optim)

  } else if (mode == "MLE") {
    #-----------------------------------------------------------------------------
    # Optimize parameters using MLE
    # For MLE we optimize over densities, so scale summits so they are the
    # relative mixture between distributions rather then the number of events
    starts[summit.names] <- as.numeric(starts[summit.names]) / sum(as.numeric(starts[summit.names]))
    upper[summit.names]  <- 1

    final.par <- opt.prolif.model.mle(log10(trace),
                                      starts=starts,
                                      upper=upper,
                                      lower=lower,
                                      fixed=fixed,
                                      peak.stats=peak.stats,
                                      plot=plot,
                                      opt.env=opt.env)

    total.events <- length(trace)
    binwidth <- (x.mids[2] - x.mids[1])

    # Scale the density to so they match counts as the prediction
    y.pred.optim <- prolif.model.wrapper(final.par,
                                         nrow(peak.stats),
                                         x=x.mids,
                                         type="density") #* (total.events / length(x.mids))

    y.pred.optim <- (y.pred.optim*binwidth) * total.events
  } else {
    stop(simpleError("No valid mode provided, must be LS or MLE."))
  }

  mtext(plot.main, side = 3, line = -2, outer = TRUE)
  par(mfrow=c(1,1))

  # Clear environment with optimization results
  rm("opt.env")

  #-----------------------------------------------------------------------------
  if (plot) {

    ####
    plot(x.mids, cur.hist$counts,
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


  # Add the optimized values to the table
  peak.stats$opt_mean    <- as.numeric(final.par[mean.names])
  peak.stats$opt_summit  <- as.numeric(final.par[summit.names])
  peak.stats$opt_peak_sd <- as.numeric(final.par[sd.names])

  #-----------------------------------------------------------------------------
  # Determine the area of each generation
  # In this case as we fix sd above, could also take the ratio of the summits
  for (i in 1:nrow(peak.stats)) {
    int.res <- integrate(prolif.single.peak,
                         lower=min(log10(trace)),
                         upper=max(log10(trace)),
                         mean=as.numeric(final.par[mean.names[i]]),
                         summit=as.numeric(final.par[summit.names[i]]),
                         sd=as.numeric(final.par[sd.names]))

    peak.stats[i,"peak_area"]       <- int.res$value
    peak.stats[i,"peak_area_error"] <- int.res$abs.error
  }

  #-----------------------------------------------------------------------------
  peak.stats$peak_area_prop   <- (peak.stats$peak_area / sum(peak.stats$peak_area)) *100
  peak.stats$peak_summit_prop <- (peak.stats$opt_summit / sum(peak.stats$opt_summit)) *100

  # Determine total number of events in a peak
  peak.stats$peak_events <- total.events * (peak.stats$peak_area_prop/100)

  # Determine number of ancestors in peak
  tmp <- vector()
  for(i in 1:nrow(peak.stats)) {
    tmp[i] <- (peak.stats[i, "peak_events"] / 2^(i-1))
  }

  peak.stats$peak_ancestors <- tmp





  return(peak.stats)
}

#-------------------------------------------------------------------------------
#' Get proliferation statistics
#'
#' Takes output from fit.peaks and calculates proliferation statistics.
#' Proliferation statistics according to the FCS express definitions
#' https://fcsexpressdownloads.s3.amazonaws.com/manual/manual_WIN_RUO/index.html?proliferation_statistics.htm
#'
#' Can be used manually by provided a data frame with two columns
#'
#' - generation
#'
#' - peak_events
#'
#' Where each row represent the peak number and number of events in a peak.
#'
#' @param peak.stats output from fit.peaks (data frame)
#' @returns A dataframe with proliferation statistics
#' @examples
#'
#' # Simulate proliferation data
#' y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
#' y <- c(y/4, y/2, y)
#'
#' # Fit peaks
#' peaks <- fit.peaks.ls(y, 0, 0, peak.0.lower.bound=9.8)
#' get.prolif.stats(peaks)
#'
#' # Manual table
#' peaks <- data.frame(generation=c(0,1,2), peak_events=c(1000, 1000, 1000))
#' get.prolif.stats(peaks)
#' @export
get.prolif.stats <- function(peak.stats) {

  founding.pop <- sum(peak.stats$peak_events / 2^peak.stats$generation)
  divided.pop  <- sum(peak.stats$peak_events[-1] / 2^peak.stats$generation[-1])

  prol.index   <- sum(peak.stats$peak_events) / founding.pop
  div.index    <- sum(peak.stats$peak_events[-1]) / divided.pop
  perc.div     <- (divided.pop / founding.pop)*100

  vec <- c(founding.pop, divided.pop, prol.index, div.index, perc.div, nrow(peak.stats)-1)
  names(vec) <- c("founding_population", "divided_population", "proliferation_index", "division_index", "percent_divided", "num_generations")

  return(vec)
}
