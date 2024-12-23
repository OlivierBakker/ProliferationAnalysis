#' @import minpack.lm
NULL

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
#' @description
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
#' Find the next peak in a trace
#'
#' @description
#' Find the mode of the next peak using the previous peak and the average peak
#' distance. The window can be asymmetrically scaled using the scaling factors.
#' This is done as the smaller the peaks get, the closer together they tend to be.
find.next.peak <- function(x, y, prev.peak, peak.dist, window.scaling.factors=c(0.25, 0.25)) {

  #est.window      <- peak.dist/2
  est.curpeak     <- log10((10^prev.peak)/2)

  low  <- est.curpeak - (peak.dist * window.scaling.factors[1])
  high <- est.curpeak + (peak.dist * window.scaling.factors[2])

  curpeak.mean    <- find.local.max(x, y, low, high)
  curpeak.summit  <- y[nearest.index(x, curpeak.mean)]

  # calculate enrichment
  fold.change  <- find.enrichment(x, y, curpeak.mean, peak.dist/2)

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
#' @description
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
#'
prolif.model.density <- function(means, sd, summits, x, verbose=F, log=F, opt.env=NULL) {
  density <- 0
  summits <- summits/sum(summits)

  if (verbose) {
    cat("[INFO] means: ",   means, "\n")
    cat("[INFO] sd: ",      sd, "\n")
    cat("[INFO] summits: ", summits, "\n")
    cat("[INFO] --------------------------------------------- \n")
  }

  if (log){
    # https://stats.stackexchange.com/questions/105602/example-of-how-the-log-sum-exp-trick-works-in-naive-bayes
    # Non-vectorized. Works, but inefficient
    #if (!is.null(opt.env)) {
    #  if (!is.null(opt.env[["densities"]])) {
    #    opt.env[["densities"]] <- matrix(nrow=length(means), ncol=length(x))
    #  }
    #  densities <- opt.env[["densities"]]
    #} else {
    #}

    #-------------
    # Implementation from Laplaces Demon
    #density <- dnormm(x, summits, means, sd, log=T)
    #----
    # Manual implementation
    densities <- matrix(nrow=length(means), ncol=length(x))
    for (i in 1:length(means)) {
      cur.d         <- dnorm(x, mean=means[i], sd=sd[i], log=T)
      densities[i,] <- log(summits[i]) + cur.d
    }

    density <- apply(densities, 2, function(cur.dens) {
      a        <- max(cur.dens)
      cur.dens <- exp(cur.dens - a)
      return(a + log(sum(cur.dens)))
    })

    #---------------------
    # Perplexity optimized of above, not much quicker
    # Preallocate memory for matrices
    #densities <- matrix(0, nrow = length(means), ncol = length(x))
    #density <- numeric(length(x))

    # Calculate densities
    #for (i in 1:length(means)) {
    #  cur.d <- dnorm(x, mean = means[i], sd = sd[i], log = TRUE)
    #  densities[i, ] <- log(summits[i]) + cur.d
    #}

    # Calculate density using vectorized operations
    #for (j in 1:length(x)) {
    #  a <- max(densities[, j])
    #  cur_dens <- exp(densities[, j] - a)
    #  density[j] <- a + log(sum(cur_dens))
    #}
    #---------------------


    # # Find max value to scale by, its quicker to re-compute then to store this
    # a         <- c()
    # for (i in 1:length(means)) {
    #   cur.dens    <- log(summits[i]) + dnorm(x, mean=means[i], sd=sd[i], log=T)
    #   if (i==1) {
    #     a         <- cur.dens
    #   } else {
    #     a[a<cur.dens] <- cur.dens[a<cur.dens]
    #   }
    # }
    #
    # # Calculate density in log
    # density   <- c()
    # for (i in 1:length(means)) {
    #   cur.dens    <- log(summits[i]) + dnorm(x, mean=means[i], sd=sd[i], log=T)
    #   if (i==1) {
    #     density   <- exp(cur.dens - a)
    #   } else {
    #     density   <- density + exp(cur.dens - a)
    #   }
    # }
    # density <- a + log(density)

  } else {
    density <- 0
    for (i in 1:length(means)) {
      cur.d   <- dnorm(x, mean=means[i], sd=sd[i])
      density  <- density + (summits[i] * cur.d)
    }

    # Hack to deal with numeric precision limit on very small densities
    if (sum(density==0) >0) {
      msg <- "[DO NOT IGNORE] Zero densities detected. Consider using log=T. This is due to a numeric precision limit in dnorm().
    As a hack setting zeroes to 5e-324. This does invalidate the PDF!!"
      warning(simpleWarning(msg))
      density[density==0] <- 5e-324
    }
  }

  return(density)
}

#-------------------------------------------------------------------------------
#' Negative log likelihood of mixture distribution of Gaussian
#'
prolif.model.nll <- function(means, sd, summits, x, verbose=F, log=T, invert=F, opt.env=NULL) {
  density <- prolif.model.density(means, sd, summits, x, verbose=verbose, log=log, opt.env=opt.env)

  if (log) {
    nll <- -sum(density)
  } else {
    nll <- -sum(log(density))
  }

  if (invert) {
    nll <- -nll
  }

  if (verbose) {
    cat("[INFO] nll: ", nll, "\n")
  }
  return(nll)
}

#-------------------------------------------------------------------------------
#' Wrapper that returns a gaussian proliferation model.
#'
#' @description
#' Takes a vector of paremeters par, number of peaks, and data X
prolif.model.wrapper <- function(par, n.peaks, x, fixed=NULL, type="prolif_model", verbose=F, opt.env=NULL, invert=F, names=NULL, log=T) {

  if (!is.null(names(par))) {
    names      <- names(par)
  }
  par        <- as.numeric(par)
  names(par) <- names

  mean.names   <- paste0("gen", 0:(n.peaks-1), ".mean")
  summit.names <- paste0("gen", 0:(n.peaks-1), ".summit")

  #genX.sd      <- paste0("gen", n.peaks, ".sd")

  if (!is.null(fixed)) {
    par <- c(par, fixed)
  }

  mean.vec   <- as.numeric(par[mean.names])
  summit.vec <- as.numeric(par[summit.names])
  sd.vec     <- as.numeric(par["peak.sd"])

  if (length(sd.vec)==1) {
    sd.vec <- rep(sd.vec, length(mean.vec))
  }

  if ("genX.sd" %in% names(par)){
    if (!is.null(par[["genX.sd"]])) {
      sd.vec[n.peaks] <- par[["genX.sd"]]
    }
  }

  # Log the values to the optimization environment
  if (!is.null(opt.env)) {
    opt.env[["means"]]   <- rbind(opt.env[["means"]], mean.vec)
    opt.env[["sd"]]      <- rbind(opt.env[["sd"]], sd.vec)
    opt.env[["summits"]] <- rbind(opt.env[["summits"]], summit.vec)
  }


  if (verbose) {
    cat("#------------------------------------------------------\n")
    cat("# New block\n")
    cat("means - wrap: ", mean.vec, "\n")
    cat("sd - wrap: ", sd.vec, "\n")
    cat("summits - wrap: ", summit.vec, "\n")
  }

  if (type == "prolif_model") {
    return(prolif.model(mean.vec, sd.vec, summit.vec, x))

  } else if (type == "density") {
    return(prolif.model.density(mean.vec, sd.vec, summit.vec, x, verbose=verbose, log=log, opt.env=opt.env))

  } else if (type == "neg_log_likelihood") {
    score <- prolif.model.nll(mean.vec, sd.vec, summit.vec, x, verbose=verbose, invert=invert, log=log, opt.env=opt.env)
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
#' @param window.scaling.factors Scaling factors for peak differences. See details (default c(0.25, 0.25))
#'
#' @details
#' `window.scaling.factors`
#' window.scaling.factors control the size of the window left and right around the next peak position
#' to find the mode of the next peak. The smaller this value is, the closer to the half intensity the
#' peak estimates will be.
#'
#' I.e. if the estimated distance between peaks is 0.5 and the scaling
#' factors are c(0.25, 0.5) and the current peak is 1.5 the next peak mode will be estimates as:
#' lower <- 0.5 * 0.25
#' upper <- 0.5 * 0.5
#'
#' # Where should the next peak be based on half it's intensity
#' est.peak.pos <- log10(10^1.5/2)
#'
#' Then find the max value in the window est.peak.pos-lower, est.peak.pos+upper
#'
#' @returns data frame with initial peak estimates

find.initial.peaks <- function(x, y, peak.0.lower.bound, peak.thresh.enrich=1, peak.thresh.summit=0.05, peak.max=12, plot=F, window.scaling.factors=c(0.25,0.25)) {
  peak0.mean       <- find.local.max(x, y, peak.0.lower.bound, max(x))

  if (is.na(peak0.mean)) {
    cat("[WARN] No values found in peak.0.lower.bound. Returning NA.\n")
    return(NULL)
  }

  peak0.summit     <- y[nearest.index(x, peak0.mean)]

  # Estimate position of next peak
  est.peak1     <- log10(10^peak0.mean/2)
  est.peak.dist <- peak0.mean - est.peak1
  #est.window    <- est.peak.dist/2
  peak0.enrichment <- find.enrichment(x, y, peak0.mean, est.peak.dist/2)

  prev.peak.mean <- peak0.mean
  peak.stats     <- data.frame()
  lim            <- min(x)

  for (i in 1:peak.max) {

    if (prev.peak.mean-est.peak.dist < lim) {
      next;
    }
    res            <- find.next.peak(x, y, prev.peak.mean, est.peak.dist, window.scaling.factors=window.scaling.factors)
    peak.stats     <- rbind(peak.stats, res)
    prev.peak.mean <- res[1]

  }

  peak.stats[is.na(peak.stats)] <- 0

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

  return(peak.stats)
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
opt.prolif.model.mle <- function(x, starts, upper, lower, fixed, peak.stats, plot=F, opt.env=NULL, verbose=F, log=T) {

  # Find optimal parameters using L-BFGS-B

  #system.time({
  res <- optim(
    par = starts,
    fn = prolif.model.wrapper,
    x = x,
    fixed=fixed,
    n.peaks = nrow(peak.stats),
    upper=upper,
    lower=lower,
    opt.env=opt.env,
    log=log,
    type="neg_log_likelihood",
    method="L-BFGS-B",
    verbose=verbose

  )

  #})
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
          exp(mod.d),
          col="blue")

    mod.d <- prolif.model.wrapper(final.par, nrow(peak.stats), x = hd$mids, type="density")
    lines(hd$mids,
          exp(mod.d),
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
  res <- minpack.lm::nls.lm(
    par = starts,
    fn = prolif.resid,
    y = y,
    x = x,
    fixed=fixed,
    n.peaks = nrow(peak.stats),
    upper=upper,
    lower=lower,
    opt.env=opt.env,
    control = minpack.lm::nls.lm.control(
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
#' Fit a proliferation model using least squares or maximum likelihood
#'
#' Fit peaks on a CTV or CFSE or other tracking dye trace using a binned
#' or maximum likelihood approach. The trace should be raw per event data from
#' the FCS file on the CTV+ population.
#'
#' @param trace raw FACS intensities
#' @param peak.0.lower.bound the value on log10 scale where the first valley is
#' @param peak.thresh.enrich fold change over valley to call peak in initial
#' estimation (default 1)
#' @param peak.thresh.summit minimum height of a peak in percentage of total
#' heights (default 0.05)
#' @param peak.max the maximum number of peaks to search for (default 12)
#' @param bins number of bins for fitting (default 250)
#' @param smoothing.window number of values up and downstream for the NN
#' smoother (default 2)
#' @param window.scaling.factors Controls the window in which initial peaks are found. Vector of length 2
#' @param plot should results be plotted. Overrides plot.optim & plot.final
#' @param plot.optim should optimization plots be plotted
#' @param plot.final should final plot of optimized fit be plotted
#' @param plot.main title for the plot
#' @param opt.peak.pos.dev limits the range in x where peak positions can be
#' optimized within. See @details
#' @param opt.trim.left.tail  Should values that fall outside the initial model
#' range be removed from optimization. See @details
#' @param opt.trim.right.tail Should values that fall outside the initial model
#' range be removed from optimization. See @details
#' @param mode either LS for least squares or 'MLE' for maximum likelihood
#' estimation based optimizer. See @details
#' @param full.out If true output a list with everything needed to easily
#' re-make plots. Useful if you want to regen plots in ggplot for instance.
#' @param peak.thresh.enrich Fold change enrichment of a peak over its valleys
#' @param peak.thresh.summit Relative height of a peak with respect to trace mode
#' @param peak.x.model Should an extra peak be fit that models a trace negative
#' population of cells. See @details
#' @param peak.x.upper.bound If peak.x.model=T, specify the hard limit for
#' finding the mode of peak x. Set to NULL for no limit. (default NULL)
#' @param peak.x.thresh.summit Controls the automated detection of the mode of
#' peak x. Only peaks with a % of the total data mode of this are considered.
#' Should be valued between 0-1. (default 0.05)
#' @param peak.x.thresh.enrich Similar to peak.thresh.enrich but for peak x.
#' (default = 1.1)
#' @param window.scaling.factors Scaling factors for peak finding See details (default c(0.25, 0.25))
#'
#' @details
#' This fits a proliferation model on a trace of raw FACS intensities.
#' The fitting happens in two stages:
#' 1. Initial estimation of number of peaks and their positions
#' 2. Optimization of the model to find final values
#'
#' Initial estimation is done based on thresholds. First the 0 peak
#' is estimated by finding the mode of the data > peak.0.lower.bound.
#'
#' The next peak is then found by finding the mode, at half of the intenstiy
#' of the mode of peak 0. The range where the mode is looked for can be
#' controlled by tweaking `window.scaling.factors`.
#'
#' Peaks are only called if they exceed `peak.thresh.enrich`, the relative
#' enrichment between the valley's either side of the peak and the peak summit
#' and `peak.thresh.summit`, the height of the peak with respect to the mode of
#' the whole trace. So this takes a value between 0-1 describing a percentage of
#' the mode the peak must have. To fit all possible peaks, set both of these to 0
#'
#' Peaks must also be adjacent. So if peak 2 does not pass `peak.thresh.enrich`
#' but peak 3 does, only 2 peaks are fit. This is done to avoid fitting many
#' peaks in the marginal count range.
#'
#' After estimating starting values, a Gaussian mixture distribution is fit to
#' the trace (a sum of the PDF of individuals Gaussian). The starting values
#' for the peak standard deviations are estimated from the trace in the bound
#' of peak0. These values are input into the Levenberg-Marquardt algo to
#' optimize with respect to the residual sum of squares between the model
#' and the smoothed trace in mode 'LS' or using MLE and R optim() when in mode
#' 'MLE'.
#'
#' # `opt.peak.pos.dev`
#' opt.peak.pos.dev controls the range peaks are allowed to vary in position
#' during optimization. Defaults to +- half the estimated peak distance (NULL).
#' To ignore this, set to Inf. This setting avoids the optimizer putting peaks
#' in the left tail that sometimes might be present or separating peaks into
#' nonsensical distances. To contol the allowable range of initial peak
#' estimates see `window.scaling.factors`
#'
#' # `opt.trim.left.tail` / `opt.trim.right.tail`
#' Removes values (MLE) or sets count to zero (LS) of values that fall outside
#' 1sd of the initial model space. This essentially functions as an auto gate
#' after identifying initial peak positions.
#'
#' # `mode`
#' The package offers two optimization schemes:
#'
#' Non linear least squares (LS) with mode=“LS”. Here the data is first binned
#' and smoothed, then parameters are optimized over the smoothed trace using non
#' linear least squares. This method is quick and robust, but relies on binning
#' and some other data processing so is technically less accurate.
#'
#' Maximum likelihood estimation (MLE) based with mode=“MLE”. This mode is a
#' little more proper in that it does not rely on binning or other data tricks
#' to optimize, and takes the full data set to optimize on directly. Initial
#' parameter estimation is still done on a binned smoothed version of the data,
#' but optimization is not.
#'
#' In practice I have found very little difference between these two when the
#' setup is performed correctly.
#'
#' # `window.scaling.factors`
#' window.scaling.factors control the size of the window left and right around the next peak position
#' to find the mode of the next peak. The smaller this value is, the closer to the half intensity the
#' peak estimates will be.
#'
#' I.e. if the estimated distance between peaks is 0.5 and the scaling
#' factors are c(0.25, 0.5) and the current peak is 1.5 the next peak mode will be estimates as:
#' lower <- 0.5 * 0.25
#' upper <- 0.5 * 0.5
#'
#' # Where should the next peak be based on half it's intensity
#' est.peak.pos <- log10(10^1.5/2)
#'
#' Then find the max value in the window est.peak.pos-lower, est.peak.pos+upper
#'
#' @returns A data frame with peak statistics or list of objects if full.out=T
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
#'
#' @export
fit.peaks  <- function(trace,
                      peak.0.lower.bound,
                      bins=250,
                      smoothing.window=2,
                      plot.optim=T,
                      plot.final=T,
                      plot=T,
                      plot.main = "Proliferation model",
                      opt.peak.pos.dev=NULL,
                      window.scaling.factors=c(0.25,0.25),
                      opt.trim.left.tail=F,
                      opt.trim.right.tail=F,
                      mode="LS",
                      full.out=F,
                      peak.x.model=F,
                      peak.x.upper.bound=NULL,
                      peak.x.position=NULL,
                      peak.x.fixed=F,
                      peak.x.thresh.summit=0.05,
                      peak.x.thresh.enrich=1.2,
                      verbose=F,
                      log=T,
                      ...) {

  # Plotting
  if (!plot) {
    plot.optim <- F
    plot.final <- F
  }
  if (plot.optim) {par(mfrow=c(2,3), mar=c(5,5,5,1))}

  # Construct histogram
  cur.hist <- hist(log10(trace), breaks=bins, plot=F)

  # Nearest neighbor smoother
  y.smth     <- nn.smoother(cur.hist$counts, window=smoothing.window)
  y.raw      <- log10(trace)

  # Find the initial peak estimates
  peak.stats <- find.initial.peaks(cur.hist$mids, y.smth, peak.0.lower.bound, plot = plot.optim, window.scaling.factors=window.scaling.factors,...)
  #peak.stats <- find.initial.peaks(cur.hist$mids, y.smth, peak.0.lower.bound, plot = plot.optim)

  #-----------------------------------------------------------------------------
  # Roughly estimate peak X and update stats table
  if (peak.x.model) {
    if (!peak.x.fixed) {

      if (!is.null(peak.x.position)) {
        peak.x.mode <- peak.x.position
      } else {
        peak.x.mode <- find.peak.x.approx.mode(y.raw,
                                               peak.x.thresh.summit=peak.x.thresh.summit,
                                               peak.x.thresh.enrich=peak.x.thresh.enrich,
                                               peak.x.upper.bound=peak.x.upper.bound)
      }

      if (!is.null(peak.x.mode)) {

        #if (peak.x.mode < peak.stats[nrow(peak.stats), "est_mean"]) {
          est.smt       <- y.smth[nearest.index(cur.hist$mids, peak.x.mode)]
          est.peak.dist <- peak.stats[1,1] - peak.stats[2,1]
          est.enrich    <- find.enrichment(cur.hist$mids, y.smth, peak.x.mode, est.peak.dist/2)

          # Only keep peaks with estimated peak sd from the mode
          #peak.stats    <- peak.stats[peak.stats$est_mean > (peak.x.mode + est.peak.dist),]

          # Only keep peaks with are > mode of peak.x
          peak.stats    <- peak.stats[peak.stats$est_mean > peak.x.mode,]

          peak.stats    <- rbind(peak.stats,
                                 c(peak.x.mode, est.smt, est.smt, nrow(peak.stats)+1, nrow(peak.stats), 0))

          peak.stats$est_summit_percentage <- peak.stats$est_summit / sum(peak.stats$est_summit)
        #}

      } else {
        peak.x.mode = NULL
        peak.x.model = FALSE
        warning("No valid mode for peak.x found, skipping fitting. Adjust parameters if it should be there.")
      }

    } else {
      peak.x.mode = peak.x.position
      peak.x.model = TRUE
    }
  } else {
    peak.x.mode = NULL
    peak.x.model = FALSE
  }

  # Plot
  if (plot.optim) {
   opt.plot.estimates(cur.hist, y.smth, peak.stats, peak.0.lower.bound,peak.x.model, peak.x.mode, peak.x.upper.bound)
  }

  # If just generation zero is left, do not proceed to optimization
  #if (nrow(peak.stats)==1) {
  #  dummy <- data.frame('opt_mean'=NA, 'opt_summit'=NA, 'opt_peak_sd'=NA, 'peak_area'=NA, 'peak_area_error'=NA, 'peak_area_prop'=NA, 'peak_events'=NA, 'peak_ancestors'=NA)
  #  return(cbind(peak.stats, dummy))
  #}

  #-----------------------------------------------------------------------------
  # Estimate starting parameters
  #-----------------------------------------------------------------------------
  x.mids        <- cur.hist$mids
  est.gen0.mean <- peak.stats[1, "est_mean"] # Starting generation mean
  est.gen1.mean <- peak.stats[2, "est_mean"]

  if (is.na(est.gen1.mean)) {
    est.gen1.mean <- log10((10^est.gen0.mean)/2)
  }

  est.peak.dist <- est.gen0.mean - est.gen1.mean # Distance between peaks

  #-----------------------------------------------------------------------------
  # Filter data that cannot be captured by a peak
  # Usually this does not give issues, but if there is a strong CTV- bleed
  # the optimization can reach a wrong optimum by increasing the peak sd too
  # much, so we only want to optimize around the peaks
  if (opt.trim.left.tail) {
    lower.lim                  <- peak.stats[nrow(peak.stats),"est_mean"] -  (est.peak.dist/2)
    y.smth[x.mids < lower.lim] <- 0
    y.raw                      <- y.raw[y.raw > lower.lim]
  }

  if (opt.trim.right.tail) {
    upper.lim                  <- peak.stats[1,"est_mean"] +  (est.peak.dist/2)
    y.smth[x.mids > upper.lim] <- 0
    y.raw                      <- y.raw[y.raw < upper.lim]
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
  names(starts)   <- c(mean.names, summit.names)
  starts$peak.sd  <- est.peak.sd

  if (peak.x.model) {
    sd.names <- c(rep("peak.sd", nrow(peak.stats)-1), "genX.sd")
  } else {
    sd.names <- rep("peak.sd", nrow(peak.stats))
  }

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
    opt.peak.pos.dev <- est.peak.dist/2 #starts$peak.sd/2
  }

  par.means        <- grep("mean", names(starts))
  upper[par.means] <- as.numeric(starts[par.means]) + (opt.peak.pos.dev * window.scaling.factors[2])
  lower[par.means] <- as.numeric(starts[par.means]) - (opt.peak.pos.dev * window.scaling.factors[1])

  # Set starting lower limits on sd (can never be zero)
  lower["peak.sd"] <-1e-16

  # Create environment to store iterations
  opt.env <- opt.new.env(nrow(peak.stats))

  # Add the constraints and start for genX
  if (peak.x.model) {
    starts[["genX.sd"]] <- starts$peak.sd*2
    lower[["genX.sd"]]  <- lower["peak.sd"]
    upper[["genX.sd"]]  <- upper["peak.sd"]

    if (!is.null(peak.x.upper.bound)) {
      upper[[par.means[length(par.means)]]] <- peak.x.upper.bound
    } else {
      upper[[par.means[length(par.means)]]] <- upper[[par.means[length(par.means)]]]+starts$peak.sd*2

    }

    if (peak.x.fixed) {
      if (!is.null(peak.x.position)) {
        if (verbose) {cat("[INFO] Fixing peak x position")}
        starts[[par.means[length(par.means)]]] <- peak.x.position
        upper[[par.means[length(par.means)]]] <- peak.x.position
        lower[[par.means[length(par.means)]]] <- peak.x.position
      } else {
        stop("peak.x.fixed=T but no position provided")
      }
    }

  }

  #-----------------------------------------------------------------------------

  if (verbose) {

    cat("[INFO] starting params: ", "\n")
    print(starts)
    cat("[INFO] starting upper:  ", "\n")
    print(upper)
    cat("[INFO] starting lower:  ", "\n")
    print(lower)
    cat("[INFO] number of bins:  ", length(x.mids), "\n")

  }

  #-----------------------------------------------------------------------------
  # Run optimization
  if (mode == "LS") {
    # Optimize parameters using least squares
    final.par <- opt.prolif.model.ls(x.mids,
                                     y.smth,
                                     starts=starts,
                                     upper=upper,
                                     lower=lower,
                                     fixed=fixed,
                                     peak.stats=peak.stats,
                                     plot=plot.optim,
                                     opt.env=opt.env)

    # total events captured by the model
    #total.events <- sum(y[x > peak.stats[nrow(peak.stats),"opt_mean"] - 3*(peak.stats[nrow(peak.stats),"opt_peak_sd"]) ])

    # Get the model prediction
    y.pred.optim <- prolif.model.wrapper(final.par,
                                         nrow(peak.stats),
                                         x=x.mids)

    # Total number of events being modeled
    total.events <- sum(y.smth)

    # Add the optimized values to the table
    peak.stats$opt_summit        <- as.numeric(final.par[summit.names])

    # Determine the area under the peak of each generation
    for (i in 1:nrow(peak.stats)) {
      int.res <- integrate(prolif.single.peak,
                           lower=min(y.raw),
                           upper=max(y.raw),
                           mean=as.numeric(final.par[mean.names[i]]),
                           summit=as.numeric(final.par[summit.names[i]]),
                           sd=as.numeric(final.par[sd.names[i]]))

      peak.stats[i,"peak_area"]       <- int.res$value
      peak.stats[i,"peak_area_error"] <- int.res$abs.error
    }

  } else if (mode == "MLE") {
    # Optimize parameters using MLE
    # For MLE we optimize over densities, so scale summits so they are the
    # relative mixture between distributions rather then the number of events
    # NOTE: this interacts with the SD, to the summits are actually proportional
    # to the area (probability) under the individual peak as opposed to the LS
    # model where counts are modeled explicitly and we need to integrate
    # to find the area.
    starts[summit.names] <- as.numeric(starts[summit.names]) / sum(as.numeric(starts[summit.names]))
    upper[summit.names]  <- 1

    final.par <- opt.prolif.model.mle(y.raw,
                                      starts=starts,
                                      upper=upper,
                                      lower=lower,
                                      fixed=fixed,
                                      peak.stats=peak.stats,
                                      plot=plot.optim,
                                      opt.env=opt.env,
                                      verbose=verbose,
                                      log=log)

    total.events <- length(y.raw)
    binwidth     <- (x.mids[2] - x.mids[1])

    # Scale the density to so they match counts as the prediction
    y.pred.optim <- exp(prolif.model.wrapper(final.par,
                                         nrow(peak.stats),
                                         x=x.mids,
                                         type="density"))

    y.pred.optim <- (y.pred.optim*binwidth) * total.events

    # Predict the values in the histogram of the summit positions
    y.pred.summit <- exp(prolif.model.wrapper(final.par,
                                nrow(peak.stats),
                                x=as.numeric(final.par[mean.names]),
                                type="density"))
    y.pred.summit <- (y.pred.summit*binwidth) * total.events

    peak.stats$opt_summit        <- y.pred.summit
    peak.stats$peak_area         <- final.par[summit.names]
    peak.stats$peak_area_error   <- 0
  } else {
    stop(simpleError("No valid mode provided, must be LS, MIXDIST or MLE."))
  }

  # Add the optimized values to the table
  peak.stats$opt_mean          <- as.numeric(final.par[mean.names])
  peak.stats$opt_peak_sd       <- as.numeric(final.par[sd.names])

  # Add the plot title
  if (plot.optim) {
    mtext(plot.main, side = 3, line = -2, outer = TRUE)
    par(mfrow=c(1,1))
  }

  # Clear environment with optimization results
  rm("opt.env")

  #-----------------------------------------------------------------------------
  # Determine total number of events in a peak
  peak.stats$peak_area_prop   <- (peak.stats$peak_area / sum(peak.stats$peak_area)) *100
  peak.stats$peak_summit_prop <- (peak.stats$opt_summit / sum(peak.stats$opt_summit)) *100
  peak.stats$peak_events      <- total.events * (peak.stats$peak_area_prop/100)

  # Determine number of ancestors in peak
  tmp <- vector()
  for(i in 1:nrow(peak.stats)) {
    tmp[i] <- (peak.stats[i, "peak_events"] / 2^(i-1))
  }
  peak.stats$peak_ancestors <- tmp


  if (peak.x.model) {
    peak.stats[nrow(peak.stats), "generation"] <- "x"
  }

  #-----------------------------------------------------------------------------
  # Should a separate plot be made with he final result on the binned data
  if (plot.final) {
    #opt.plot.final(x.mids, cur.hist$counts, y.pred.optim, main=plot.main)
    opt.plot.final.pp(x.mids, cur.hist$counts, peak.stats, main=plot.main)
  }

  if (full.out) {
    return(list(peak.stats=peak.stats,
                hist=cur.hist,
                x=x.mids,
                y=y.pred.optim,
                par=final.par))
  } else {
    return(peak.stats)
  }
}

#-------------------------------------------------------------------------------
#' Find the mode of the first left sided peak in the data, given some
#' constraints.
#'
#' @export
find.peak.x.approx.mode <- function(trace, peak.x.thresh.summit, peak.x.thresh.enrich, peak.x.upper.bound=NULL, ...) {

  if (is.null(peak.x.upper.bound)) {
    peak.x.upper.bound <- max(trace)
  }

  # Calculate the kernel density estimate of the trace
  dens <- density(trace, ...)

  # Scale by the mode of the trace, so the mode is 1
  d     <- dens$y / max(dens$y)

  # Filter scaled density estimate on peak height
  dd    <- d[d > peak.x.thresh.summit]

  # If there are too few values passing peak.x.thresh.summit return
  if (length(dd) < 3) {
    return(NULL)
  }

  # Keep track of the x values matching dd
  xx    <- dens$x[d >peak.x.thresh.summit]

  # Find the max value and exit after some level of enrichment has been reached
  max.d <- 0
  for (i in 1:length(dd)) {
    if (dd[i] > max.d) {
      max.d <- dd[i]
    }
    if ((max.d/dd[i]) >= peak.x.thresh.enrich) {
      break
    }
  }

  # Which x value matches the estimated mode
  max.x <- xx[which(dd == max.d)]

  # If the max value falls outside the hard limit return the x where the largest
  # density is found within that limit
  if (max.x > peak.x.upper.bound) {
    msg <- "No peak found at specified thresholds. Returning the mode witin hard limit set by peak.x.upper.bound"
    warning(simpleWarning(msg))

    return(NULL)
    # TODO: Check if this can be overidden
    tmp <- dd[xx < peak.x.upper.bound]

    if (length(tmp) > 2) {
      return(xx[which(dd ==  max(tmp))])
    } else {
      return(NULL)
    }

  } else {
    return(max.x)
  }
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

  total.events <- sum(peak.stats[,"peak_events"])
  gen.x.events <- peak.stats[peak.stats$generation=="x","peak_events"]

  if (length(gen.x.events) == 0) {
    gen.x.events <- NA
  }
  #gen.x <- peak.stats[peak.stats$generation != "x",]

  # Filter peak X
  peak.stats <- peak.stats[peak.stats$generation != "x",]
  peak.stats$generation <- as.numeric(peak.stats$generation)

  founding.pop <- sum(peak.stats$peak_events / 2^peak.stats$generation)
  divided.pop  <- sum(peak.stats$peak_events[-1] / 2^peak.stats$generation[-1])

  prol.index   <- sum(peak.stats$peak_events) / founding.pop
  div.index    <- sum(peak.stats$peak_events[-1]) / divided.pop
  perc.div     <- (divided.pop / founding.pop)*100

  vec <- c(founding.pop, divided.pop, prol.index, div.index, perc.div, nrow(peak.stats)-1, total.events, gen.x.events)
  names(vec) <- c("founding_population", "divided_population", "proliferation_index", "division_index", "percent_divided", "num_generations", "total_events", "events_peak_x")

  return(vec)
}
