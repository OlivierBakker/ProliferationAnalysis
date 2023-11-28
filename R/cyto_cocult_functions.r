#-------------------------------------------------------------------------------
# Find the maximum value of y within a range of x values
find.local.max <- function(x, y, lower.x, upper.x){
  y.f <- y[x >= lower.x & x <= upper.x]
  x.f <- x[x >= lower.x & x <= upper.x]
  ymax <- max(y.f)
  return(mean(x.f[y.f==ymax]))
}

#-------------------------------------------------------------------------------
# Find the minimum value y of x within a range of x values
find.local.min <- function(x, y, lower.x, upper.x){
  y.f <- y[x >= lower.x & x <= upper.x]
  x.f <- x[x >= lower.x & x <= upper.x]
  ymin <- min(y.f)
  return(mean(x.f[y.f==ymin]))
}

#-------------------------------------------------------------------------------
# Smooth a trace using a kernel of size window
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
# Find the index of the value in x that is closest to value
nearest.index <- function(x, value) {
  which.min(abs(x-value))
}

#-------------------------------------------------------------------------------
# Find the mode of the next peak using the previous peak and the average peak
# distance
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
# Find the relative enrichment over a valley 
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
# Gaussian of a single peak
gaussian.prolif.single.peak <- function(x, mean, sd, summit) {
  # Gaussian density
  cur.density <- dnorm(x, mean=mean, sd=sd)
  
  # Density so the density on the mode / mean is equal to one
  cur.density <- cur.density / dnorm(mean, mean=mean, sd=sd)
  
  # Scale by the summit
  return(cur.density*summit)
}

#-------------------------------------------------------------------------------
# Gaussian prolif model, it is a mixture distribution scaled by relative
# peak heights.
gaussian.prolif.model <- function(means, sd, summits, x) {
  
  if (length(sd)==1) {
    sd <- rep(sd, length(means))
  }
  
  y.pred <- 0
  
  for (i in 1:length(means)) {
    cur.pred <- gaussian.prolif.single.peak(x, means[i], sd[i], summits[i])
    y.pred   <- y.pred + cur.pred
  }
  
  return(y.pred)
}

#-------------------------------------------------------------------------------
# Wrapper that returns a gaussian proliferation model.
# Takes a vector of paremeters par, number of peaks, and data X
gaussian.prolif.model.wrapper <- function(par, n.peaks, x, fixed=NULL) {
  
  mean.names   <- paste0("gen", 0:(n.peaks-1), ".mean")
  summit.names <- paste0("gen", 0:(n.peaks-1), ".summit")
  
  if (!is.null(fixed)) {
    par <- c(par, fixed)
  }
  
  mean.vec   <- as.numeric(par[mean.names])
  summit.vec <- as.numeric(par[summit.names])
  sd         <- par$peak.sd
  
  return(gaussian.prolif.model(mean.vec, sd, summit.vec, x))
} 

#-------------------------------------------------------------------------------
# Residuals between observed CTV trace and gaussian prolif model
gaussian.prolif.resid <- function(par, n.peaks, x, y, fixed=NULL) {
  return(y - gaussian.prolif.model.wrapper(par, n.peaks, x, fixed))
}

#-------------------------------------------------------------------------------
# Fit peaks on a CTV or CFSE or other tracking dye trace using a binned approach.
# The trace should be raw per event data from the FCS file on the CTV+ population.

# As a proliferation model, a guassian mixture distribution is fit. 
# Initial parameters (number of peaks, peak means, peak summits) are estimated
# using a simple approach where first the generation zero peak is estimated based
# on 'peak.0.lower.bound' on smoothed data. The bins are smoothed using a
# nearest neighbour smoother in the window provided.

# Peaks are only called if they exceed peak.threshold, the relative enrichment
# between the valley's either side of the peak and the peak summit. Peaks
# must also be adjacent. So if peak 2 does not pass peak.threshold but peak
# 3 does, only 2 peaks are fit. This is done to avoid fitting many peaks in the
# marginal count range.

# Using the gen0 peak the subsequent peak position is estimated, after which
# the summit is found and the peak mean updated to the local summit. This is
# repeated using each subsequent peak.

# After estimating starting values, a gaussian mixture distribution is fit to
# the trace (a sum of the PDF of individuals gaussians). The starting values 
# for the peak standard deviations are estimated from the trace in the bound 
# of peak0. These values are input into the Levenberg-Marquardt algo to
# optimise with respect to the resdiual sum of squares between the model
# and the smoothed trace.

# opt.peak.pos.dev sets the maximum deviation allowed after the last peak
# during optimization. Defaults to 1sd of the estimated peak width.
# This avoids issues with gating, where a tail might be present stemming
# from bleedthrough of a CTV- population. To ignore this, set to Inf. 
fit.peaks <- function(cur.trace, bins, smoothing.window, peak.thresh.enrich,
                      peak.thresh.summit, peak.0.lower.bound, max.peaks=12,
                      plot=T,
                      plot.main = "Proliferation model",
                      opt.peak.pos.dev=NULL) {
  
  # Construct histogram
  cur.hist <- hist(log10(cur.trace), breaks=bins, plot=F)
  
  # Kolmogorov–Zurbenko filter to smooth the counts
  #y <- kz(cur.hist$counts, m=smoothing.window)
  
  # Nearest neighbour smoother
  y <- nn.smoother(cur.hist$counts, window=smoothing.window)
  
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
  
  for (i in 1:max.peaks) {
    
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
  
  # Always force starting generation
  if (n.peaks==0) {
    n.peaks <- 1
  }
  
  peak.stats <- peak.stats[1:n.peaks,, drop=F]
  
  #####
  # Plot
  if (plot) {
    par(mfrow=c(2,2), mar=c(5,5,5,1))
    
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
    
    mtext(plot.main, side = 3, line = -2, outer = TRUE)
    
    
    legend("topleft",
           legend=c("Smoothed data", "Gen0"),
           fill=c("blue", "red"),
           bty="n")
    
  }
  #####
  
  # If just generation zero is left, do not proceed to optimization
  if (n.peaks==1) {
    dummy <- data.frame('opt_mean'=NA, 'opt_summit'=NA, 'opt_peak_sd'=NA, 'peak_area'=NA, 'peak_area_error'=NA, 'peak_area_prop'=NA, 'peak_events'=NA, 'peak_ancestors'=NA)
    return(cbind(peak.stats, dummy))
  }
  
  #-------------------------------------------------------------------------------
  
  # Starting parameters
  x             <- cur.hist$mids
  y             <- y
  est.gen0.mean <- peak.stats[1, "est_mean"] # Starting generation mean
  est.peak.dist <- est.gen0.mean - peak.stats[2, "est_mean"] # Distance between peaks
  
  #-------------------------------------------------------------------------------
  # Filter data that cannot be captured by a peak
  # Usually this does not give issues, but if there is a strong CTV- bleed
  # the optimization can reach a wrong optimum by increasing the peak sd too
  # much, so we only want to optimise arround the peaks
  
  lower.lim <- peak.stats[nrow(peak.stats),"est_mean"] -  (est.peak.dist/2)
  #y <- y[x > lower.lim]
  #x <- x[x > lower.lim]
  y[x < lower.lim] <- 0
  #-------------------------------------------------------------------------------
  # Estimate the SD of peak zero
  h.dist <- est.peak.dist / 2
  yy <- y[(x > est.gen0.mean-h.dist) & (x < est.gen0.mean + h.dist)]
  xx <- x[(x > peak0.mean-h.dist) & (x < peak0.mean+h.dist)]
  est.peak.sd <-  sd(xx)
  
  # Starting parameter list
  starts          <- c(as.list(peak.stats$est_mean), as.list(peak.stats$est_summit))
  mean.names      <- paste0("gen", 0:(n.peaks-1), ".mean")
  summit.names    <- paste0("gen", 0:(n.peaks-1), ".summit")
  names(starts)   <- c(mean.names, summit.names)
  starts$peak.sd  <- est.peak.sd
  
  # Fix the position of peak 0
  #fixed <- starts[c("gen0.mean")]
  #fixed <- starts[grep("mean", names(starts))]
  fixed  <- NULL
  starts <- starts[!names(starts) %in% names(fixed)]
  
  # Deterimine the limits of optimisation of peak positions 
  upper        <- rep(Inf, length(starts))
  names(upper) <- names(starts)
  lower        <- rep(-Inf, length(starts))
  names(lower) <- names(starts)
  
  # Allow initial peak positions to vary by opt.peak.pos.dev
  # Defaults to the estimated sd of the peak
  # Prevents peaks shifting too much
  if (is.null(opt.peak.pos.dev)) {
    opt.peak.pos.dev <- starts$peak.sd
  }
  
  par.means <- grep("mean", names(starts))
  upper[par.means] <- as.numeric(starts[par.means]) + opt.peak.pos.dev
  lower[par.means] <- as.numeric(starts[par.means]) - opt.peak.pos.dev
  
  # Find optimal model parameters using Levenberg–Marquardt algo
  # This algo tries to find the minimal least squares solution
  res <- nls.lm(
    par = starts,
    fn = gaussian.prolif.resid,
    y = y,
    x = x,
    fixed=fixed,
    n.peaks = n.peaks,
    upper=upper,
    lower=lower,
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
  
  #####
  # See how residual sum of squares changes
  if (plot) {
    plot(log2(res$rsstrace), 
         xlab="Iteration",
         ylab="log2 - Residual sum of squares",
         main="Optimization result",
         bty="n",
         pch=20)
    lines(log2(res$rsstrace))
    
    # Plot estimated parameters to optimal fit
    plot(x, y,
         main="Estimated fit vs optimal fit",
         xlab="Bin average log10(CTV)",
         ylab="Count",
         pch=20,
         col="grey", bty="n")
    
    y.pred <- sapply(x, function(x){gaussian.prolif.model(means=peak.stats$est_mean,
                                                          sd=est.peak.sd,
                                                          summits=peak.stats$est_summit, x)})
    
    lines(y.pred ~ x, lwd=2, col="blue")
    
    y.pred.optim <- sapply(x, function(x){gaussian.prolif.model(means=as.numeric(final.par[mean.names]),
                                                                sd=final.par$peak.sd,
                                                                summits=as.numeric(final.par[summit.names]), x)})
    
    lines(y.pred.optim ~ x, lwd=2, col="red")
    legend("topleft",
           legend=c("Initial fit", "Optimized fit", "Smoothed data"),
           fill=c("blue", "red", "grey"),
           bty="n")
    
    
    ####
    plot(x, cur.hist$counts,
         main="Optimal fit vs data",
         xlab="Bin average log10(CTV)",
         ylab="Count",
         pch=20,
         col="grey", bty="n")
    lines(y.pred.optim ~ x, lwd=2, col="red")
    
    legend("topleft",
           legend=c("Data", "Optimized fit"),
           fill=c("grey", "red"),
           bty="n")
    
  }
  #####
  
  # Add the optimized values to the table
  peak.stats$opt_mean    <- as.numeric(final.par[mean.names])
  peak.stats$opt_summit  <- as.numeric(final.par[summit.names])
  peak.stats$opt_peak_sd <- as.numeric(final.par$peak.sd)
  
  # Determine the area of each generation
  for (i in 1:n.peaks) {
    int.res <- integrate(gaussian.prolif.single.peak,
                         lower=min(x),
                         upper=max(x),
                         mean=as.numeric(final.par[mean.names[i]]),
                         summit=as.numeric(final.par[summit.names[i]]),
                         sd=as.numeric(final.par$peak.sd))
    
    peak.stats[i,"peak_area"]       <- int.res$value
    peak.stats[i,"peak_area_error"] <- int.res$abs.error
  }
  
  peak.stats$peak_area_prop <- (peak.stats$peak_area / sum(peak.stats$peak_area)) *100
  
  # Events within the peak range
  #total.events <- sum(y[x > peak.stats[nrow(peak.stats),"opt_mean"] - 3*(peak.stats[nrow(peak.stats),"opt_peak_sd"]) ])
  total.events <- sum(gaussian.prolif.model(peak.stats$opt_mean, peak.stats$opt_peak_sd, peak.stats$opt_summit, x))
  
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
# Takes output from fit.peaks and calculates proliferation statistics
# FCS express definitions
# https://fcsexpressdownloads.s3.amazonaws.com/manual/manual_WIN_RUO/index.html?proliferation_statistics.htm
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
