##############################################
# KDE function for single observation sequence #
##############################################

# load packages
require(moments, quietly = TRUE)
require(doParallel, quietly = TRUE)


# get names for the columns of the table characterizing the peaks
col.names.peaks <- function()
{
  v <- c("Number of modes", "Number of modes (5% excluded)",
         "Number of modes (10% excluded)", "Skewness",
         "Mode skewness", "Nonparametric skewness", "Q50 skewness",
         "Absolute Q50 mode skewness", "Absolute Q80 mode skewness",
         "Peak 1", "Probability Mass 1",
         "Peak 2", "Probability Mass 2",
         "Peak 3", "Probability Mass 3",
         "Peak 4", "Probability Mass 4",
         "Peak 5", "Probability Mass 5",
         "Peak 6", "Probability Mass 6",
         "Peak 7", "Probability Mass 7",
         "Peak 8", "Probability Mass 8",
         "Peak 9", "Probability Mass 9",
         "Peak 10", "Probability Mass 10", 
         "Warning close modes",
         "Number close modes", 
         "Modes (close modes excluded)",
         "SD", "IQR 80", "IQR 90",
         "Total number of sequences")
  return(v)
}


# get names for the columns of the table containing the boostrap results
col.names.bs <- function()
{
  v <- c("Number of modes (NM)", 
        "% of samples with same NM",
        "% of samples with more NM",
        "% of samples with less NM",
        "no. of samples with same NM",
        "% BS samples excluded by prob. mass crit.",
        "Warning CI")
  return(v)
}


# plot KDE for one set up CpGo/e ratios
# obs: CpGo/e ratios
# bs.cis: Is bootstrap done?
# t.name: Title of the plot
# t.sub: Text that is added below the title
# t.legend: Is a legend printed?
plot.KDE <- function(obs, t.name, bs.cis = FALSE, bstrp.reps = 1500, conf.lev = 0.95, t.sub = NULL, t.legend = TRUE, min.dist = 0.2, mode.mass = 0.01, band.width = 1.06)
{ 
  # determine directory where functions are located
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  str <- "--file="
  match <- grep(str, cmdArgs)
  if (length(match) == 0) { 
    FCTN.DIR <- "../Galaxy/Functions"
  } else {
    path <- normalizePath( sub(str, "", cmdArgs[match]) )
    FCTN.DIR <- file.path(dirname(path), "Functions")
  }
  
 # part 1: initialize parameters etc
  # ---------------------------------
  
  # table with names and number of peaks
  v <- col.names.peaks()
  tab1.m <- data.frame(matrix(NA, nrow = 1, ncol = length(v)))
  names(tab1.m) <- col.names.peaks()
  # parameters to set in any case
  num.points <- 10 ^ 4 # number of points where to estimate the density
  p.bw <- "nrd" # algorithm for the bandwidth selection, "nrd" for Scott's bandwith
  use.seed <- TRUE
  threshold.modes <- mode.mass 
  threshold.bs.ci <- 0.2 # only changes of +- threshold.bs.ci% in prob. mass is allowed for entering the CI calculation
  # bootstrap with optional parameters
  if (bs.cis) {
    ncpus = max(1, detectCores()) # number of processsors
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)
    # table for the bootstrap
    tab2.m <- data.frame(matrix(NA, nrow = 1, ncol = 7))
    names(tab2.m) <- col.names.bs()
  }


  # part 2: estimate the KD and calculate probability masses
  # --------------------------------------------------------
 
  source( file.path(FCTN.DIR, "density_pm.R") )
  estimate <- density_pm(obs, num.points, p.bw = band.width * sd(obs) / length(obs)^0.2, threshold.modes = threshold.modes)
  ker <- estimate$ker # kernel density
  p <- estimate$peaks
  v <- estimate$valleys
  pm <- estimate$pm
  
  # select for each peak a line type
  num.pm <- length(p) # number of peaks
  p5 <- estimate$p5
  p10 <- estimate$p10
  lty = rep(1, num.pm) #pm > 0.10
  if (length(p10) > 0) { #pm < 0.10
    lty[p10] <- 2
  }
  if (length(p5) > 0) { #pm < 0.05
    lty[p5] <- 4
  }
  p.lty <- lty
  
  # part 3: bootstrap
  # -----------------
  
  if (bs.cis == TRUE) { 
    # do bootstrapping
    estimateB = foreach(j = 1 : bstrp.reps ) %dopar% {
      if (use.seed == TRUE) {set.seed(j)}
      source( file.path(FCTN.DIR, "density_pm.R") )   # As doParallel is used, the source code has to be included here
      obs.boot = obs[sample(seq(obs), replace = TRUE)]
      density_pm_boot(obs.boot, num.points = num.points, p.bw = band.width * sd(obs) / length(obs)^0.2, threshold.modes = threshold.modes)
    }
    stopCluster(cl)
    
    # calculate CIs based on samples where the number of peaks of the KDE coincide with the original data
    # and the probability mass of the peaks of the sample don not divert to strongly from the original data
    # ... extract modes and pm for samples
    num.pmB <- sapply(lapply(estimateB, "[[", "peaks"), length)   # no of peaks before cleaning
    num.peaks.ok <- num.pmB == num.pm    # samples with same number of peaks before cleaning
    pB <- t(sapply( estimateB[num.peaks.ok], "[[", "peaks") )
    pmB <- t(sapply( estimateB[num.peaks.ok], "[[", "pm") )
    if (num.pm == 1) { # put into right matrix form
      pB <- t(pB)
      pmB <- t(pmB)
    }
    
    # ... check if prob. mass changes too strongly for any mode
    t.pmB <- t(pmB)
    keep.ub <- !apply(pm * (1 + threshold.bs.ci) < t.pmB, 2, any)
    keep.lb <- !apply(pm * (1 - threshold.bs.ci) > t.pmB, 2, any)
    mass.ok <- keep.ub & keep.lb
    pB.cl <- as.matrix( pB[mass.ok, ] )
    
    # ... determine CIs
    q <- (1 - conf.lev) / 2
    p.CI <- apply(pB.cl, 2, quantile, probs = c(q, 1 - q))
  }
  
  
  # part 4: plots
  # -------------
    
  #t.breaks <- seq(0, max(obs)*1.05, by = 0.03)
  t.breaks = 50
  hist_data <- hist(obs, breaks = t.breaks, plot = FALSE)
  hist(obs, breaks = t.breaks, prob = TRUE, main = t.name,
       # sub = paste("Gaussian kernel with band width", band.width),
	   sub = "CpGo/e Ratio",
       col = grey(0.9), border = grey(0.6))
  if (!is.null(t.sub)) {
    mtext(t.sub)
  }
  if (bs.cis == TRUE) { # CI
    j <- 1 : ncol(p.CI)
    rect(ker$x[p.CI[1, j]], 0, ker$x[p.CI[2, j]],
         15, density = 20, angle = 45 + (j - 1) * 90, col = "blue") 
  }
  lines(ker, col = "red", lwd = 2) # density
  
  # vertical lines at peaks
  x.pos = ker$x[p]
  ok <- c()
  close <- c()
  for (i in 1:length(x.pos)) {
    b <- TRUE
    if (i > 1) {
      if (x.pos[i] - x.pos[i - 1] < min.dist) {
        b <- FALSE
      }
    } 
    if (i < length(x.pos)) {
      if (x.pos[i + 1] - x.pos[i] < min.dist) {
        b <- FALSE
      }
    } 
     
    if (b) {
      ok <- c(ok, i)
    } else {
      close <- c(close, i)
    }
  }
  # peaks
  if (length(ok) > 0) {
    abline(v = ker$x[p][ok], col = "blue", lwd = 3, lty = p.lty) 
  }
  if (length(close) > 0) {
    abline(v = ker$x[p][close], col = "orange", lwd = 3, lty = p.lty) 
  }
  # valleys
  abline(v = ker$x[v], col = "black", lwd = 1) 
  
  # legend
  if(t.legend == TRUE) {
    t.sym <- expression(""<="", ""<"", "">="")
    thr = threshold.modes
    if (thr >= 0.1) {
      md.labs <- substitute(paste("Mode with PM ", sym3, " ", thr, sep = ""), list(thr = thr, sym3 = t.sym[[3]]))
    } else {
      md.labs <- substitute(paste("Mode with PM ", sym3, " 0.1", sep = ""), list(sym3 = t.sym[[3]]))
      if (thr >= 0.05) {
        md.labs <- c( md.labs, substitute(paste("Mode with ", thr, " ", sym1, " PM ", sym2, " 0.1", sep = ""), list(thr = thr, sym1 = t.sym[[1]], sym2 = t.sym[[2]])) )     
      } else {
        md.labs <- c( md.labs, substitute(paste("Mode with 0.05 ", sym1, " PM ", sym2, " 0.1", sep = ""), list(sym1 = t.sym[[1]], sym2 = t.sym[[2]])),
                               substitute(paste("Mode with ", thr, sym1, " PM ", sym2, " 0.05"), list(thr = thr, sym1 = t.sym[[1]], sym2 = t.sym[[2]])))    
      }
    }
    legend("topright", 
           c(expression("Kernel Density"), md.labs),
           lty = c(1, 1, 2, 4), lwd = c(2, 3, 3, 3),
           col = c("red", "blue", "blue", "blue"), bg = "white") 
  }
  
  
  # part 5: return results
  # ----------------------
  
  # part 5 a): results table 1
  # tbd (maybe): introduce maximum number of modes (10 right now)
  tab1.m[1, "Number of modes"] <- num.pm
  tab1.m[1, 2] = num.pm - length(estimate$p5)
  tab1.m[1, 3] = num.pm - length(estimate$p10)
  for(j in 1 : 10){
    if(num.pm < j){
      tab1.m[1, j * 2 + 8] = c(" ")
      tab1.m[1, j * 2 + 9] = c(" ")
    } else{
      tab1.m[1, j * 2 + 8] = ker$x[p][j]
      tab1.m[1, j * 2 + 9] = pm[j]
    }
  }
  
  # fill table 1 with descriptives
  tab1.m[1, "Skewness"] <- skewness(obs)
  mode <- ker$x[ which.max(ker$y) ]
  tab1.m[1, "Mode skewness"] <- (mean(obs) - mode) / sd(obs)
  tab1.m[1, "Nonparametric skewness"] <- (mean(obs) - median(obs)) / sd(obs)
  q <- quantile(obs, c(0.25, 0.5, 0.75))
  tab1.m[1, "Q50 skewness"] <- (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1])
  tab1.m[1, "Absolute Q50 mode skewness"] <- (q[3] + q[1]) / 2 - mode
  q <- quantile(obs, c(0.1, 0.5, 0.9))
  tab1.m[1, "Absolute Q80 mode skewness"] <- (q[3] + q[1]) / 2 - mode
  tab1.m[1, "SD"] <- sd(obs)
  tab1.m[1, "IQR 80"] <- diff(quantile(obs, c(0.1, 0.9)))
  tab1.m[1, "IQR 90"] <- diff(quantile(obs, c(0.05, 0.95)))
  tab1.m[1, "Total number of sequences"] = length(obs)
  
  # check if any peak is closer than a given threshold to any other
  num.close.modes <- sum(diff(ker$x[p]) < min.dist)
  if ( any(diff(ker$x[p]) < min.dist) && (num.pm > 1) ) {
    tab1.m[1, "Warning close modes"] <- "Modes too close"
    tab1.m[1, "Number close modes"] <- num.close.modes
    tab1.m[1, "Modes (close modes excluded)"] <- num.pm - num.close.modes
  } else {
    tab1.m[1, "Modes (close modes excluded)"] <- num.pm 
  }
  
  # part 5 b): results table 2  
  if (bs.cis == TRUE) {
    ker <- lapply( estimateB, "[[", "ker")
    peaks <- lapply( estimateB, "[[", "peaks")
    num.peaks <- c()
    for (i in 1:length(peaks)) {
      curr.peaks <- ker[[i]]$x[ peaks[[i]] ]
      num.cl <- sum(diff(curr.peaks) < min.dist)  
      num.peaks <- c(num.peaks, length(curr.peaks) - num.cl)
    }
    
    # fill table 2 with stats on number of modes in bs samples 
    num <- num.pm - num.close.modes
    tab2.m[1, "Number of modes (NM)"] <- num
    tab2.m[1, "% of samples with same NM"] <- 100 * sum(num.peaks == num) / bstrp.reps  # equal
    tab2.m[1, "% of samples with more NM"] <- 100 * sum(num.peaks > num) / bstrp.reps  # more
    tab2.m[1, "% of samples with less NM"] <- 100 * sum(num.peaks < num) / bstrp.reps  # less
    if (num.pm > 1) {
      tab2.m[1, "Warning CI"] <- "CI's may be unreliable"
    }
    
    tab2.m[1, "no. of samples with same NM"] <- sum(num.peaks == num)
    tab2.m[1,  "% BS samples excluded by prob. mass crit."] <- (1 - sum(mass.ok) / sum(num.peaks.ok)) * 100
  }

  # return the results
  if (bs.cis){
    return(list(tab.des = tab1.m, tab.bs = tab2.m))
  } else {
    return(list(tab.des = tab1.m))
  }
}



