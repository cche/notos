# Carry out analysis of CpGo/E data for Galaxy module
# Ingo Bulla
# 27 Jan 16

# load packages
pckg <- c("methods", "optparse")
for (p in pckg) {
  if (!(p %in% rownames(installed.packages()))) {
    stop( paste("R package", p , "is not installed"), call. = FALSE)
  }
}
require(methods, quietly = TRUE)
require(optparse, quietly = TRUE)

# determine directory where functions are located
cmdArgs <- commandArgs(trailingOnly = FALSE)
str <- "--file="
match <- grep(str, cmdArgs)
if (length(match) == 0) {
  stop("notos.r not set up to be called from R console")
}
path <- normalizePath( sub(str, "", cmdArgs[match]) )
FCTN.DIR <- file.path(dirname(path), "Functions")

source( file.path( FCTN.DIR, "Kernel_function_form.R") )


MAX.CPGOE <- 10   # maximum value for CpGo/e ratios


# process outliers and return quantities characterizing the distribution
# obs: CpGo/e ratios
proc.outliers <- function(obs, frac.outl) {
  ret <- list()

  # remove all zeros from sample
  no.obs.raw <- length(obs)
  ret[["prop.zero"]] <- sum(obs == 0) / no.obs.raw
  obs <- obs[obs != 0]
  if (length(obs) < 3) {
    ret[["valid"]] <- FALSE
    return(ret)
  }

  # replace very large values by a maximum value
  obs <- sapply(obs, function(x) min(x, MAX.CPGOE))

  # defining variables
  # ... mean, median and standard deviation
  ret[["mu.obs"]] <- mu.obs <- mean(obs)
  ret[["me.obs"]] <- me.obs <- median(obs)
  sd.obs <- sd(obs)
  iqr.obs <- IQR(obs)

  # ... uppper and lower limits, based on mean +- k * sd, med. +- k * iqr, k = 2, ..., 4
  ul.mu <- mu.obs + (2 : 5) * sd.obs
  ll.mu <- mu.obs - (2 : 5) * sd.obs
  ul.me <- quantile(obs, 0.75) + (2 : 5) * iqr.obs
  ll.me <- quantile(obs, 0.25) - (2 : 5) * iqr.obs
  names(ul.mu) <- names(ll.mu) <- 2 : 5
  names(ul.me) <- names(ll.me) <- 2 : 5
  ret[["ul.mu"]] <- ul.mu
  ret[["ll.mu"]] <- ll.mu
  ret[["ul.me"]] <- ul.me
  ret[["ll.me"]] <- ll.me

  # summary statistics and data output
  # ... calculate proportion of data excluded when using different ranges
  ret[["prop2"]] <- prop2 <- length(obs[obs < ll.me["2"] | ul.me["2"] < obs]) / no.obs.raw
  ret[["prop3"]] <- prop3 <- length(obs[obs < ll.me["3"] | ul.me["3"] < obs]) / no.obs.raw
  ret[["prop4"]] <- prop4 <- length(obs[obs < ll.me["4"] | ul.me["4"] < obs]) / no.obs.raw
  ret[["prop5"]] <- prop5 <- length(obs[obs < ll.me["5"] | ul.me["5"] < obs]) / no.obs.raw
  # ... choose k in Q1 / Q3 +- k * IQR such that no more than 1% of the data are excluded
  v <- c(prop2, prop3, prop4, prop5) < frac.outl

  if (any(v)) {
    excl.crit <- min(which(v))
	ret[["obs.nz"]] <- obs
    ret[["obs.cl"]] <- obs[!(obs < ll.me[excl.crit] | ul.me[excl.crit] < obs)]
    ret[["used"]] <- paste(2 : 5, "iqr", sep = "")[excl.crit]
  } else {
	ret[["obs.nz"]] <- obs
    ret[["obs.cl"]] <- obs
    ret[["used"]] <- "too few values"
  }
  ret[["valid"]] <- TRUE
  return(ret)
}


# Read CpGo/e ratios from file
# warn: issue warning if necessary
read.CpGoe <- function(fname, warn) {
	# read input file line by line, split by whitespaces, assign last substring to CpGo/e ratios
	# ... remove comments and trailing whitespaces
	v <- read.table(fname,  fill = TRUE, col.names = c("seq", "val"))
	obs <- v$val

	obs <- obs[!is.na(obs)]
	return(obs)
}


# process command line arguments
# expected arguments:
# - names of the species (each as a separate argument)
# - names of CpGo/e files of the species (each as a separate argument)
# ... parse arguments
option_list <- list(make_option(c("-o", "--frac-outl"), type = "double", default = 0.01,
                                help = "maximum fraction of CpGo/e ratios excluded as outliers [default %default]"),
                    make_option(c("-d", "--min-dist"), type = "double", default = 0.2,
                                help = "minimum distance between modes, modes that are closer are joined [default %default]"),
                    make_option(c("-c", "--conf-level"), type = "double", default = 0.95,
                                help = "level of the confidence intervals of the mode positions [default %default]"),
                    make_option(c("-m", "--mode-mass"), type = "double", default = 0.05,
                                help = "minimum probability mass of a mode [default %default]"),
                    make_option(c("-b", "--band-width"), type = "double", default = 1.06,
                                help = "bandwidth constant for kernels [default %default]"),
                    make_option(c("-B", "--bootstrap"), action="store_true", default = FALSE,
                                help = "calculate confidence intervals of mode positions using bootstrap [default %default]"),
                    make_option(c("-r", "--bootstrap-reps"), type = "integer", default = 1500,
                                help = "number of bootstrap repetitions [default %default]"),
                    make_option(c("-H", "--outlier-hist-file"), type = "character", default = "outliers_hist.pdf",
                                help = "name of the output file for the outlier histograms [default %default]"),
                    make_option(c("-C", "--cutoff-file"), type = "character", default = "outliers_cutoff.csv",
                                help = "name of the output file for the outlier cutoff [default %default]"),
                    make_option(c("-k", "--kde-file"), type = "character", default = "KDE.pdf",
                                help = "name of the output file for the KDE [default %default]"),
                    make_option(c("-p", "--peak-file"), type = "character", default = "peaks.csv",
                                help = "name of the output file describing the peaks of the KDE [default %default]"),
                    make_option(c("-s", "--bootstrap-file"), type = "character", default = "bootstrap.csv",
                                help = "name of the output file for the bootstrap results [default %default]"),
                    make_option(c("-f", "--no-warning-few-seqs"), action = "store_true", default = FALSE,
                                help = paste("suppress warning in case the input file only contains few values ",
                                             "[default %default]", sep = "")))

op <- OptionParser(usage = "notos.r [options] spc_name_1 ... spc_name_N CpGoe_file_name_1 ... CpGoe_file_name_N",
                   description = paste("\nDescription: Notos generates a histogram and a kernel density estimator from files containing CpGo/e ratios. ",
                                       "Moreover, it determines the number of modes of the CpGo/e ratio for each input file. The input files ",
                                       "can either be composed of \n",
                                       "1) CpGo/e ratios separated by linebreaks or\n",
                                       "2) sequence names and CpGo/e ratios with each sequence name put on a separate line together with its CpGo/e ratio ",
                                       "and sequence and CpGo/e being separated by whitespaces on each line.", sep = ""),
                   option_list = option_list)
args <- parse_args(op, positional_arguments = c(2, Inf))
num.args <- length(args$args)
use.bstrp <- args$options$`bootstrap`
supp.warn.few <- args$options$`no-warning-few-seqs`


# ... check number of arguments
# ... ... check number of mandatory arguments
if (num.args < 2) {
   stop("One species name and one file containing CpGo/e ratios have to be provided")
}

# ... ... check whether number of mandatory arguments is even
if (num.args %% 2 != 0) {
   stop("Number of arguments has to be even")
}

# ... ... check maximum fraction of CpGo/e ratios excluded as outliers
frac.outl <- args$options$`frac-outl`
if ((frac.outl < 0) || (frac.outl >= 1)) {
   stop("The maximum fraction of CpGo/e ratios excluded as outliers has to be greater or equal to zero and less than one")
}
if (frac.outl >= 0.2) {
   warning("The maximum fraction of CpGo/e ratios excluded as outliers has been set to a rather large value, resulting in the removal of many CpGo/e ratios")
}


# ... check numerical arguments
# ... ... check minimum distance between modes
min.dist <- args$options$`min-dist`
if (min.dist < 0) {
   stop("The minimum distance between modes has to be equal to or larger than zero")
}
if (min.dist >= 0.4) {
   warning("The minimum distance between modes has been set to a rather large value, resulting in a strong reduction of the number of modes")
}

# ... ... check confidence level
conf.lev <- args$options$`conf-level`
if ((conf.lev <= 0) || (conf.lev >= 1)) {
   stop("The level of the confidence intervals of the mode positions has to be larger than zero and smaller than one.")
}
if (conf.lev >= 0.995) {
   warning("The level of the confidence intervals of the mode positions has been set to a rather high value, resulting in very broad confidence intervals")
}

# ... ... check minimum probability mass of a mode
mode.mass <- args$options$`mode-mass`
if ((mode.mass < 0) || (mode.mass >= 1)) {
   stop("The minimum probability mass of a mode has to be larger than or equal to zero and smaller than one.")
}
if (mode.mass >= 0.3) {
   warning("The minimum probability mass of a mode has been set to a rather large value, resulting in the elemination of a high number of modes.")
}

# ... ... check bandwidth constant
band.width <- args$options$`band-width`
if (band.width <= 0) {
   stop("The bandwidth constant has to be positive")
}
if (band.width >= 5) {
   warning("The bandwidth constant has to been set to a rather large value, resulting in a strong smoothing")
}

# ... ... check number of boostrap repetitions
bstrp.reps <- args$options$`bootstrap-reps`
if (bstrp.reps != round(bstrp.reps)) {
  stop("The number of boostrap repetitions has to be a positive integer")
}
if (bstrp.reps <= 0) {
   stop("The number of boostrap repetitions has to be positive")
}
if (bstrp.reps >= 10000) {
   warning("The number of boostrap repetitions has been set to a rather large value, resulting in a long running time")
}

# ... check file name arguments
# ... ... check histogram output file name
outlier.hist.fname <- args$options$`outlier-hist-file`
if ( file.exists(outlier.hist.fname) && (file.info(outlier.hist.fname)$isdir) ) {
  stop(paste("File name for the outlier histogram output refers to a directory:", outlier.hist.fname))
}
v <- strsplit(outlier.hist.fname, split = ".", fixed = TRUE)[[1]]
if ((length(v) == 1) || (v[ length(v) ] != "pdf")) {
  warning(paste("File name for the outlier histogram output does not have a .pdf extension:", outlier.hist.fname))
}
g <- gregexpr(pattern ='/', outlier.hist.fname)[[1]]
if (as.vector(g)[1] != -1) {
  v <- as.vector(g)
  d <- substr(outlier.hist.fname, 1, v[length(v)])
  if (!file.exists(d)) {
    stop(paste("Path to file for the outlier histogram output is not valid:", outlier.hist.fname))
  }
}

# ... ... check outlier cutoff output file name
cutoff.fname <- args$options$`cutoff-file`
if ( file.exists(cutoff.fname) && (file.info(cutoff.fname)$isdir) ) {
  stop(paste("File name for the outlier cutoff table output refers to a directory:", cutoff.fname))
}
v <- strsplit(cutoff.fname, split = ".", fixed = TRUE)[[1]]
if (length(v) == 1) {
  stop(paste("File name for the outlier cutoff table output does not have a file extension:", cutoff.fname))
}
#if (v[ length(v) ] != "xlsx") {
#  warning(paste("File name for the outlier cutoff table output does not have a .xlsx extension:", cutoff.fname))
#}
g <- gregexpr(pattern ='/', cutoff.fname)[[1]]
if (as.vector(g)[1] != -1) {
  v <- as.vector(g)
  d <- substr(cutoff.fname, 1, v[length(v)])
  if (!file.exists(d)) {
    stop(paste("Path to file for the outlier cutoff is not valid:", cutoff.fname))
  }
}

# ... ... check KDE output file name
kde.fname <- args$options$`kde-file`
if ( file.exists(kde.fname) && (file.info(kde.fname)$isdir) ) {
  stop(paste("File name for the KDE output refers to a directory:", kde.fname))
}
v <- strsplit(kde.fname, split = ".", fixed = TRUE)[[1]]
if ((length(v) == 1) || (v[ length(v) ] != "pdf")) {
  warning(paste("File name for the KDE output does not have a .pdf extension:", kde.fname))
}
g <- gregexpr(pattern ='/', kde.fname)[[1]]
if (as.vector(g)[1] != -1) {
  v <- as.vector(g)
  d <- substr(kde.fname, 1, v[length(v)])
  if (!file.exists(d)) {
    stop(paste("Path to file for the KDE output is not valid:", kde.fname))
  }
}


# ... ... check peak descriptives output file name
peak.fname <- args$options$`peak-file`
if ( file.exists(peak.fname) && (file.info(peak.fname)$isdir) ) {
  stop(paste("File name for the peak descriptives refers to a directory:", peak.fname))
}
v <- strsplit(peak.fname, split = ".", fixed = TRUE)[[1]]
if ((length(v) == 1) || (v[ length(v) ] != "csv")) {
  warning(paste("File name for the peak descriptives does not have a .csv extension:", peak.fname))
}
g <- gregexpr(pattern ='/', peak.fname)[[1]]
if (as.vector(g)[1] != -1) {
  v <- as.vector(g)
  d <- substr(peak.fname, 1, v[length(v)])
  if (!file.exists(d)) {
    stop(paste("Path to file for the peak descriptives is not valid:", peak.fname))
  }
}

# ... ... check bootstrap results output file name
bstrp.fname <- args$options$`bootstrap-file`
if ( file.exists(bstrp.fname) && (file.info(bstrp.fname)$isdir) ) {
  stop(paste("File name for the bootstrap results refers to a directory:", bstrp.fname))
}
v <- strsplit(bstrp.fname, split = ".", fixed = TRUE)[[1]]
if ((length(v) == 1) || (v[ length(v) ] != "csv")) {
  warning(paste("File name for the bootstrap results does not have a .csv extension:", bstrp.fname))
}
g <- gregexpr(pattern ='/', bstrp.fname)[[1]]
if (as.vector(g)[1] != -1) {
  v <- as.vector(g)
  d <- substr(bstrp.fname, 1, v[length(v)])
  if (!file.exists(d)) {
    stop(paste("Path to file for the bootstrap results is not valid:", bstrp.fname))
  }
}


# ... ... check CpGo/e input file names
num.spec <- num.args / 2
spec.names <- args$args[1:num.spec]
cpgoe.fnames <- args$args[(num.spec + 1):num.args]
for (i in 1:length(cpgoe.fnames)) {
  if (!file.exists(cpgoe.fnames[i])) {
    stop(paste("CpGo/e file does not exist:", cpgoe.fnames[i]))
  }
  if (file.info(cpgoe.fnames[i])$isdir) {
    stop(paste("CpGo/e file name refers to a directory:", cpgoe.fnames[i]))
  }
}


# remove outliers and output histograms
# ... set up table with cutoff quantities
tab.des <- data.frame(matrix(NA, nrow = num.spec, ncol = 6))
names(tab.des) <- c("prop.zero", "prop.out.2iqr", "prop.out.3iqr",
                    "prop.out.4iqr", "prop.out.5iqr", "used")
rownames(tab.des) <- spec.names

# ... set up figure
t.height <- 6
t.width <- 20
pdf(outlier.hist.fname, height = t.height,width = t.width, paper = "special")
par(mfrow = c(1, 3), mgp = c(2, 0.5, 0), mar = c(4.0, 3.0, 1.5, 1))
tmp.fnames <- c()

# ... iterate through species
for (i in 1:num.spec) {
  fname <- cpgoe.fnames[i]

  obs <- read.CpGoe(fname, TRUE)

  # check CpGo/e ratios
  for (j in 1:length(obs)) {
    # is format legal?
    val <- as.numeric( obs[j] )
    err.str <- paste("Observation", i, "in", fname)
    if (!is.finite(val)) {
      stop(paste(err.str, "could not be converted to a number:", obs[j]))
    }

    # is ratio too small / large?
    if (val < 0) {
      stop(paste(err.str, "is negative:", val))
    } else {
      if (val > MAX.CPGOE) {
        warning(paste(err.str   , "is suspiciously large:", val, "\nthis value is replaced by", MAX.CPGOE))
      }
    }
  }

  # process outliers and store the results
  obs.org <- obs
  l <- proc.outliers(obs, frac.outl)
  if (!l[["valid"]]) {
    stop( paste("Too few values in", fname, "(less than 3) after removal of zeros"), call. = FALSE )
  }
  tab.des[i, "prop.zero"] <- l[["prop.zero"]]
  mu.obs <- l[["mu.obs"]]
  me.obs <- l[["me.obs"]]
  ul.mu <- l[["ul.mu"]]
  ll.mu <- l[["ll.mu"]]
  ul.me <- l[["ul.me"]]
  ll.me <- l[["ll.me"]]
  tab.des[i, "prop.out.2iqr"] <- l[["prop2"]]
  tab.des[i, "prop.out.3iqr"] <- l[["prop3"]]
  tab.des[i, "prop.out.4iqr"] <- l[["prop4"]]
  tab.des[i, "prop.out.5iqr"] <- l[["prop5"]]
  obs.cl <- l[["obs.cl"]]
  obs.nz <- l[["obs.nz"]]
  tab.des[i, "used"] <- l[["used"]]
  usedindex <- substr(l[["used"]],1,1)
  # Histograms
  # ... histogram 1: original data with zeros
  t.breaks <- seq(0, max(obs.org) + 1, by = 0.03)
  t.xlim <- c(0, ul.me["5"] + 0.1)
  hist(obs.org, breaks = t.breaks, xlim = t.xlim, xlab = "CpG o/e", main = "",
      sub = "Original data", prob = TRUE,
	  col = grey(0.9), border = grey(0.6))
  mtext(paste(spec.names[i]), side = 3, adj = 0)


  # ... histogram 3: median / iqr based
  t.lty <- rep(3, 4)
  t.lty[usedindex] <- 1

  hist(obs.nz, breaks = t.breaks, xlim = t.xlim, xlab = "CpG o/e", main = "",
      sub = "Data without zeros, Q1/3 +- k*IQR, k=2,...,4", prob = TRUE,
	  col = grey(0.9), border = grey(0.6))
  abline(v = me.obs, col = 'blue', lwd = 2)
  abline(v = c(ll.me, ul.me), col = "red", lty = rep(t.lty, 2))

  # ... histogram 4: cleaned data
  hist(obs.cl, breaks = t.breaks, xlim = t.xlim, xlab = "CpG o/e", main = "",
      sub = "Cleaned data", prob = TRUE,
	  col = grey(0.9), border = grey(0.6))
  abline(v = me.obs, col = 'blue', lwd = 2)
  abline(v = c(ll.me[usedindex], ul.me[usedindex]), col = "red")
}
invisible(dev.off())

# output cutoff quantities
write.table(tab.des, file = cutoff.fname, sep = "\t", col.names=NA)

# plot KDE and output quantities characterizing the peaks and the bootstrap results
# ... table with quantities characterizing the peaks
v <- col.names.peaks()
tab1.m <- data.frame(matrix(NA, nrow = num.spec, ncol = length(v)))
names(tab1.m) <- col.names.peaks()
rownames(tab1.m) <- spec.names

# ... table for the bootstrap
tab2.m <- data.frame(matrix(NA, nrow = num.spec, ncol = 7))
names(tab2.m) <- col.names.bs()
rownames(tab2.m) <- spec.names


# ... plotting
t.height <- 6
t.width <- 20
pdf(kde.fname, height = t.height,width = t.width, paper = "special")
for (i in 1:num.spec) {
  # read in GcGo/e ratios
  obs <- read.CpGoe(cpgoe.fnames[i], FALSE)
  l <- proc.outliers(obs, frac.outl)
  obs.cl <- l[["obs.cl"]]

  # check number of values
  fname <- cpgoe.fnames[i]
  if (length(obs.cl) < 3) {
    stop( paste("Too few values in", fname, "(less than 3) after removal of outliers and zeros"), call. = FALSE )
  }
  if (!supp.warn.few & length(obs.cl) < 250) {
    warning( paste(fname, " contains only few values (", length(obs.cl), ") after removal of outliers and zeros, which may lead to unreliable results", sep = ""), call. = FALSE )
  }

  # plotting
  l <- plot.KDE(obs.cl, t.name = spec.names[i], bs.cis = use.bstrp, bstrp.reps = bstrp.reps, conf.lev = conf.lev,
                min.dist = min.dist, mode.mass = mode.mass, band.width = band.width)
  tab1.m[i, ] <- l$tab.des
  if (use.bstrp) {
    tab2.m[i, ] <- l$tab.bs
  }
}
invisible(dev.off())
#sessionInfo()

# ... output quantities in tables
write.table(tab1.m, file = peak.fname, sep = "\t", col.names=NA)
if (use.bstrp) {
    write.table(tab2.m, file = bstrp.fname, sep = "\t", col.names=NA)
}
