#| label: Figure_1_LOSS.R
# ======================================================
# Spore Survival — version without ENmisc (robust)
# ======================================================

## ---------- Helpers (no dependencies) ----------

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) approx(x = cw, y = x, xout = p, rule = 2, ties = "ordered")$y)
}

wbox_stats_one <- function(y, w) {
  qs <- weighted_quantile(y, w, c(0.25, 0.5, 0.75))
  q1 <- qs[1]; med <- qs[2]; q3 <- qs[3]
  if (!is.finite(q1) || !is.finite(q3)) {
    rng <- range(y, na.rm = TRUE)
    return(list(q1 = rng[1], med = median(y, na.rm = TRUE), q3 = rng[2],
                lo = rng[1], hi = rng[2]))
  }
  iqr <- q3 - q1
  if (!is.finite(iqr) || iqr == 0) {
    lo <- min(y, na.rm = TRUE); hi <- max(y, na.rm = TRUE)
    return(list(q1 = q1, med = med, q3 = q3, lo = lo, hi = hi))
  }
  data_min <- min(y, na.rm = TRUE); data_max <- max(y, na.rm = TRUE)
  lo <- max(data_min, q1 - 1.5 * iqr)
  hi <- min(data_max, q3 + 1.5 * iqr)
  list(q1 = q1, med = med, q3 = q3, lo = lo, hi = hi)
}

draw_wbox <- function(stats, xpos, boxw = 1.2, col = NA, border = "black", lwd = 2) {
  rect(xpos - boxw/2, stats$q1, xpos + boxw/2, stats$q3, col = col, border = border, lwd = lwd)
  segments(xpos - boxw/2, stats$med, xpos + boxw/2, stats$med, lwd = lwd, col = border)
  segments(xpos, stats$q3, xpos, stats$hi, lwd = lwd, col = border)
  segments(xpos, stats$q1, xpos, stats$lo, lwd = lwd, col = border)
  segments(xpos - boxw/4, stats$hi, xpos + boxw/4, stats$hi, lwd = lwd, col = border)
  segments(xpos - boxw/4, stats$lo, xpos + boxw/4, stats$lo, lwd = lwd, col = border)
}

wtd.boxplot <- function(formula, data = NULL, weights, at = NULL, add = FALSE,
                        xlim = NULL, ylim = NULL, xaxt = par("xaxt"),
                        xlab = "", ylab = "", col = NA, outcol = "black",
                        cex.lab = par("cex.lab"), boxw = 1.2, ...) {
  mf <- model.frame(formula, data = data)
  y <- mf[[1]]; x <- mf[[2]]
  
  if (!is.null(data) && is.character(weights)) weights <- data[[weights]]
  w <- weights; if (is.null(w)) w <- rep(1, length(y))
  
  x_key <- x
  grp_levels <- sort(unique(x_key))
  n_groups <- length(grp_levels)
  
  if (is.null(at)) pos <- if (is.numeric(grp_levels)) grp_levels else seq_len(n_groups)
  else {
    pos <- at
    if (length(pos) != n_groups)
      stop(sprintf("Length of 'at' (%d) must match number of groups (%d).", length(pos), n_groups))
  }
  
  if (!add) {
    if (is.null(xlim)) xlim <- range(pos) + c(-boxw, boxw)
    if (is.null(ylim)) ylim <- range(y, finite = TRUE)
    plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
         xaxt = xaxt, yaxt = par("yaxt"), type = "n", cex.lab = cex.lab)
  }
  
  if (length(col) == 1)    col    <- rep(col,    n_groups)
  if (length(outcol) == 1) outcol <- rep(outcol, n_groups)
  
  for (i in seq_len(n_groups)) {
    idx <- is.finite(y) & is.finite(w) & (w > 0)
    if (is.factor(x_key))       idx <- idx & (x_key == grp_levels[i])
    else if (is.character(x_key)) idx <- idx & (x_key == grp_levels[i])
    else                         idx <- idx & (round(as.numeric(x_key),6) == round(as.numeric(grp_levels[i]),6))
    
    if (!any(idx)) next
    stats <- wbox_stats_one(y[idx], w[idx])
    draw_wbox(stats, xpos = pos[i], boxw = boxw, col = col[i], border = outcol[i], lwd = 2)
  }
  invisible(pos)
}

mgp.axis.labels <- function(mgp_vec = c(2.5, 0.7, 0), type = c("x","y")) {
  type <- match.arg(type); old <- par(mgp = mgp_vec); invisible(old)
}
mgp.axis <- function(side, at, labels, cex.axis = 1, tick = FALSE, font = 1,
                     col.axis = "black", mar = NULL, ...) {
  if (!is.null(mar)) { old <- par(mar = mar); on.exit(par(old), add = TRUE) }
  axis(side = side, at = at, labels = labels, cex.axis = cex.axis,
       tick = tick, font = font, col.axis = col.axis, ...)
}

wmean_indexed <- function(x, y, w) {
  lv <- sort(unique(x))
  sapply(seq_along(lv), function(i) {
    idx <- is.finite(y) & is.finite(w) & (w > 0) &
      (if (is.factor(x)) x == lv[i]
       else if (is.character(x)) x == lv[i]
       else round(as.numeric(x),6) == round(as.numeric(lv[i]),6))
    if (!any(idx)) return(NA_real_)
    sum(y[idx] * w[idx]) / sum(w[idx])
  })
}

## ---------- Data ----------

JB32_MMN_R3   <- read.csv2("data/Spore_Survival_RAW_Values.csv")
JB32_MMN_R3F1 <- read.csv2("data/Spore_Survival_RAW_Values_F1.csv")

numify <- function(v) { if (is.factor(v)) as.numeric(as.character(v)) else suppressWarnings(as.numeric(v)) }
JB32_MMN_R3[,1:7]   <- lapply(JB32_MMN_R3[,1:7], numify)
JB32_MMN_R3F1[,1:7] <- lapply(JB32_MMN_R3F1[,1:7], numify)

JB32_MMN_R3_names_bold   <- c("11","24","30","60")
JB32_MMN_R3_names        <- c("\n\nn=212", "\n\nn=208", "\n\nn=208", "\n\nn=208")
JB32_MMN_R3F1_names_bold <- c("13","26","32","62")
JB32_MMN_R3F1_names      <- c("\n\nn=212","\n\nn=196", "\n\nn=204", "\n\nn=184")

at_R3   <- c(11,24,30,60)
at_R3F1 <- c(13,26,32,62)

## ---------- Plot Function ----------

draw_spore_survival_raw <- function() {
  par(mar = c(6.5, 5.0, 4.5, 2.1))
  
  # R3 boxplots (blanc)
  wtd.boxplot(JB32_MMN_R3[,5] ~ JB32_MMN_R3[,1],
              weights = JB32_MMN_R3[,6],
              xaxt = "n", cex.lab = 1.6,
              ylim = c(0.65, 1.1), xlim = c(10, 63),
              col = rgb(255,255,255, alpha=30, maxColorValue=255),
              outcol = rgb(0,0,0,0.35), ylab = "Spore Survival", xlab = "Days",
              at = at_R3, boxw = 1.2)
  
  title(main = "Spore Survival Raw Values", cex.main = 2.6, font.main = 2)
  
  # desplaçament de parelles al x
  offset_R3   <- -0.45
  offset_R3F1 <- +0.45
  
  mgp.axis.labels(c(3.0, 0.1, 0), type = "x")
  # Numeral gran a sobre
  mgp.axis(1, at = at_R3 + offset_R3, labels = JB32_MMN_R3_names_bold,
           cex.axis = 1.5, tick = FALSE, font = 2,
           col.axis = grey(0.65), line = 2)
  # n=... just sota
  mgp.axis(1, at = at_R3 + offset_R3, labels = JB32_MMN_R3_names,
           cex.axis = 0.9, tick = FALSE,
           col.axis = grey(0.80), line = 3)
  
  # punts bruts R3
  points(JB32_MMN_R3[,5] ~ JB32_MMN_R3[,7], col = rgb(0,0,0,0.25), pch = 16, cex = 0.75)
  
  ## --- Mitjanes R3 (manuals si existeixen; si no, ponderades)
  manual_R3 <- c(
    get0("MMN_JB32_R3_T2", ifnotfound = NA_real_),
    get0("MMN_JB32_R3_T3", ifnotfound = NA_real_),
    get0("MMN_JB32_R3_T4", ifnotfound = NA_real_),
    get0("MMN_JB32_R3_T5", ifnotfound = NA_real_)
  )
  if (all(is.finite(manual_R3))) {
    mean_R3 <- manual_R3
  } else {
    mean_R3 <- wmean_indexed(JB32_MMN_R3[,1], JB32_MMN_R3[,5], JB32_MMN_R3[,6])
  }
  lines(at_R3, mean_R3, type = "l", lwd = 3, col = "grey70")
  points(at_R3, mean_R3, pch = 4, cex = 1.8, lwd = 2, col = "darkblue")
  
  # R3F1 boxplots (beix)
  wtd.boxplot(JB32_MMN_R3F1[,5] ~ JB32_MMN_R3F1[,1],
              weights = JB32_MMN_R3F1[,6],
              xaxt = "n", cex.lab = 1.6,
              ylim = c(0.65, 1.1),
              col = "#F5F5DC", outcol = "black",
              ylab = "Spore Survival", xlab = "Days",
              at = at_R3F1, add = TRUE, boxw = 1.2)
  
  mgp.axis(1, at = at_R3F1 + offset_R3F1, labels = JB32_MMN_R3F1_names_bold,
           cex.axis = 1.5, tick = FALSE, font = 2,
           col.axis = grey(0.30), line = 2)
  mgp.axis(1, at = at_R3F1 + offset_R3F1, labels = JB32_MMN_R3F1_names,
           cex.axis = 0.9, tick = FALSE,
           col.axis = grey(0.45), line = 3)
  
  # punts bruts R3F1
  points(JB32_MMN_R3F1[,5] ~ JB32_MMN_R3F1[,7], col = rgb(0,0,0,0.35), pch = 1, cex = 0.85)
  
  ## --- Mitjanes R3F1 (manuals si existeixen; si no, ponderades)
  manual_R3F1 <- c(
    get0("MMN_JB32_R3F1_T2", ifnotfound = NA_real_),
    get0("MMN_JB32_R3F1_T3", ifnotfound = NA_real_),
    get0("MMN_JB32_R3F1_T4", ifnotfound = NA_real_),
    get0("MMN_JB32_R3F1_T5", ifnotfound = NA_real_)
  )
  if (all(is.finite(manual_R3F1))) {
    mean_R3F1 <- manual_R3F1
  } else {
    mean_R3F1 <- wmean_indexed(JB32_MMN_R3F1[,1], JB32_MMN_R3F1[,5], JB32_MMN_R3F1[,6])
  }
  lines(at_R3F1, mean_R3F1, type = "l", lwd = 3, lty = 3)
  points(at_R3F1, mean_R3F1, pch = 4, cex = 1.8, lwd = 2, col = "darkblue")
  
  # connectors verds
  for (i in seq_along(at_R3)) {
    if (is.finite(mean_R3[i]) && is.finite(mean_R3F1[i])) {
      lines(c(at_R3[i], at_R3F1[i]), c(mean_R3[i], mean_R3F1[i]),
            type = "l", lwd = 3, col = "#008000")
    }
  }
  
  # bracket de significació
  segments(61.5, 1.05, 62.5, 1.05)
  segments(61.5, 1.05, 61.5, 1.025)
  segments(62.5, 1.05, 62.5, 1.025)
  text(62, 1.07, labels = "**")
}
#Show in the RStudio window
draw_spore_survival_raw()