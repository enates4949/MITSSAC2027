#!/usr/bin/env Rscript
#
# Nates_SSAC_rel_timing.R
#
# Analyzes BOTH force and gate angle signals, but each
#     analyzed over the window that's actually appropriate for it: gate
#     ANGLE over the FULL STROKE CYCLE (catch to next catch, drive +
#     recovery -- angle is a real, well-behaved motion throughout), and
#     gate FORCE-X restricted to the DRIVE PHASE ONLY (catch to finish --
#     during recovery there is no real load on the gate, so force there
#     is mostly sensor noise, same reasoning as the original pipeline).
#
# Method (Hilbert-transform CRP, SRVF warping with center_warpings=FALSE,
# RMS-based synchrony magnitude, log-scale Welch/Wilcoxon/Hedges' g group
# comparison) is otherwise identical to rel_timing_simulated.R.
#
# Output goes to mixed_window/ (figures/ and results/ subfolders)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  match_idx <- grep(file_flag, cmd_args)
  if (length(match_idx) > 0) return(dirname(normalizePath(sub(file_flag, "", cmd_args[match_idx[1]]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path)) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
}

BASE_DIR <- get_script_dir()
DATA_DIR <- file.path(BASE_DIR, "data_simulated")
# This copy lives inside the self-contained full_stroke/ submission
# package, so figures/results are direct siblings, not nested under an
# extra "mixed_window" subfolder the way the original (in sim/ and
# to_submit/ top level) is.
OUT_DIR  <- file.path(BASE_DIR, "figures")
RESULTS_DIR <- file.path(BASE_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(DATA_DIR)) stop("data_simulated/ not found next to this script -- run simulate_data.R first.")

# ---------------------------------------------------------------------------
# Load + quality screen (identical to rel_timing_simulated.R)
# ---------------------------------------------------------------------------

list_athlete_files <- function(group) sort(list.files(file.path(DATA_DIR, group), pattern = "\\.csv$", full.names = TRUE))
anon_id <- function(fp) tools::file_path_sans_ext(basename(fp))
find_stroke_indices <- function(norm_time) {
  n <- length(norm_time); idx <- c(1L)
  if (n > 1) for (i in 2:n) if (isTRUE(norm_time[i] < norm_time[i - 1])) idx <- c(idx, i)
  c(idx, n + 1L)
}
load_and_check <- function(fp) {
  athlete <- anon_id(fp)
  d <- tryCatch(read.csv(fp), error = function(e) NULL)
  if (is.null(d)) return(list(status = "load_failed", athlete = athlete))
  sg_a <- grep("^Stroke_GateAngle_1_", colnames(d), value = TRUE)[1]
  bg_a <- grep("^Bow_GateAngle_1_", colnames(d), value = TRUE)[1]
  sg_f <- grep("^Stroke_GateForceX_1_", colnames(d), value = TRUE)[1]
  bg_f <- grep("^Bow_GateForceX_1_", colnames(d), value = TRUE)[1]
  nt_col <- grep("^Normalised.Time$", colnames(d), value = TRUE)[1]
  required <- c(sg_a, bg_a, sg_f, bg_f, nt_col)
  if (any(is.na(required))) return(list(status = "missing_columns", athlete = athlete))
  nt <- suppressWarnings(as.numeric(d[[nt_col]]))
  angle_s <- suppressWarnings(as.numeric(d[[sg_a]])); angle_b <- suppressWarnings(as.numeric(d[[bg_a]]))
  force_s <- suppressWarnings(as.numeric(d[[sg_f]])); force_b <- suppressWarnings(as.numeric(d[[bg_f]]))
  keep <- is.finite(nt) & is.finite(angle_s) & is.finite(angle_b) & is.finite(force_s) & is.finite(force_b)
  nt <- nt[keep]; angle_s <- angle_s[keep]; angle_b <- angle_b[keep]; force_s <- force_s[keep]; force_b <- force_b[keep]
  if (length(nt) < 100) return(list(status = "too_short", athlete = athlete))
  if (identical(angle_s, angle_b)) return(list(status = "duplicate_angle_columns", athlete = athlete))
  if (identical(force_s, force_b)) return(list(status = "duplicate_force_columns", athlete = athlete))
  stroke_idx <- find_stroke_indices(nt); n_strokes <- length(stroke_idx) - 1
  if (n_strokes < 3) return(list(status = "too_few_strokes", athlete = athlete))
  list(status = "ok", athlete = athlete, n_strokes = n_strokes, stroke_idx = stroke_idx,
       angle_s = angle_s, angle_b = angle_b, force_s = force_s, force_b = force_b)
}

FS <- 50

open_png_retry <- function(path, ..., max_tries = 3, wait_sec = 2) {
  for (attempt in seq_len(max_tries)) {
    ok <- tryCatch({ grDevices::png(path, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    if (attempt < max_tries) Sys.sleep(wait_sec)
  }
  stop("could not open PNG device for ", path, " after ", max_tries, " attempts")
}

# ---------------------------------------------------------------------------
# CRP (Hilbert) + warping (SRVF) 
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(signal)
  library(fdasrvf)
})

MIN_STROKES <- 20
DESIRED_LENGTH <- 100

butter_lowpass_filter <- function(data, cutoff = 6, fs = FS, order = 4) {
  bf <- signal::butter(order, cutoff / (0.5 * fs), type = "low")
  as.numeric(signal::filtfilt(bf, data))
}
hilbert_transform <- function(x) {
  n <- length(x); X <- fft(x); h <- rep(0, n)
  if (n %% 2 == 0) { h[1] <- 1; h[n / 2 + 1] <- 1; h[2:(n / 2)] <- 2 } else { h[1] <- 1; h[2:((n + 1) / 2)] <- 2 }
  fft(X * h, inverse = TRUE) / n
}
compute_hilbert_phase <- function(sig) {
  centered <- sig - mean(sig, na.rm = TRUE)
  an <- hilbert_transform(centered)
  atan2(-Im(an), Re(an))
}
wrap_pi <- function(x) (x + pi) %% (2 * pi) - pi
offset_summary <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(mae = NA_real_, rms = NA_real_, ce = NA_real_, ve = NA_real_, dominance = NA_real_))
  ce <- mean(x)
  c(mae = mean(abs(x)), rms = sqrt(mean(x^2)), ce = ce, ve = sqrt(mean((x - ce)^2)), dominance = mean(sign(x)))
}
register_pair <- function(curve1, curve2, mean_time) {
  warp <- suppressMessages(fdasrvf::time_warping(cbind(curve1, curve2), mean_time, center_warpings = FALSE))
  list(gamma1 = warp$warping_functions[, 1], gamma2 = warp$warping_functions[, 2])
}

#' Drive-phase catch/finish windows -- used directly for FORCE (recovery
#' excluded, since it's sensor noise there), and as the source of robust
#' catch positions for the full-stroke ANGLE windows below.
find_drive_windows <- function(r, n_strokes = MIN_STROKES, edge_tolerance = 0.15, catch_pad = 5) {
  angle_ref <- (butter_lowpass_filter(r$angle_s) + butter_lowpass_filter(r$angle_b)) / 2
  windows <- list(); n_scanned <- 0L; n_skipped <- 0L; i <- 2L
  while (length(windows) < n_strokes && i < length(r$stroke_idx)) {
    nominal_start <- r$stroke_idx[i]; seg_end <- r$stroke_idx[i + 1] - 1L; n_scanned <- n_scanned + 1L
    pad_lo <- max(1L, nominal_start - catch_pad); pad_hi <- min(seg_end, nominal_start + catch_pad)
    if (pad_hi > pad_lo && (seg_end - nominal_start + 1L) >= 10) {
      catch_abs <- pad_lo + which.min(angle_ref[pad_lo:pad_hi]) - 1L
      drive_region <- angle_ref[catch_abs:seg_end]; idx_max <- which.max(drive_region)
      finish_abs <- catch_abs + idx_max - 1L; seg_len <- seg_end - catch_abs + 1L
      ok <- idx_max > 1 && idx_max <= (1 - edge_tolerance) * seg_len
      if (ok) windows[[length(windows) + 1]] <- c(catch_abs, finish_abs) else n_skipped <- n_skipped + 1L
    } else n_skipped <- n_skipped + 1L
    i <- i + 1L
  }
  list(windows = windows, n_scanned = n_scanned, n_skipped = n_skipped)
}

#' FULL stroke cycle windows (catch to next catch -- drive + recovery),
#' used for ANGLE only. Built from n_strokes+1 real catch positions
#' (reusing find_drive_windows()'s own robust catch detection).
find_full_stroke_windows <- function(r, n_strokes = MIN_STROKES) {
  dw <- find_drive_windows(r, n_strokes = n_strokes + 1L)
  if (length(dw$windows) < n_strokes + 1L) {
    return(list(windows = list(), n_scanned = dw$n_scanned, n_skipped = dw$n_skipped))
  }
  catches <- sapply(dw$windows, function(w) w[1])
  windows <- lapply(seq_len(n_strokes), function(k) c(catches[k], catches[k + 1] - 1L))
  list(windows = windows, n_scanned = dw$n_scanned, n_skipped = dw$n_skipped)
}

#' Shared CRP+warping computation, unchanged from rel_timing_simulated.R
#' -- takes windows as a parameter, so the SAME function works whether
#' those windows are drive-phase (force) or full-stroke (angle).
compute_athlete_offsets <- function(r, windows, signal = c("angle", "force")) {
  signal <- match.arg(signal)
  sig_s <- if (signal == "angle") r$angle_s else r$force_s
  sig_b <- if (signal == "angle") r$angle_b else r$force_b
  n_show <- length(windows)
  sig_s_filt <- butter_lowpass_filter(sig_s); sig_b_filt <- butter_lowpass_filter(sig_b)
  phase_s_full <- compute_hilbert_phase(sig_s_filt); phase_b_full <- compute_hilbert_phase(sig_b_filt)

  mean_time <- seq(0, 1, length.out = DESIRED_LENGTH)
  crp_metric_mat <- matrix(NA_real_, n_show, 5, dimnames = list(NULL, names(offset_summary(1))))
  warp_metric_mat <- crp_metric_mat
  crp_curves <- matrix(NA_real_, DESIRED_LENGTH, n_show); warp_curves <- matrix(NA_real_, DESIRED_LENGTH, n_show)

  for (k in seq_len(n_show)) {
    start <- windows[[k]][1]; end <- windows[[k]][2]
    crp_native <- wrap_pi(phase_s_full[start:end] - phase_b_full[start:end])
    crp_metric_mat[k, ] <- offset_summary(crp_native)
    crp_curves[, k] <- approx(seq(0, 1, length.out = length(crp_native)), crp_native, xout = mean_time)$y

    c1 <- approx(seq(0, 1, length.out = end - start + 1), sig_s[start:end], xout = mean_time)$y
    c2 <- approx(seq(0, 1, length.out = end - start + 1), sig_b[start:end], xout = mean_time)$y
    reg <- register_pair(c1, c2, mean_time)
    warp_curves[, k] <- reg$gamma1 - reg$gamma2
    warp_metric_mat[k, ] <- offset_summary(warp_curves[, k])
  }
  list(crp_metrics = colMeans(crp_metric_mat, na.rm = TRUE), warp_metrics = colMeans(warp_metric_mat, na.rm = TRUE),
       crp_curves = crp_curves, warp_curves = warp_curves, mean_time = mean_time, n_strokes = n_show)
}

# ---------------------------------------------------------------------------
# Visualization helpers
# ---------------------------------------------------------------------------

plot_raw_grid <- function(good_list, windows_list, signal_label, window_label, group_label, out_file, sig_col_s, sig_col_b) {
  n <- length(good_list); if (n == 0) return(invisible(NULL))
  ncol_grid <- ceiling(sqrt(n * 1.4)); nrow_grid <- ceiling(n / ncol_grid)
  open_png_retry(out_file, width = 3.2 * ncol_grid, height = 2.6 * nrow_grid, units = "in", res = 130)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nrow_grid, ncol_grid), mar = c(2, 2, 2.6, 0.5), oma = c(0, 0, 3, 0))
  for (r in good_list) {
    windows <- windows_list[[r$athlete]]$windows
    win <- unlist(lapply(windows, function(w) w[1]:w[2]))
    catch_positions <- sapply(windows, function(w) w[1])
    rate_spm <- if (length(catch_positions) > 1) 60 * FS / mean(diff(catch_positions)) else NA
    title_txt <- sprintf("%s\n%.1f spm (%d strokes shown)", r$athlete, rate_spm, length(windows))
    plot(seq_along(win), r[[sig_col_s]][win], type = "l", col = "steelblue", xlab = "", ylab = "", main = title_txt, cex.main = 0.8)
    lines(seq_along(win), r[[sig_col_b]][win], col = "darkorange")
  }
  mtext(paste0(group_label, ": raw Stroke (blue) vs Bow (orange), ", signal_label, " -- ", window_label), outer = TRUE, cex = 1.1, font = 2)
}

extract_metric <- function(off_list, method = c("crp", "warp"), metric = "rms") {
  method <- match.arg(method); key <- paste0(method, "_metrics")
  vals <- sapply(off_list, function(m) m[[key]][metric]); names(vals) <- names(off_list); vals
}
cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * sd(x)^2 + (ny - 1) * sd(y)^2) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp
}
hedges_g_ci <- function(x, y, conf_level = 0.95) {
  nx <- length(x); ny <- length(y); df <- nx + ny - 2
  d <- cohens_d(x, y); J <- 1 - 3 / (4 * df - 1); g <- d * J
  se_g <- sqrt(J^2 * ((nx + ny) / (nx * ny) + d^2 / (2 * df)))
  z <- qnorm(1 - (1 - conf_level) / 2)
  c(g = g, lo = g - z * se_g, hi = g + z * se_g)
}
compare_groups <- function(elite_vals, junior_vals, label) {
  log_elite <- log(elite_vals); log_junior <- log(junior_vals)
  tt <- t.test(log_elite, log_junior)
  wt <- suppressWarnings(wilcox.test(elite_vals, junior_vals, conf.int = TRUE))
  g_ci <- hedges_g_ci(log_elite, log_junior)
  log_ratio_est <- unname(tt$estimate[1] - tt$estimate[2])
  data.frame(
    comparison = label, n_elite = length(elite_vals), n_junior = length(junior_vals),
    mean_elite = mean(elite_vals), sd_elite = sd(elite_vals), mean_junior = mean(junior_vals), sd_junior = sd(junior_vals),
    welch_t = unname(tt$statistic), welch_df = unname(tt$parameter), welch_p = tt$p.value,
    geomean_ratio = exp(log_ratio_est), geomean_ratio_ci_lo = exp(tt$conf.int[1]), geomean_ratio_ci_hi = exp(tt$conf.int[2]),
    wilcox_W = unname(wt$statistic), wilcox_p = wt$p.value, wilcox_hl_estimate = unname(wt$estimate),
    wilcox_ci_lo = unname(wt$conf.int[1]), wilcox_ci_hi = unname(wt$conf.int[2]),
    hedges_g = unname(g_ci["g"]), hedges_g_ci_lo = unname(g_ci["lo"]), hedges_g_ci_hi = unname(g_ci["hi"]),
    stringsAsFactors = FALSE
  )
}

#' Compact, publication-ready version of the comparison table: CRP rows
#' converted from radians to degrees (the more common CRP convention,
#' and what the companion figure already uses) -- ratio, p, and Hedges'
#' g are left unconverted since they're unit-independent (ratios and
#' standardized effect sizes cancel units out, verified: converting the
#' underlying values to degrees changes mean_elite/sd_elite/etc. but
#' leaves geomean_ratio/hedges_g numerically identical). Includes
#' Hedges' g's 95% CI in the saved CSV either way -- whether that CI
#' also appears in the abstract text itself is a separate space/clutter
#' call, not a computation one.
save_abstract_table <- function(comparison_df, out_file) {
  rad2deg <- function(x) x * 180 / pi
  is_crp <- grepl("^CRP", comparison_df$comparison)
  df <- comparison_df
  for (col in c("mean_elite", "sd_elite", "mean_junior", "sd_junior")) df[[col]][is_crp] <- rad2deg(df[[col]][is_crp])
  out <- data.frame(
    Comparison = df$comparison,
    Elite_mean = round(df$mean_elite, 3), Elite_SD = round(df$sd_elite, 3),
    Junior_mean = round(df$mean_junior, 3), Junior_SD = round(df$sd_junior, 3),
    Ratio = round(df$geomean_ratio, 3), Ratio_CI_lo = round(df$geomean_ratio_ci_lo, 3), Ratio_CI_hi = round(df$geomean_ratio_ci_hi, 3),
    p = signif(df$welch_p, 3),
    Hedges_g = round(df$hedges_g, 3), Hedges_g_CI_lo = round(df$hedges_g_ci_lo, 3), Hedges_g_CI_hi = round(df$hedges_g_ci_hi, 3)
  )
  write.csv(out, out_file, row.names = FALSE)
  cat("  Saved:", out_file, "\n")
}

plot_offset_overlay_grid <- function(off_list, curve_field, group_label, signal_label, window_label, method_label, out_file) {
  ids <- names(off_list); n <- length(ids); if (n == 0) return(invisible(NULL))
  ncol_grid <- ceiling(sqrt(n * 1.4)); nrow_grid <- ceiling(n / ncol_grid)
  open_png_retry(out_file, width = 3.2 * ncol_grid, height = 2.6 * nrow_grid, units = "in", res = 130)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nrow_grid, ncol_grid), mar = c(2, 2, 2, 0.5), oma = c(0, 0, 3, 0))
  curves_key <- paste0(curve_field, "_curves")
  for (id in ids) {
    m <- off_list[[id]]; curves <- m[[curves_key]]
    matplot(m$mean_time, curves, type = "l", lty = 1, col = grDevices::adjustcolor("steelblue", 0.3), xlab = "", ylab = "", main = id, cex.main = 0.85)
    lines(m$mean_time, rowMeans(curves, na.rm = TRUE), col = "darkorange", lwd = 2)
    abline(h = 0, lty = 2, col = "grey50")
  }
  mtext(paste0(group_label, ": ", method_label, " offset (Stroke - Bow), ", signal_label, " -- ", window_label,
               " -- ", off_list[[1]]$n_strokes, " strokes overlaid per athlete, mean in orange"), outer = TRUE, cex = 1.0, font = 2)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  elite_loaded  <- lapply(list_athlete_files("elite"), load_and_check)
  junior_loaded <- lapply(list_athlete_files("junior"), load_and_check)
  filter_good <- function(loaded) {
    good <- Filter(function(r) r$status == "ok", loaded)
    names(good) <- sapply(good, function(r) r$athlete); good
  }
  elite_good  <- filter_good(elite_loaded)
  junior_good <- filter_good(junior_loaded)
  cat(sprintf("ELITE: %d/%d loaded OK | JUNIOR: %d/%d loaded OK\n",
              length(elite_good), length(elite_loaded), length(junior_good), length(junior_loaded)))

  # ANGLE: full-stroke windows (catch to next catch). FORCE: drive-phase
  # windows (catch to finish) -- each signal gets the window that's
  # actually appropriate for it, computed independently.
  elite_angle_windows  <- lapply(elite_good,  find_full_stroke_windows)
  junior_angle_windows <- lapply(junior_good, find_full_stroke_windows)
  elite_force_windows  <- lapply(elite_good,  find_drive_windows)
  junior_force_windows <- lapply(junior_good, find_drive_windows)

  report_windows <- function(windows_list, label) {
    ok_all <- all(sapply(windows_list, function(w) length(w$windows) == MIN_STROKES))
    if (ok_all) cat(sprintf("  %s: all athletes found %d/%d windows.\n", label, MIN_STROKES, MIN_STROKES))
    else for (id in names(windows_list)) {
      w <- windows_list[[id]]
      if (length(w$windows) != MIN_STROKES) cat(sprintf("  %s %-16s found %d/%d windows\n", label, id, length(w$windows), MIN_STROKES))
    }
  }
  cat("=== Window detection ===\n")
  cat(" -- angle (full stroke) --\n"); report_windows(elite_angle_windows, "ELITE"); report_windows(junior_angle_windows, "JUNIOR")
  cat(" -- force (drive phase) --\n"); report_windows(elite_force_windows, "ELITE"); report_windows(junior_force_windows, "JUNIOR")

  plot_raw_grid(elite_good,  elite_angle_windows,  "Gate Angle", "FULL STROKE (catch to next catch)", "Simulated Elite W1X",
                file.path(OUT_DIR, "mixed_window_raw_elite_angle.png"), "angle_s", "angle_b")
  plot_raw_grid(junior_good, junior_angle_windows, "Gate Angle", "FULL STROKE (catch to next catch)", "Simulated Junior W1X",
                file.path(OUT_DIR, "mixed_window_raw_junior_angle.png"), "angle_s", "angle_b")
  plot_raw_grid(elite_good,  elite_force_windows,  "Gate Force X", "DRIVE PHASE ONLY (catch to finish)", "Simulated Elite W1X",
                file.path(OUT_DIR, "mixed_window_raw_elite_force.png"), "force_s", "force_b")
  plot_raw_grid(junior_good, junior_force_windows, "Gate Force X", "DRIVE PHASE ONLY (catch to finish)", "Simulated Junior W1X",
                file.path(OUT_DIR, "mixed_window_raw_junior_force.png"), "force_s", "force_b")

  cat("\n=== CRP + warping ===\n")
  elite_angle_off  <- lapply(names(elite_good),  function(id) compute_athlete_offsets(elite_good[[id]],  elite_angle_windows[[id]]$windows,  "angle"))
  names(elite_angle_off) <- names(elite_good)
  junior_angle_off <- lapply(names(junior_good), function(id) compute_athlete_offsets(junior_good[[id]], junior_angle_windows[[id]]$windows, "angle"))
  names(junior_angle_off) <- names(junior_good)
  elite_force_off  <- lapply(names(elite_good),  function(id) compute_athlete_offsets(elite_good[[id]],  elite_force_windows[[id]]$windows,  "force"))
  names(elite_force_off) <- names(elite_good)
  junior_force_off <- lapply(names(junior_good), function(id) compute_athlete_offsets(junior_good[[id]], junior_force_windows[[id]]$windows, "force"))
  names(junior_force_off) <- names(junior_good)
  cat(sprintf("  angle (full stroke): %d elite + %d junior | force (drive phase): %d elite + %d junior\n",
              length(elite_angle_off), length(junior_angle_off), length(elite_force_off), length(junior_force_off)))

  # Cache full offset objects (curves included) for the companion
  # plot_group_mean_offset_mixed_window.R figure script.
  saveRDS(list(elite_angle_off = elite_angle_off, junior_angle_off = junior_angle_off,
               elite_force_off = elite_force_off, junior_force_off = junior_force_off),
          file.path(RESULTS_DIR, "mixed_window_offsets.rds"))

  to_df <- function(off_list, group) {
    do.call(rbind, lapply(names(off_list), function(id) {
      m <- off_list[[id]]
      row <- data.frame(athlete = id, group = group, stringsAsFactors = FALSE)
      for (nm in names(m$crp_metrics))  row[[paste0("crp_",  nm)]] <- m$crp_metrics[nm]
      for (nm in names(m$warp_metrics)) row[[paste0("warp_", nm)]] <- m$warp_metrics[nm]
      row
    }))
  }
  write.csv(rbind(to_df(elite_angle_off, "elite"), to_df(junior_angle_off, "junior")),
            file.path(RESULTS_DIR, "mixed_window_angle_summary.csv"), row.names = FALSE)
  write.csv(rbind(to_df(elite_force_off, "elite"), to_df(junior_force_off, "junior")),
            file.path(RESULTS_DIR, "mixed_window_force_summary.csv"), row.names = FALSE)
  cat("  Saved: mixed_window_angle_summary.csv, mixed_window_force_summary.csv\n")

  cat("\n=== Elite vs junior RMS comparison ===\n")
  comparison_df <- rbind(
    compare_groups(extract_metric(elite_angle_off, "crp"),  extract_metric(junior_angle_off, "crp"),  "CRP x angle (full stroke)"),
    compare_groups(extract_metric(elite_force_off, "crp"),  extract_metric(junior_force_off, "crp"),  "CRP x force (drive phase)"),
    compare_groups(extract_metric(elite_angle_off, "warp"), extract_metric(junior_angle_off, "warp"), "Warp x angle (full stroke)"),
    compare_groups(extract_metric(elite_force_off, "warp"), extract_metric(junior_force_off, "warp"), "Warp x force (drive phase)")
  )
  comparison_df$welch_p_holm <- p.adjust(comparison_df$welch_p, method = "holm")
  comparison_df$wilcox_p_holm <- p.adjust(comparison_df$wilcox_p, method = "holm")
  print(comparison_df, row.names = FALSE, digits = 4)
  comparison_out <- file.path(RESULTS_DIR, "mixed_window_group_comparison.csv")
  write.csv(comparison_df, comparison_out, row.names = FALSE)
  cat("  Saved:", comparison_out, "\n")
  save_abstract_table(comparison_df, file.path(RESULTS_DIR, "mixed_window_abstract_table.csv"))

  plot_offset_overlay_grid(elite_angle_off,  "crp",  "Simulated Elite W1X",  "Gate Angle",   "FULL STROKE",  "CRP (Hilbert)", file.path(OUT_DIR, "mixed_window_crp_elite_angle.png"))
  plot_offset_overlay_grid(junior_angle_off, "crp",  "Simulated Junior W1X", "Gate Angle",   "FULL STROKE",  "CRP (Hilbert)", file.path(OUT_DIR, "mixed_window_crp_junior_angle.png"))
  plot_offset_overlay_grid(elite_force_off,  "crp",  "Simulated Elite W1X",  "Gate Force X", "DRIVE PHASE",  "CRP (Hilbert)", file.path(OUT_DIR, "mixed_window_crp_elite_force.png"))
  plot_offset_overlay_grid(junior_force_off, "crp",  "Simulated Junior W1X", "Gate Force X", "DRIVE PHASE",  "CRP (Hilbert)", file.path(OUT_DIR, "mixed_window_crp_junior_force.png"))
  plot_offset_overlay_grid(elite_angle_off,  "warp", "Simulated Elite W1X",  "Gate Angle",   "FULL STROKE",  "Warping",       file.path(OUT_DIR, "mixed_window_warp_elite_angle.png"))
  plot_offset_overlay_grid(junior_angle_off, "warp", "Simulated Junior W1X", "Gate Angle",   "FULL STROKE",  "Warping",       file.path(OUT_DIR, "mixed_window_warp_junior_angle.png"))
  plot_offset_overlay_grid(elite_force_off,  "warp", "Simulated Elite W1X",  "Gate Force X", "DRIVE PHASE",  "Warping",       file.path(OUT_DIR, "mixed_window_warp_elite_force.png"))
  plot_offset_overlay_grid(junior_force_off, "warp", "Simulated Junior W1X", "Gate Force X", "DRIVE PHASE",  "Warping",       file.path(OUT_DIR, "mixed_window_warp_junior_force.png"))

  cat("\nDone. Figures in:", OUT_DIR, "\nResults in:", RESULTS_DIR, "\n")
  invisible(list(elite_angle_off = elite_angle_off, junior_angle_off = junior_angle_off,
                 elite_force_off = elite_force_off, junior_force_off = junior_force_off))
}

if (sys.nframe() == 0) main()
