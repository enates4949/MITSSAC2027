#!/usr/bin/env Rscript
#
# simulate_data.R
#
# Generates a fully synthetic, GitHub-safe dataset for rel_timing_simulated.R, 
# using data in data_deindentified/  (not included for data security reasons).

# Every output stroke is resampled to a new duration, given a per-athlete random
# scale/offset, and has independent Gaussian noise injected on top.
#
# SCOPE: W1X only (elite + junior female), matching rel_timing_female.R
#
# Design:
#   1. Load the same data_deidentified/{elite,junior} real recordings used by
#      rel_timing_female.R (same load_and_check()/find_drive_windows(), kept
#      standalone/duplicated here per this project's one-script-one-file
#      convention -- see rel_timing_female.R for the validated originals).
#   2. For each real "good" athlete, extract 20 FULL stroke cycles (catch to
#      next catch, i.e. drive + recovery) with each stroke's own catch/finish
#      position -- not just the drive-phase sub-window -- so the simulated
#      files are ordinary continuous recordings the existing loader can
#      re-segment itself, exactly like a real file.
#   3. Pool all real full-stroke cycles within a level (elite: 14 athletes x
#      20 = 280 strokes; junior: 15 x 20 = 300 strokes) into two banks.
#   4. For each of the 2 groups (elite, junior), generate 25 synthetic athletes: 20 strokes each,
#      bootstrap-resampled (with replacement) from the matching level's bank,
#      time-resampled to a per-athlete target rate ~N(30, 0.6) spm (matching
#      the real spread verified earlier in this project) with small per-
#      stroke jitter, given one athlete-level random scale + offset per
#      signal type (applied identically to Stroke and Bow, so the real
#      Stroke-Bow offset relationship inside each borrowed stroke is
#      preserved, not distorted), and finally perturbed with independent
#      per-sample Gaussian noise per channel.
#   5. Write each synthetic athlete as a minimal CSV (Normalised.Time +
#      the 4 gate columns, same naming convention as the real files) to
#      data_simulated/<group>/.
#   6. Re-load every written file with the SAME load_and_check()/
#      find_drive_windows() and require exactly 20/20 valid drive phases,
#      zero skips -- hard-verified, not assumed.
#
# Required packages: signal (for the Butterworth filter find_drive_windows()
# uses internally -- same as rel_timing_female.R).

set.seed(2025)  # fixed seed -- reruns reproduce byte-identical simulated data

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
REAL_DATA_DIR <- file.path(BASE_DIR, "data_deidentified")
SIM_DATA_DIR  <- file.path(BASE_DIR, "data_simulated")

if (!dir.exists(REAL_DATA_DIR)) stop("data_deidentified/ not found next to this script -- needed as the real-shape source for simulation.")

FS <- 50  # sampling rate, Hz -- verified throughout this project

# ---------------------------------------------------------------------------
# Load + quality screen -- duplicated verbatim from rel_timing_female.R
# (kept standalone per this project's one-script-one-file convention).
# ---------------------------------------------------------------------------

list_athlete_files <- function(group) {
  sort(list.files(file.path(REAL_DATA_DIR, group), pattern = "\\.csv$", full.names = TRUE))
}

anon_id <- function(fp) tools::file_path_sans_ext(basename(fp))

find_stroke_indices <- function(norm_time) {
  n <- length(norm_time)
  idx <- c(1L)
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
  if (any(is.na(required))) {
    missing <- c("Stroke_GateAngle", "Bow_GateAngle", "Stroke_GateForceX", "Bow_GateForceX", "Normalised.Time")[is.na(required)]
    return(list(status = "missing_columns", athlete = athlete, detail = paste(missing, collapse = ", ")))
  }

  nt <- suppressWarnings(as.numeric(d[[nt_col]]))
  angle_s <- suppressWarnings(as.numeric(d[[sg_a]])); angle_b <- suppressWarnings(as.numeric(d[[bg_a]]))
  force_s <- suppressWarnings(as.numeric(d[[sg_f]])); force_b <- suppressWarnings(as.numeric(d[[bg_f]]))
  keep <- is.finite(nt) & is.finite(angle_s) & is.finite(angle_b) & is.finite(force_s) & is.finite(force_b)
  nt <- nt[keep]; angle_s <- angle_s[keep]; angle_b <- angle_b[keep]; force_s <- force_s[keep]; force_b <- force_b[keep]
  if (length(nt) < 100) return(list(status = "too_short", athlete = athlete, detail = paste0(length(nt), " usable rows")))

  if (identical(angle_s, angle_b)) return(list(status = "duplicate_angle_columns", athlete = athlete, detail = "export duplication fault"))
  if (identical(force_s, force_b)) return(list(status = "duplicate_force_columns", athlete = athlete, detail = "export duplication fault"))

  stroke_idx <- find_stroke_indices(nt)
  n_strokes <- length(stroke_idx) - 1
  if (n_strokes < 3) return(list(status = "too_few_strokes", athlete = athlete, detail = paste0(n_strokes, " strokes detected")))

  list(status = "ok", athlete = athlete, n_strokes = n_strokes,
       stroke_idx = stroke_idx, angle_s = angle_s, angle_b = angle_b, force_s = force_s, force_b = force_b)
}

suppressPackageStartupMessages(library(signal))

butter_lowpass_filter <- function(data, cutoff = 6, fs = FS, order = 4) {
  bf <- signal::butter(order, cutoff / (0.5 * fs), type = "low")
  as.numeric(signal::filtfilt(bf, data))
}

find_drive_windows <- function(r, n_strokes = 20L, edge_tolerance = 0.15, catch_pad = 5) {
  angle_ref <- (butter_lowpass_filter(r$angle_s) + butter_lowpass_filter(r$angle_b)) / 2
  windows <- list()
  n_scanned <- 0L; n_skipped <- 0L
  i <- 2L  # skip stroke_idx[1]:stroke_idx[2]-1, the fractional first segment
  while (length(windows) < n_strokes && i < length(r$stroke_idx)) {
    nominal_start <- r$stroke_idx[i]; seg_end <- r$stroke_idx[i + 1] - 1L
    n_scanned <- n_scanned + 1L
    pad_lo <- max(1L, nominal_start - catch_pad)
    pad_hi <- min(seg_end, nominal_start + catch_pad)
    if (pad_hi > pad_lo && (seg_end - nominal_start + 1L) >= 10) {
      catch_abs <- pad_lo + which.min(angle_ref[pad_lo:pad_hi]) - 1L
      drive_region <- angle_ref[catch_abs:seg_end]
      idx_max <- which.max(drive_region)
      finish_abs <- catch_abs + idx_max - 1L
      seg_len <- seg_end - catch_abs + 1L
      ok <- idx_max > 1 && idx_max <= (1 - edge_tolerance) * seg_len
      if (ok) windows[[length(windows) + 1]] <- c(catch_abs, finish_abs) else n_skipped <- n_skipped + 1L
    } else {
      n_skipped <- n_skipped + 1L
    }
    i <- i + 1L
  }
  list(windows = windows, n_scanned = n_scanned, n_skipped = n_skipped)
}

# ---------------------------------------------------------------------------
# Step 1: build the real per-level stroke-template banks
# ---------------------------------------------------------------------------

#' For one real "good" athlete, extract 20 FULL stroke cycles (catch to next
#' catch -- drive + recovery), each tagged with its own drive-phase end
#' (finish_rel, relative to the stroke's own start). Needs 21 valid drive-
#' window catches (to bound 20 full cycles) -- verified below to always be
#' available given the real data's minimum stroke counts (elite min 23,
#' junior min 33 total strokes; find_drive_windows() found 0 skips out of
#' 20 requested for every real athlete earlier in this project, so 21 has
#' ample headroom). Hard stop()s rather than silently truncating if not.
extract_full_stroke_templates <- function(r) {
  windows21 <- find_drive_windows(r, n_strokes = 21L)$windows
  if (length(windows21) < 21L) {
    stop("extract_full_stroke_templates: ", r$athlete, " only yielded ", length(windows21),
         " valid drive windows (need 21 to bound 20 full cycles) -- investigate before simulating from it.")
  }
  lapply(1:20, function(k) {
    full_start <- windows21[[k]][1]
    full_end   <- windows21[[k + 1]][1] - 1L
    finish_rel <- windows21[[k]][2] - full_start + 1L  # 1-indexed position of finish within this segment
    list(angle_s = r$angle_s[full_start:full_end], angle_b = r$angle_b[full_start:full_end],
         force_s = r$force_s[full_start:full_end], force_b = r$force_b[full_start:full_end],
         finish_rel = finish_rel, n = full_end - full_start + 1L)
  })
}

build_bank <- function(group) {
  loaded <- lapply(list_athlete_files(group), load_and_check)
  good <- Filter(function(r) r$status == "ok", loaded)
  cat(sprintf("  %s real template pool: %d/%d athletes usable\n", group, length(good), length(loaded)))
  templates <- unlist(lapply(good, extract_full_stroke_templates), recursive = FALSE)
  cat(sprintf("  %s bank: %d real full-stroke templates (%d athletes x 20)\n", group, length(templates), length(good)))
  templates
}

cat("=== Building real stroke-template banks from data_deidentified/ ===\n")
elite_bank  <- build_bank("elite")
junior_bank <- build_bank("junior")
cat("\n")

# Pooled per-channel SD, used to scale injected noise/offsets relative to
# real variability rather than picking arbitrary absolute units.
pooled_sd <- function(bank, fields) sd(unlist(lapply(bank, function(t) unlist(t[fields]))))
elite_angle_sd  <- pooled_sd(elite_bank,  c("angle_s", "angle_b"))
elite_force_sd  <- pooled_sd(elite_bank,  c("force_s", "force_b"))
junior_angle_sd <- pooled_sd(junior_bank, c("angle_s", "angle_b"))
junior_force_sd <- pooled_sd(junior_bank, c("force_s", "force_b"))

# ---------------------------------------------------------------------------
# Step 2: simulate one synthetic athlete
# ---------------------------------------------------------------------------

N_STROKES_ANALYSIS <- 20L
# find_drive_windows() ALWAYS discards the first detected stroke segment,
# assuming (correctly, for real recordings) that it's fractional -- the
# recording started mid-stroke. Our simulated files start exactly at a
# catch, so their first stroke is actually clean -- and would still get
# thrown away by that same unconditional skip. Write one extra "throwaway"
# leading stroke so the standard, UNCHANGED find_drive_windows(n_strokes=20)
# still ends up with exactly 20 usable strokes, matching the real-data
# convention exactly rather than special-casing simulated files.
N_STROKES <- N_STROKES_ANALYSIS + 1L
TARGET_RATE_MEAN <- 30    # spm -- matches this project's verified real elite spread (~30.2 +/- 0.6)
TARGET_RATE_SD <- 0.6     # between-athlete spread
STROKE_RATE_JITTER_SD <- 0.3  # small within-athlete, stroke-to-stroke jitter around that athlete's own rate
NOISE_FRAC <- 0.03        # per-sample injected noise SD, as a fraction of the pooled real channel SD
ATHLETE_SCALE_SD <- 0.04  # per-athlete multiplicative amplitude jitter (applied identically to Stroke+Bow)
ATHLETE_OFFSET_FRAC <- 0.06  # per-athlete additive offset, as a fraction of the pooled channel SD

#' Bootstrap-resample 20 strokes from `bank` into one synthetic athlete's
#' continuous recording. Per-athlete scale/offset is drawn ONCE per signal
#' type (angle, force) and applied identically to Stroke and Bow, so the
#' real within-stroke Stroke-Bow relationship each borrowed stroke carries
#' is preserved, not distorted -- only the independent per-sample noise
#' differs between Stroke and Bow, as real measurement noise would.
simulate_one_athlete <- function(bank, angle_sd, force_sd) {
  idxs <- sample(seq_along(bank), N_STROKES, replace = TRUE)
  athlete_rate <- rnorm(1, TARGET_RATE_MEAN, TARGET_RATE_SD)
  angle_scale <- rnorm(1, 1, ATHLETE_SCALE_SD); angle_offset <- rnorm(1, 0, ATHLETE_OFFSET_FRAC * angle_sd)
  force_scale <- rnorm(1, 1, ATHLETE_SCALE_SD); force_offset <- rnorm(1, 0, ATHLETE_OFFSET_FRAC * force_sd)

  out_nt <- numeric(0); out_as <- numeric(0); out_ab <- numeric(0); out_fs <- numeric(0); out_fb <- numeric(0)
  stroke_rates_used <- numeric(N_STROKES)

  for (k in seq_len(N_STROKES)) {
    tmpl <- bank[[idxs[k]]]
    stroke_rate <- max(26, min(34, rnorm(1, athlete_rate, STROKE_RATE_JITTER_SD)))
    stroke_rates_used[k] <- stroke_rate
    n_target <- max(20L, round(FS * 60 / stroke_rate))

    old_t <- seq(0, 1, length.out = tmpl$n); new_t <- seq(0, 1, length.out = n_target)
    as_v <- approx(old_t, tmpl$angle_s, xout = new_t)$y; ab_v <- approx(old_t, tmpl$angle_b, xout = new_t)$y
    fs_v <- approx(old_t, tmpl$force_s, xout = new_t)$y; fb_v <- approx(old_t, tmpl$force_b, xout = new_t)$y

    as_v <- as_v * angle_scale + angle_offset + rnorm(n_target, 0, NOISE_FRAC * angle_sd)
    ab_v <- ab_v * angle_scale + angle_offset + rnorm(n_target, 0, NOISE_FRAC * angle_sd)
    fs_v <- fs_v * force_scale + force_offset + rnorm(n_target, 0, NOISE_FRAC * force_sd)
    fb_v <- fb_v * force_scale + force_offset + rnorm(n_target, 0, NOISE_FRAC * force_sd)

    nt_seg <- (0:(n_target - 1)) * (100 / n_target)  # matches real Normalised.Time convention: 0 up to just under 100, resets next stroke

    out_nt <- c(out_nt, nt_seg); out_as <- c(out_as, as_v); out_ab <- c(out_ab, ab_v); out_fs <- c(out_fs, fs_v); out_fb <- c(out_fb, fb_v)
  }
  # (achieved rate is verified properly below, from reloading the written
  # file through the real catch-detection logic -- not estimated here.)
  list(nt = out_nt, angle_s = out_as, angle_b = out_ab, force_s = out_fs, force_b = out_fb)
}

# ---------------------------------------------------------------------------
# Step 3: generate all 4 groups (elite, junior x force, angle), write CSVs, verify
# ---------------------------------------------------------------------------

N_PER_GROUP <- 25L
GROUPS <- list(
  elite  = list(bank = elite_bank,  angle_sd = elite_angle_sd,  force_sd = elite_force_sd),
  junior = list(bank = junior_bank, angle_sd = junior_angle_sd, force_sd = junior_force_sd)
)

cat("=== Generating simulated athletes ===\n")
generated_files <- character(0)
for (group in names(GROUPS)) {
  cfg <- GROUPS[[group]]
  group_dir <- file.path(SIM_DATA_DIR, group)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)
  for (n in seq_len(N_PER_GROUP)) {
    athlete_id <- sprintf("sim_%s_%02d", group, n)
    sim <- simulate_one_athlete(cfg$bank, cfg$angle_sd, cfg$force_sd)
    df <- setNames(
      data.frame(sim$nt, sim$angle_s, sim$angle_b, sim$force_s, sim$force_b),
      c("Normalised.Time",
        paste0("Stroke_GateAngle_1_", athlete_id), paste0("Bow_GateAngle_1_", athlete_id),
        paste0("Stroke_GateForceX_1_", athlete_id), paste0("Bow_GateForceX_1_", athlete_id))
    )
    out_fp <- file.path(group_dir, paste0(athlete_id, ".csv"))
    write.csv(df, out_fp, row.names = FALSE)
    generated_files <- c(generated_files, out_fp)
  }
  cat(sprintf("  %-14s %d athletes written to %s\n", group, N_PER_GROUP, group_dir))
}
cat(sprintf("\nTotal: %d simulated athletes across %d groups.\n\n", length(generated_files), length(GROUPS)))

# ---------------------------------------------------------------------------
# Step 4: hard verification -- reload every generated file with the SAME
# loader/window-detector the real pipeline uses, require exactly 20/20
# valid catch-to-finish drive phases with zero skips. Never assumed.
# ---------------------------------------------------------------------------

cat("=== Verifying every simulated file: reload + 20/20 drive-phase detection ===\n")
achieved_rates <- numeric(0)
failures <- character(0)
for (fp in generated_files) {
  r <- load_and_check(fp)
  if (r$status != "ok") { failures <- c(failures, sprintf("%s: load_and_check status=%s", basename(fp), r$status)); next }
  w <- find_drive_windows(r, n_strokes = N_STROKES_ANALYSIS)
  if (length(w$windows) != N_STROKES_ANALYSIS || w$n_skipped != 0) {
    failures <- c(failures, sprintf("%s: found %d/%d windows, %d skipped", basename(fp), length(w$windows), N_STROKES_ANALYSIS, w$n_skipped))
    next
  }
  catch_positions <- sapply(w$windows, function(x) x[1])
  achieved_rates <- c(achieved_rates, 60 * FS / mean(diff(catch_positions)))
}

if (length(failures) > 0) {
  cat(sprintf("  %d file(s) FAILED verification:\n", length(failures)))
  for (f in failures) cat("    ", f, "\n")
  stop("simulate_data.R: verification failed for ", length(failures), " file(s) -- see above. Not safe to proceed to rel_timing_simulated.R.")
} else {
  cat(sprintf("  All %d simulated files: exactly %d/%d valid drive phases, zero skips.\n", length(generated_files), N_STROKES_ANALYSIS, N_STROKES_ANALYSIS))
}
cat(sprintf("  Achieved stroke rate across all simulated athletes: mean=%.2f spm, sd=%.2f, range=%.1f-%.1f\n",
            mean(achieved_rates), sd(achieved_rates), min(achieved_rates), max(achieved_rates)))
cat("\nDone. Simulated data in:", SIM_DATA_DIR, "\n")
