## Overview

This folder contains data and code for an abstract submission to the 2027 MIT
Sloan Sports Analytics Conference (SSAC). This piece of work is regarding the
development of relative timing measures to provide a new level of analysis of
rowing biomechanical data.

## AI-assistance disclaimer

Claude (Anthropic's AI assistant, via Claude Code) was used throughout this project to:
- write and debug R scripts
- iterate the data-simulation method
- create figures
- run statistical analyses
- draft and improve conciseness of written text
All methodological choices, code, and results were directed,
reviewed, and verified step-by-step by the project author throughout;
Claude did not independently design the study or select its own outputs.
This disclaimer is included for transparency, consistent with disclosing
AI-assistance in research work.

## Data

The data analyzed are simulated from real athlete data due to organizational data 
governance requirements. The biomechanical data that were collected on elite and 
junior single scullers is not contained in this folder, however similar analyses 
have been successfully applied to this dataset. This folder contains the complete pipeline, 
from data simulation through to the final elite-vs-junior synchrony analysis.

## Analysis

Female junior vs elite single scull data comparison of relative timing 
between stroke (athlete right) and bow (athlete left) sides.

Two biomechanical measures:

- **Gate angle**: analyzed over the **full stroke cycle** (catch, drive, finish, 
recovery, catch). This is the position of the gate in degrees.
- **Gate force-X**: analyzed over the **drive phase only** (catch to
  finish) because the gate is only loaded during the drive so the recovery phase is
  dominated by sensor noise rather than a meaningful biomechanical signal. This 
  is the propulsive force measure originally recorded in kgf.

Synchrony magnitude is quantified via two methids:

- **Continuous Relative Phase (CRP)**, via the Hilbert transform.
- **Elastic registration (SRVF warping)**, via `fdasrvf::time_warping()`.

For both methods, the offset is defined as Stroke minus Bow: **negative
means Stroke reaches a landmark first (leads); positive means Bow does**.
The RMS of that offset (in degrees for CRP, unitless normalized-time for
warping) is the scalar synchrony-magnitude measure: smaller RMS means
tighter synchronization.

## Contents

- **`simulate_data.R`** -- generates `data_simulated/`. Included for
  transparency of method; will not run as-is here (see
  **Reproducibility**).
- **`data_simulated/{elite,junior}/`** -- 25 + 25 synthetic single-scull
  athletes, 20 analyzable strokes each, ~30 strokes/min.
- **`Nates_SSAC_rel_timing.R`** -- the analysis: loads
  `data_simulated/`, computes CRP and warping for both signals over their
  respective windows, runs the elite-vs-junior group comparison, and
  produces all figures. Fully self-contained and runnable as-is:
  ```
  Rscript rel_timing_simulated_mixed_window.R


## How the simulated data was made

`simulate_data.R` builds each synthetic athlete's strokes by
bootstrap-resampling real stroke *shapes* (with replacement) from a pool of
real elite/junior recordings, time-rescaling each to hit a per-athlete
target stroke rate (~N(30, 0.6) spm), then adding a small per-athlete
random scale/offset plus independent per-sample Gaussian noise. No real
recording is reproduced verbatim. Every generated file is re-verified
(reloaded through the same catch detector the analysis script uses)
before being written.

**Scope note**: only real female (W1X) data was available as a source, so
both simulated groups ("elite" and "junior") are drawn from real
elite/junior **female** single-scull recordings. This dataset does not
speak to any real sex-based comparison. Future analyses aim to include male
athletes.

## Reproducibility

`simulate_data.R` reads real de-identified recordings from a
`data_deidentified/` folder that is **not included** here (real,
non-public athlete data, even after de-identification) -- it documents
and verifies the generation method, it is not meant to be re-executed.
`Nates_SSAC_rel_timing.R`
has no such dependency and run standalone against the included
`data_simulated/`.

## Interpreting the results here

The elite-vs-junior comparison reproduces the same *direction* of effect
expected from the real analysis (tighter synchrony in the elite group),
since the real level-based stroke pool structure is preserved through the
resampling. However, because each synthetic athlete's strokes are
bootstrapped from the *whole* level's pooled real strokes rather than tied
to one real individual, real between-athlete variability is averaged away
-- so the effect sizes and p-values here are **larger and more
significant** than the real study's, and should not be read as estimates
of true population effect size or power. Treat this dataset and its
results as a demonstration of the analysis pipeline and its mechanics,
not as the scientific finding in their own right. However, real data analyses
can be presented on and published, although the raw biomechanical data cannot
be released.

