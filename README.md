# Inverse Target Trial Emulation (ITTE)

Reproducibility code for the revised manuscript:
**“Inverse Target Trial Emulation: A Method to Simulate Realistic Observational Data.”**

## Structure

- `HELP/main/` — final HELP analysis
- `HELP/supplementary/` — HELP validation and overlap diagnostics
- `HELP/src/` — shared HELP functions and DGP
- `GBSG2/main/` — final GBSG2 confounding analysis
- `GBSG2/immortal_time/` — final immortal-time analysis
- `GBSG2/supplementary/` — validation, fixed-final-n and q sensitivity
- `GBSG2/src/` — shared GBSG2 functions and DGP

The HELP analysis data file is included under `HELP/`.
GBSG2 is loaded directly from the `TH.data` R package.

Each analysis folder contains one `00_run_*.R` entry point.
