#!/bin/bash
# Seed 232373165, present->2075 disturbance, protection STAGED over the 30 steps
# (161 new cells per step) instead of allocated one-shot at t=0.
#
#   bash scripts/run_future_staged.sh
#
# Same everything else as outputs/captain3_seed232373165_corrected_future, so the
# two are directly comparable and the only difference is the schedule.
set -u
PROJ="/Users/georgiagoldsmith/Library/CloudStorage/GoogleDrive-georgiagoldsmith@ucsb.edu/My Drive/Bren/CI Internship/Exploratory/case_study_ugf"
PY="/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain3preview-main/.venv/bin/python"
cd "$PROJ" || exit 1

export V2_DIR="/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain_testing/captain2-main/ugf_data_3km_corrected"
export SEED=232373165
export FUTURE_DISTURBANCE_FILE=disturbance_future_2075.tif
export TOTAL_TARGET_CELLS=13201
export CELLS_PER_STEP=161          # (13201 - 8374 locked-in PAs) / 30 steps
export HIDDEN_DIM=16,8
export RESULTS_DIR="$PROJ/outputs/captain3_seed232373165_corrected_future_staged"

"$PY" scripts/captain3_ugf_seed8052026.py || { echo "TRAIN FAILED"; exit 1; }
"$PY" scripts/captain3_ugf_predict.py     || { echo "PREDICT FAILED"; exit 1; }
echo "STAGED RUN DONE"
