#!/usr/bin/env bash
# Explode AutoDock DLGs (all 10 runs/ligand) into per-pose SDFs for the 3 control targets.
set -euo pipefail
ROOT=/home/joe/projects/drug_repurposing
PY=$HOME/miniconda3/envs/diffdock_nmdn/bin/python
cd "$ROOT"
for T in "$@"; do
  echo "=== $T $(date +%H:%M:%S)"
  "$PY" multipose_experiment/scripts/explode_dlg_poses.py "$T" 10
done
echo "=== all done $(date +%H:%M:%S)"
