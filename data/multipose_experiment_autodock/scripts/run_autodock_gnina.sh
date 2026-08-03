#!/usr/bin/env bash
# GNINA multi-pose rescoring over the AutoDock poses, one target after another.
# Waits for any in-flight GNINA sweep first so containers don't oversubscribe the GPU.
# Usage: [JOBS=6] run_autodock_gnina.sh TP53 ESR1_ant PPARG
set -uo pipefail
ROOT=/home/joe/projects/drug_repurposing
JOBS=${JOBS:-6}
cd "$ROOT"

# Don't start while another sweep is mid-flight.
while ps -eo cmd | grep -q "run_gnina_max_throughput_sdf_[1]"; do
  echo "waiting for in-flight GNINA sweep... $(date +%H:%M:%S)"
  sleep 60
done

for T in "$@"; do
  n=$(find "multipose_experiment/work_autodock/$T/allposes" -name '*.sdf' 2>/dev/null | wc -l)
  done_n=$(find "multipose_experiment/work_autodock/$T/scored" -name '*_scored.sdf' 2>/dev/null | wc -l)
  if [ "$done_n" -ge "$n" ] && [ "$n" -gt 0 ]; then
    echo "=== $T already scored ($done_n/$n), skipping"; continue
  fi
  echo "=== GNINA $T  ($done_n/$n done)  $(date +%H:%M:%S)"
  WORKDIR=work_autodock bash multipose_experiment/scripts/run_gnina_multipose.sh "$T" "$JOBS" 2>&1 | tail -3
done
echo "=== gnina chain finished $(date +%H:%M:%S)"
