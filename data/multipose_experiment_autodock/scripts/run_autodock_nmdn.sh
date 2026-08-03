#!/usr/bin/env bash
# Drive NMDN multi-pose scoring over the AutoDock poses until done or budget spent.
# nmdn_drive.py is resumable, so this just keeps re-entering it per target.
# Usage: BUDGET=55 run_autodock_nmdn.sh TP53 ESR1_ant PPARG
set -uo pipefail
ROOT=/home/joe/projects/drug_repurposing
PY=$HOME/miniconda3/envs/diffdock_nmdn/bin/python
BUDGET=${BUDGET:-55}
BATCH=${BATCH:-2000}
POSEDIR=${POSEDIR:-pose_sdfs_fixed}   # NMDN gets the valence-repaired poses, as production does   # keep modest: predict.py puts the whole batch in ONE GPU batch
cd "$ROOT"
END=$(( $(date +%s) + BUDGET * 60 ))

for T in "$@"; do
  while :; do
    LEFT=$(( (END - $(date +%s)) / 60 ))
    if [ "$LEFT" -lt 2 ]; then
      echo "=== budget spent, stopping before $T"; exit 0
    fi
    echo "=== $T  ($LEFT min left)  $(date +%H:%M:%S)"
    OUT=$("$PY" multipose_experiment/scripts/nmdn_drive.py "$T" \
            --work work_autodock --posedir "$POSEDIR" --minutes "$LEFT" --batch "$BATCH" 2>&1 | tail -4)
    echo "$OUT"
    case "$OUT" in
      *COMPLETE*) echo "=== $T done"; break ;;
    esac
  done
done
echo "=== finished $(date +%H:%M:%S)"
