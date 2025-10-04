#!/bin/bash
set -euo pipefail

# Test matrix for PRECAL_DIFF vs UPDATE_FROM relative to correlation time
# Generates OU trajectories with known tau and runs MetaD with and without PRECAL_DIFF
# Then extracts the first hill sigma for comparison
# OU process parameters: tau=0.5 ps, sigma=0.15 nm timestep=0.005 ps
# metaD parameters: ADAPTIVE=DIFF, sigma=4000 (20ps), height=1.0, pace=5, temp=300, biasfactor=5 (PRECAL_DIFF)
# UPDATE_FROM use steps: 0 50 100 200 400 800 1000 2000 4000 8000
# correspoinding UPDATE_FROM use times: 0 0.25 0.5 1 2 4 5 10 20 40 ps

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

export LD_LIBRARY_PATH="/home/baizh/soft/plumed2/src/lib:${LD_LIBRARY_PATH:-}"
PLUMED_BIN="/home/baizh/soft/plumed2/src/lib/plumed"

gen_traj() {
  local tau=$1   # ps
  local out=$2
  python3 generate_trajectory.py --tau "$tau" --nsteps 10000 --dt 0.005 --sigma 0.15 --tau 0.5 --out "$out"
}

write_inputs() {
  local upd=$1
  local base=plumed_base_${upd}.dat
  local on=plumed_precal_on_${upd}.dat
  cat > "$base" <<EOF
d1: DISTANCE ATOMS=1,2

METAD ...
  ARG=d1
  ADAPTIVE=DIFF
  SIGMA=4000
  HEIGHT=1.0
  PACE=5
  TEMP=300
  BIASFACTOR=5
  UPDATE_FROM=${upd}
  FILE=HILLS_OFF_${upd}
  FMT=%14.8f
  LABEL=metad
... METAD

PRINT STRIDE=1 ARG=d1,metad.bias FILE=COLVAR_OFF_${upd} FMT=%14.8f

FLUSH STRIDE=50
EOF
  cat > "$on" <<EOF
d1: DISTANCE ATOMS=1,2

METAD ...
  PRECAL_DIFF
  ARG=d1
  ADAPTIVE=DIFF
  SIGMA=4000
  HEIGHT=1.0
  PACE=5
  TEMP=300
  BIASFACTOR=5
  UPDATE_FROM=${upd}
  FILE=HILLS_ON_${upd}
  FMT=%14.8f
  LABEL=metad
... METAD

PRINT STRIDE=1 ARG=d1,metad.bias FILE=COLVAR_ON_${upd} FMT=%14.8f

FLUSH STRIDE=50
EOF
}

run_pair() {
  local traj=$1
  local upd=$2
  "$PLUMED_BIN" driver --plumed "plumed_base_${upd}.dat" --timestep 0.005 --ixyz "$traj" > "run_off_${upd}.log" 2>&1 || true
  "$PLUMED_BIN" driver --plumed "plumed_precal_on_${upd}.dat" --timestep 0.005 --ixyz "$traj" > "run_on_${upd}.log" 2>&1 || true
}

first_sigma() {
  local file=$1
  awk '!/^#/ {print $3; exit}' "$file"
}

summarize_pair() {
  local upd=$1
  local off="HILLS_OFF_${upd}"
  local on="HILLS_ON_${upd}"
  local s_off="NA"
  local s_on="NA"
  if [[ -s "$off" ]]; then s_off=$(first_sigma "$off"); fi
  if [[ -s "$on" ]]; then s_on=$(first_sigma "$on"); fi
  printf "UPDATE_FROM=%-6s  OFF_sigma=%-10s  ON_sigma=%-10s\n" "$upd" "$s_off" "$s_on"
}

main() {
  echo "Running PRECAL_DIFF test suite..."
  echo ""
  # Define tau and two regimes of UPDATE_FROM relative to tau
  local tau_ps=0.2
  local traj="trajectory_ou_tau${tau_ps}.xyz"
  gen_traj "$tau_ps" "$traj"

  # UPDATE_FROM values to scan: < tau, ~= tau, > tau
  local update_steps=(0 50 100 200 400 800 1000 2000 4000 8000)
  local updates=()
  for step in "${update_steps[@]}"; do
    updates+=($(echo "$step * 0.005" | bc -l))
  done

  echo "Writing inputs..."
  for upd in "${updates[@]}"; do
    write_inputs "$upd"
  done

  echo "Running pairs (OFF vs ON)..."
  for upd in "${updates[@]}"; do
    rm -f HILLS_OFF_${upd} HILLS_ON_${upd}
    run_pair "$traj" "$upd"
  done
  echo ""
  echo "Summary (first hill sigma):"
  echo "---------------------------------------------"
  printf "%-16s %-16s %-16s\n" "UPDATE_FROM" "OFF_sigma" "ON_sigma"
  echo "---------------------------------------------"
  for upd in "${updates[@]}"; do
    summarize_pair "$upd"
  done | sed 's/^/  /'
  echo "---------------------------------------------"
  rm HILLS* COLVAR* run_*.log trajectory_ou_tau${tau_ps}.xyz plumed_*.dat
}

main "$@"


