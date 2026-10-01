#!/bin/bash
#BSUB -J roadies97_accurate_96c_1T
#BSUB -q Q96C1T_X12
#BSUB -n 88
#BSUB -R "span[hosts=1]"
#BSUB -o /ds3200_1/users_root/liuliangmin/workspace/PHD/run97/logs/roadies97_accurate_96c_1T.%J.out
#BSUB -e /ds3200_1/users_root/liuliangmin/workspace/PHD/run97/logs/roadies97_accurate_96c_1T.%J.err

set -euo pipefail

WORKDIR="/ds3200_1/users_root/liuliangmin/workspace/PHD/run97"
CFG="${WORKDIR}/config_97sp_accurate.yaml"

OUT_DIR="${WORKDIR}/accurate_output_97sp"
ALL_OUT_DIR="${WORKDIR}/accurate_converge_files"
LOG_DIR="${WORKDIR}/logs"

mkdir -p "${LOG_DIR}"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate roadies_env

ROADIES_DIR="${CONDA_PREFIX}/ROADIES"
test -d "${ROADIES_DIR}" || { echo "[ERROR] ROADIES repo not found: ${ROADIES_DIR}"; exit 1; }

cd "${ROADIES_DIR}"

echo "[INFO] start : $(date)"
echo "[INFO] host  : $(hostname)"
echo "[INFO] env   : ${CONDA_PREFIX}"
echo "[INFO] repo  : $(pwd)"
echo "[INFO] cfg   : ${CFG}"

test -f "${CFG}" || { echo "[ERROR] missing config: ${CFG}"; exit 1; }

GENOME_DIR=$(python - <<'PY' "${CFG}"
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
print(cfg["GENOMES"])
PY
)

NUM_INSTANCES=$(python - <<'PY' "${CFG}"
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
print(cfg["NUM_INSTANCES"])
PY
)

test -d "${GENOME_DIR}" || { echo "[ERROR] missing genome dir: ${GENOME_DIR}"; exit 1; }

N_GENOMES=$(find -L "${GENOME_DIR}" -maxdepth 1 -type f \
  \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fasta" -o -name "*.fasta.gz" \) \
  | wc -l | awk '{print $1}')

echo "[INFO] genome dir         : ${GENOME_DIR}"
echo "[INFO] genome files found : ${N_GENOMES}"

if [[ "${N_GENOMES}" -ne 97 ]]; then
    echo "[ERROR] expected 97 genome files, found ${N_GENOMES}"
    echo "[INFO] detected files:"
    find -L "${GENOME_DIR}" -maxdepth 1 -type f \
      \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fasta" -o -name "*.fasta.gz" \) \
      | sort
    exit 1
fi

CORES="${LSB_DJOB_NUMPROC:-88}"
echo "[INFO] cores requested         : ${CORES}"
echo "[INFO] NUM_INSTANCES          : ${NUM_INSTANCES}"
echo "[INFO] approx threads/instance: $(( CORES / NUM_INSTANCES ))"

rm -rf "${OUT_DIR}" "${ALL_OUT_DIR}"
mkdir -p "${OUT_DIR}" "${ALL_OUT_DIR}"

snakemake --unlock >/dev/null 2>&1 || true

python run_roadies.py \
  --cores "${CORES}" \
  --mode accurate \
  --config "${CFG}"

echo "[INFO] done: $(date)"

echo "[INFO] iteration directories:"
find "${ALL_OUT_DIR}" -maxdepth 1 -type d -name "iteration_*" | sort -V || true

LATEST=$(find "${ALL_OUT_DIR}" -maxdepth 1 -type d -name "iteration_*" | sort -V | tail -n 1 || true)
if [[ -n "${LATEST}" ]]; then
    echo "[INFO] latest iteration: ${LATEST}"
    find "${LATEST}" -maxdepth 1 -type f | sort || true
fi