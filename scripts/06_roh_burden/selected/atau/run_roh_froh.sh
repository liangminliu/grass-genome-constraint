#!/usr/bin/env bash
set -euo pipefail

############################################################
# Atau PLINK ROH/FROH final ROH/FROH calculation
#
# Purpose:
#   1. Do NOT rerun VCF filtering.
#   2. Reuse existing filtered VCF / PLINK bed.
#   3. Rerun only PLINK ROH with ROH parameters:
#        --homozyg-kb 100
#        --homozyg-gap 100
#        --homozyg-window-missing 5
#   4. Exclude extreme individual BW_01192 by default.
#
# Root:
#   /gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/6.pop/05.FROH
#
# Required existing files:
#   ROH_pipeline/Atau/01_filter/Atau141.filtered.bi-allelic.snp.vcf.gz
#   ROH_pipeline/Atau/04_froh/Lauto.from_filtered_SNPs.tsv
#
# Existing PLINK bed will be reused if present:
#   ROH_pipeline/Atau/02_plink/Atau141.bed
#   ROH_pipeline/Atau/02_plink/Atau141.bim
#   ROH_pipeline/Atau/02_plink/Atau141.fam
#
# New outputs:
#   ROH_pipeline/Atau/04_froh/plink_kb100_gap100_winmiss5_noBW01192/
############################################################

WORKDIR="/gpfs/hpc/home/kib_unit/lidezhu/liuliangmin/6.pop/05.FROH"
SPECIES="Atau"
PREFIX="Atau141"

# 默认去掉极端个体 BW_01192。
# 如果想保留所有个体，把这一行改成：
# EXCLUDE_IID=""
EXCLUDE_IID="BW_01192"

BASE="${WORKDIR}/ROH_pipeline/${SPECIES}"
FILTER_DIR="${BASE}/01_filter"
PLINK_DIR="${BASE}/02_plink"
FROH_DIR="${BASE}/04_froh"
LOG_DIR="${BASE}/logs"
SCRIPT_DIR="${BASE}/scripts"

FILTERED_VCF="${FILTER_DIR}/${PREFIX}.filtered.bi-allelic.snp.vcf.gz"
LAUTO_FILE="${FROH_DIR}/Lauto.from_filtered_SNPs.tsv"
PLINK_PREFIX="${PLINK_DIR}/${PREFIX}"

# 敏感性分析：最短 ROH 长度改为 100 kb
BASE_TAG="kb100_gap100_winmiss5"

if [[ -n "${EXCLUDE_IID}" ]]; then
    SAMPLE_TAG="no${EXCLUDE_IID//_/}"
else
    SAMPLE_TAG="all_samples"
fi

RUN_TAG="${BASE_TAG}_${SAMPLE_TAG}"
PLINK_OUT_DIR="${FROH_DIR}/plink_${RUN_TAG}"

mkdir -p "${PLINK_DIR}" "${PLINK_OUT_DIR}" "${LOG_DIR}" "${SCRIPT_DIR}"

DATE_TAG=$(date +%Y%m%d_%H%M%S)
LOG="${LOG_DIR}/${SPECIES}.rerun_PLINK_ROH.${RUN_TAG}.${DATE_TAG}.log"

exec > >(tee "${LOG}") 2>&1

echo "============================================================"
echo "[INFO] Job started: $(date)"
echo "[INFO] WORKDIR       = ${WORKDIR}"
echo "[INFO] BASE          = ${BASE}"
echo "[INFO] FILTERED_VCF  = ${FILTERED_VCF}"
echo "[INFO] PLINK_PREFIX  = ${PLINK_PREFIX}"
echo "[INFO] LAUTO_FILE    = ${LAUTO_FILE}"
echo "[INFO] RUN_TAG       = ${RUN_TAG}"
echo "[INFO] EXCLUDE_IID   = ${EXCLUDE_IID}"
echo "[INFO] Output dir    = ${PLINK_OUT_DIR}"
echo "============================================================"

############################################################
# 0. Check existing filtered data
############################################################

if [[ ! -s "${FILTERED_VCF}" ]]; then
    echo "[ERROR] Filtered VCF not found:"
    echo "        ${FILTERED_VCF}"
    echo "[ERROR] This script does not rerun VCF filtering."
    exit 1
fi

if [[ ! -s "${LAUTO_FILE}" ]]; then
    echo "[ERROR] Lauto file not found:"
    echo "        ${LAUTO_FILE}"
    echo "[ERROR] Please use the Lauto generated from the existing filtered SNPs."
    exit 1
fi

TOTAL_LAUTO=$(awk -F'\t' '$1=="TOTAL"{print $4}' "${LAUTO_FILE}")

if [[ -z "${TOTAL_LAUTO}" || "${TOTAL_LAUTO}" == "0" ]]; then
    echo "[ERROR] TOTAL Lauto is empty or 0 in:"
    echo "        ${LAUTO_FILE}"
    exit 1
fi

echo "[INFO] TOTAL_LAUTO = ${TOTAL_LAUTO}"

############################################################
# 1. Reuse existing PLINK bed, or convert from existing filtered VCF
#    This step does NOT refilter VCF.
############################################################

if [[ -s "${PLINK_PREFIX}.bed" && -s "${PLINK_PREFIX}.bim" && -s "${PLINK_PREFIX}.fam" ]]; then
    echo "[INFO] Existing PLINK bed/bim/fam found. Reusing:"
    echo "       ${PLINK_PREFIX}.bed"
    echo "       ${PLINK_PREFIX}.bim"
    echo "       ${PLINK_PREFIX}.fam"
else
    echo "[WARN] Existing PLINK bed/bim/fam not found."
    echo "[INFO] Converting existing filtered VCF to PLINK bed."
    echo "[INFO] This does NOT rerun VCF filtering."

    plink \
      --vcf "${FILTERED_VCF}" \
      --double-id \
      --allow-extra-chr \
      --keep-allele-order \
      --set-missing-var-ids '@:#:$1:$2' \
      --make-bed \
      --out "${PLINK_PREFIX}"
fi

############################################################
# 2. Prepare sample exclusion
############################################################

REMOVE_ARGS=()
REMOVE_FILE="${PLINK_OUT_DIR}/remove.${SAMPLE_TAG}.txt"
USED_SAMPLE_FILE="${PLINK_OUT_DIR}/${PREFIX}.${RUN_TAG}.samples.used.tsv"

if [[ -n "${EXCLUDE_IID}" ]]; then
    awk -v ex="${EXCLUDE_IID}" 'BEGIN{OFS="\t"} $2==ex{print $1,$2}' \
        "${PLINK_PREFIX}.fam" > "${REMOVE_FILE}"

    if [[ -s "${REMOVE_FILE}" ]]; then
        REMOVE_ARGS=(--remove "${REMOVE_FILE}")
        echo "[INFO] Excluding sample:"
        cat "${REMOVE_FILE}"

        awk -v ex="${EXCLUDE_IID}" 'BEGIN{OFS="\t"} $2!=ex{print $1,$2}' \
            "${PLINK_PREFIX}.fam" > "${USED_SAMPLE_FILE}"
    else
        echo "[WARN] EXCLUDE_IID=${EXCLUDE_IID} was not found in ${PLINK_PREFIX}.fam"
        echo "[WARN] No sample will be removed."

        EXCLUDE_IID=""
        SAMPLE_TAG="all_samples"
        RUN_TAG="${BASE_TAG}_${SAMPLE_TAG}"
        PLINK_OUT_DIR="${FROH_DIR}/plink_${RUN_TAG}"
        mkdir -p "${PLINK_OUT_DIR}"

        REMOVE_ARGS=()
        REMOVE_FILE="${PLINK_OUT_DIR}/remove.${SAMPLE_TAG}.txt"
        USED_SAMPLE_FILE="${PLINK_OUT_DIR}/${PREFIX}.${RUN_TAG}.samples.used.tsv"

        awk 'BEGIN{OFS="\t"}{print $1,$2}' \
            "${PLINK_PREFIX}.fam" > "${USED_SAMPLE_FILE}"
    fi
else
    echo "[INFO] No sample exclusion requested."
    awk 'BEGIN{OFS="\t"}{print $1,$2}' \
        "${PLINK_PREFIX}.fam" > "${USED_SAMPLE_FILE}"
fi

N_USED=$(wc -l < "${USED_SAMPLE_FILE}")
echo "[INFO] Number of samples used = ${N_USED}"
echo "[INFO] Used sample file:"
echo "       ${USED_SAMPLE_FILE}"

############################################################
# 3. Rerun PLINK ROH only
############################################################

ROH_PREFIX="${PLINK_DIR}/${PREFIX}.ROH.${RUN_TAG}"

echo "[INFO] Running PLINK ROH with ROH parameterss:"
echo "       --homozyg-snp 50"
echo "       --homozyg-kb 100"
echo "       --homozyg-density 50"
echo "       --homozyg-gap 100"
echo "       --homozyg-window-snp 50"
echo "       --homozyg-window-het 1"
echo "       --homozyg-window-missing 5"
echo "       --homozyg-window-threshold 0.05"

plink \
  --bfile "${PLINK_PREFIX}" \
  --allow-extra-chr \
  "${REMOVE_ARGS[@]}" \
  --homozyg \
  --homozyg-snp 50 \
  --homozyg-kb 100 \
  --homozyg-density 50 \
  --homozyg-gap 100 \
  --homozyg-window-snp 50 \
  --homozyg-window-het 1 \
  --homozyg-window-missing 5 \
  --homozyg-window-threshold 0.05 \
  --out "${ROH_PREFIX}"

############################################################
# 4. Parse PLINK ROH and calculate FROH
############################################################

echo "[INFO] Parsing PLINK ROH and calculating FROH..."

python3 - <<PY
import pandas as pd
from pathlib import Path

plink_dir = Path("${PLINK_DIR}")
out_dir = Path("${PLINK_OUT_DIR}")
lfile = Path("${LAUTO_FILE}")

run_tag = "${RUN_TAG}"
method = "plink.${RUN_TAG}"

hom_file = Path("${ROH_PREFIX}.hom")
sample_file = Path("${USED_SAMPLE_FILE}")

# Sensitivity threshold: 100 kb
min_roh_bp = 100000
min_snp = 50

seg_out = out_dir / f"ROH.segments.plink.ge100000.minSNP50.{run_tag}.classified.tsv"
froh_out = out_dir / f"FROH.final.plink.{run_tag}.tsv"
class_out = out_dir / f"FROH.by_class.plink.{run_tag}.tsv"
summary_out = out_dir / f"FROH.summary.plink.{run_tag}.tsv"

l = pd.read_csv(lfile, sep="\t")
Lauto = int(l.loc[l["chr"] == "TOTAL", "Lauto_bp"].iloc[0])

if Lauto <= 0:
    raise ValueError("Lauto <= 0")

samples = pd.read_csv(
    sample_file,
    sep=r"\s+",
    header=None,
    names=["FID", "IID"],
    dtype=str
)

empty_cols = [
    "FID","IID","sample","chr","snp1","snp2",
    "start_bp","end_bp","length_bp","length_kb",
    "nsnp","density","roh_class","method"
]

if not hom_file.exists() or hom_file.stat().st_size == 0:
    seg = pd.DataFrame(columns=empty_cols)
else:
    hom = pd.read_csv(hom_file, sep=r"\s+")

    if hom.empty:
        seg = pd.DataFrame(columns=empty_cols)
    else:
        seg = pd.DataFrame()
        seg["FID"] = hom["FID"].astype(str)
        seg["IID"] = hom["IID"].astype(str)
        seg["sample"] = hom["IID"].astype(str)
        seg["chr"] = hom["CHR"].astype(str)
        seg["snp1"] = hom["SNP1"].astype(str)
        seg["snp2"] = hom["SNP2"].astype(str)
        seg["start_bp"] = hom["POS1"].astype(int)
        seg["end_bp"] = hom["POS2"].astype(int)

        length_bp = seg["end_bp"] - seg["start_bp"] + 1
        fallback = (hom["KB"].astype(float) * 1000).round().astype(int)
        seg["length_bp"] = length_bp.where(length_bp > 0, fallback).astype(int)
        seg["length_kb"] = seg["length_bp"] / 1000.0
        seg["nsnp"] = hom["NSNP"].astype(int)

        if "DENSITY" in hom.columns:
            seg["density"] = hom["DENSITY"]
        else:
            seg["density"] = pd.NA

        # Keep only ROH segments consistent with --homozyg-kb 100 and --homozyg-snp 50.
        seg = seg[
            (seg["length_bp"] >= min_roh_bp) &
            (seg["nsnp"] >= min_snp)
        ].copy()

        def classify(x):
            if x < 250000:
                return "100-250kb"
            elif x < 500000:
                return "250-500kb"
            elif x < 1000000:
                return "500kb-1Mb"
            else:
                return ">1Mb"

        seg["roh_class"] = seg["length_bp"].apply(classify)
        seg["method"] = method

seg.to_csv(seg_out, sep="\t", index=False)

if seg.empty:
    g = pd.DataFrame(columns=["FID","IID","n_roh","sum_roh_bp"])
else:
    g = (
        seg.groupby(["FID","IID"], as_index=False)
        .agg(
            n_roh=("length_bp", "size"),
            sum_roh_bp=("length_bp", "sum")
        )
    )

out = samples.merge(g, on=["FID","IID"], how="left")
out["n_roh"] = out["n_roh"].fillna(0).astype(int)
out["sum_roh_bp"] = out["sum_roh_bp"].fillna(0).astype(int)
out["sum_roh_kb"] = out["sum_roh_bp"] / 1000.0
out["Lauto_bp"] = Lauto
out["FROH"] = out["sum_roh_bp"] / Lauto
out["method"] = method
out.to_csv(froh_out, sep="\t", index=False)

# Keep the four-class structure for plotting compatibility.
# With --homozyg-kb 100, the 100-250kb class is now informative.
classes = ["100-250kb", "250-500kb", "500kb-1Mb", ">1Mb"]
rows = []

for _, s in samples.iterrows():
    FID = s["FID"]
    IID = s["IID"]

    if seg.empty:
        sub = seg
    else:
        sub = seg[(seg["FID"] == FID) & (seg["IID"] == IID)]

    for c in classes:
        if sub.empty:
            n = 0
            bp = 0
        else:
            ss = sub[sub["roh_class"] == c]
            n = int(ss.shape[0])
            bp = int(ss["length_bp"].sum())

        rows.append({
            "FID": FID,
            "IID": IID,
            "sample": IID,
            "roh_class": c,
            "n_roh": n,
            "sum_roh_bp": bp,
            "sum_roh_kb": bp / 1000.0,
            "Lauto_bp": Lauto,
            "FROH": bp / Lauto,
            "method": method
        })

class_df = pd.DataFrame(rows)
class_df.to_csv(class_out, sep="\t", index=False)

summary = pd.DataFrame([{
    "n_sample": int(out.shape[0]),
    "min_roh_bp": min_roh_bp,
    "min_snp": min_snp,
    "Lauto_bp": Lauto,
    "mean_n_roh": float(out["n_roh"].mean()),
    "median_n_roh": float(out["n_roh"].median()),
    "mean_sum_roh_bp": float(out["sum_roh_bp"].mean()),
    "median_sum_roh_bp": float(out["sum_roh_bp"].median()),
    "mean_FROH": float(out["FROH"].mean()),
    "median_FROH": float(out["FROH"].median()),
    "min_FROH": float(out["FROH"].min()),
    "max_FROH": float(out["FROH"].max()),
    "sd_FROH": float(out["FROH"].std(ddof=1)),
    "method": method
}])
summary.to_csv(summary_out, sep="\t", index=False)

print(f"[INFO] PLINK ROH segments: {seg_out}")
print(f"[INFO] PLINK FROH:        {froh_out}")
print(f"[INFO] PLINK class FROH:  {class_out}")
print(f"[INFO] PLINK summary:     {summary_out}")
PY

echo "============================================================"
echo "[INFO] Job finished: $(date)"
echo "[INFO] New output directory:"
echo "       ${PLINK_OUT_DIR}"
echo ""
echo "[INFO] Main outputs:"
echo "  ${PLINK_OUT_DIR}/FROH.final.plink.${RUN_TAG}.tsv"
echo "  ${PLINK_OUT_DIR}/ROH.segments.plink.ge100000.minSNP50.${RUN_TAG}.classified.tsv"
echo "  ${PLINK_OUT_DIR}/FROH.by_class.plink.${RUN_TAG}.tsv"
echo "  ${PLINK_OUT_DIR}/FROH.summary.plink.${RUN_TAG}.tsv"
echo ""
echo "[INFO] PLINK raw ROH file:"
echo "  ${ROH_PREFIX}.hom"
echo "============================================================"