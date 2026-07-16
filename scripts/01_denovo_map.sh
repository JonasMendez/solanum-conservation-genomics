#!/usr/bin/env bash
# =============================================================================
# 01_denovo_map.sh
# Run denovo_map.pl (Stacks) on all samples.
# Parameters: m3, M2, n1 (previously validated)
#
# Stacks is locally installed (no AVX). All other tools via conda.
# Locally installed path: ~/local/stacks-2.68-noavx/bin
#
# Run from project root.
# =============================================================================
set -euo pipefail

export PATH="$HOME/local/stacks-2.68-noavx/bin:$PATH"

# ── parameters ────────────────────────────────────────────────────────────
RUN_ID="run_m3_M2_n1"
THREADS=8

# ── paths ─────────────────────────────────────────────────────────────────
# FASTQs live in fastq/ in the project root (demultiplexed, named <sample_id>.fastq.gz)
SAMPLES="fastq"
POPMAP="metadata/popmaps/popmap_all_taxa_species.tsv"
OUTDIR="01_stacks/runs/${RUN_ID}"
LOGDIR="01_stacks/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

# ── verify inputs ─────────────────────────────────────────────────────────
[[ -d "${SAMPLES}" ]] || { echo "ERROR: input dir not found: ${SAMPLES}"; exit 1; }
[[ -s "${POPMAP}" ]]  || { echo "ERROR: popmap not found: ${POPMAP}"; exit 1; }

N_FASTQ=$(ls "${SAMPLES}"/*.fastq.gz 2>/dev/null | wc -l)
N_POPMAP=$(awk 'NF' "${POPMAP}" | wc -l)
echo "Samples in input dir : ${N_FASTQ}"
echo "Samples in popmap    : ${N_POPMAP}"
echo "Run ID               : ${RUN_ID}"
echo "Parameters           : m=3, M=2, n=1"
echo ""

# ── run ───────────────────────────────────────────────────────────────────
denovo_map.pl \
    --samples "${SAMPLES}" \
    --popmap  "${POPMAP}" \
    -o        "${OUTDIR}" \
    -T        "${THREADS}" \
    -m 3 -M 2 -n 1 \
    2>&1 | tee "${LOGDIR}/denovo_map_${RUN_ID}.log"

echo ""
echo "denovo_map complete: ${OUTDIR}"
echo "Review log: ${LOGDIR}/denovo_map_${RUN_ID}.log"
