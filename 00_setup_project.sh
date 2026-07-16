#!/usr/bin/env bash
# =============================================================================
# 00_setup_project.sh
# Creates the full project directory scaffold.
# Run once from the project root directory.
# =============================================================================
set -euo pipefail

echo "Creating project directory structure..."

# --- Core pipeline directories ---
mkdir -p fastq
mkdir -p metadata/popmaps
mkdir -p metadata/keeplists

mkdir -p 00_qc/fastqc

mkdir -p 01_stacks/input
mkdir -p 01_stacks/runs
mkdir -p 01_stacks/logs

mkdir -p 02_datasets/populations_runs
mkdir -p 02_datasets/logs

mkdir -p 03_analyses/pca
mkdir -p 03_analyses/admixture
mkdir -p 03_analyses/pixy
mkdir -p 03_analyses/private_alleles
mkdir -p 03_analyses/phylogeny
mkdir -p 03_analyses/kinship
mkdir -p 03_analyses/gadma

mkdir -p 04_plots

mkdir -p scripts

echo ""
echo "Directory structure created:"
find . -maxdepth 3 -type d | sort | sed 's|^\./||' | grep -v '^\.' | awk '{print "  " $0}'
echo ""
echo "Next steps:"
echo "  1. Place demultiplexed FASTQs in: fastq/"
echo "  2. Place Solanum_Metadata.csv in: metadata/"
echo "  3. Run: bash scripts/00_qc.sh"
