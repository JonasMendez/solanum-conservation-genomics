# Solanum 2bRAD Population Genomics Pipeline

### Population genomics of Hawaiian *Solanum* (*S. sandwicense from Oahu*, *S. sandwicense from Kauai--syn. S. kavaiensis*, and *S. incompletum*)

This repository contains the analysis scripts and GADMA2 parameter files associated with:

> Publication in Preparation

**Data type:** 2bRAD reduced-representation sequencing  
**Assembly approach:** De novo, Stacks v2.68  
**Study taxa:** *Solanum sandwicense* (Oahu), *S. sandwicense* (Kauai)/ *S. kavaiensis* , *S. incompletum* (Hawaii Island)

---

## Repository Structure

```
solanum-conservation-genomics/
├── README.md
├── metadata/
│   └── Solanum_Metadata.csv          # Master sample sheet
├── scripts/
│   ├── 00_setup_project.sh
│   ├── 00_qc.sh
│   ├── 00_make_metadata_files.py
│   ├── 01_denovo_map.sh
│   ├── 02_populations_runs.sh
│   ├── 02_prep_vcfs.sh
│   ├── 02_fix_phylip.py
│   ├── 02_missingness_report.sh
│   ├── 03_prefilter_vcf_popgen.sh
│   ├── 03_plink_pca.sh
│   ├── 03_plot_pca.R
│   ├── 03_pixy_pi_tajimad.sh
│   ├── 03_popgen_diversity_fst.R
│   ├── 03_kinship_betadosage.R
│   ├── 04_easysfs.sh
│   └── 04_gadma2_rescale_fix.py
└── gadma_params/
    ├── README_gadma.md
    ├── runA_Si_outgroup_primary.txt
    ├── runB_Sk_outgroup_topology_test.txt
    └── runC_Si_outgroup_convergence_mixed.txt
```

---

## Dependencies

All software dependencies are managed via conda. Four environments are used across the pipeline, corresponding to the pipeline steps described below. Full software versions are reported in the manuscript Methods section.

| Environment | Key tools |
|---|---|
| `qc` | FastQC, MultiQC |
| `vcftools` | bcftools, bgzip, tabix, vcftools |
| `analysis` | PLINK v1.9, pixy, IQ-TREE2, R (hierfstat, vcfR, ggplot2, vegan) |
| `gadma2` | GADMA2, moments, easySFS, Python 3.10 |

**Stacks v2.68** is installed locally (no-AVX build) rather than via conda due to hardware compatibility requirements. The local installation path is set at the top of `01_denovo_map.sh` and `02_populations_runs.sh` via:
```bash
export PATH="$HOME/local/stacks-2.68-noavx/bin:$PATH"
```
Update this path to match your local Stacks installation before running those scripts.

---

## Data Availability

Raw sequencing reads (demultiplexed FASTQs) are deposited at NCBI SRA under accession **[SRA accession in preparation]**.

Sample metadata is provided in `metadata/Solanum_Metadata.csv` with the following columns: `sample_id`, `species`, `species_abbrev`, `generation`, `population_code`, `location`, `island`, `source`.

---

## Pipeline Execution Order

All scripts are run from the **project root directory**. Update any hardcoded paths at the top of each script to match your local environment before running.

### Step 0 — Project setup, QC, and metadata

```bash
# Create directory scaffold
bash scripts/00_setup_project.sh

# Place demultiplexed FASTQs in fastq/ and Solanum_Metadata.csv in metadata/

# Quality control (conda env: qc)
bash scripts/00_qc.sh

# Generate all popmaps and keeplists from metadata (conda env: analysis)
python3 scripts/00_make_metadata_files.py
```

Review `00_qc/multiqc_report.html` before proceeding. Flag any samples with anomalous adapter content, low base quality, or unusual duplication rates.

### Step 1 — De novo assembly (Stacks)

```bash
bash scripts/01_denovo_map.sh
```

Runs `denovo_map.pl` on all 37 samples with parameters m=3, M=2, n=1. Input FASTQs must be named `<sample_id>.fastq.gz` and located in `fastq/`. Output goes to `01_stacks/runs/run_m3_M2_n1/`.

### Step 2 — Dataset preparation

```bash
# Run all 17 populations exports (one per analysis dataset; conda env: vcftools)
bash scripts/02_populations_runs.sh

# Fix Stacks VCF headers, bgzip, and tabix-index all VCFs (conda env: vcftools)
bash scripts/02_prep_vcfs.sh

# De-interleave Stacks PHYLIP output and fix truncated sample labels
python3 scripts/02_fix_phylip.py

# Generate per-individual and per-site missingness reports (conda env: vcftools)
bash scripts/02_missingness_report.sh
```

After `02_fix_phylip.py`, pass the four `*_relaxed_full.phy` outputs from `02_datasets/phylip_fixed/` to your IUPAC VCF converter to generate invariant-site VCFs for pixy. bgzip and tabix-index the resulting VCFs before running pixy.

### Step 3 — Analyses

```bash
# Pre-filter VCFs to remove paralogous loci (conda env: vcftools)
bash scripts/03_prefilter_vcf_popgen.sh

# PCA (conda env: vcftools + analysis)
bash scripts/03_plink_pca.sh
Rscript scripts/03_plot_pca.R

# Nucleotide diversity, FST, AMOVA, He/Ho/FIS (conda env: analysis)
Rscript scripts/03_popgen_diversity_fst.R

# Pairwise kinship and breeding recommendations (conda env: analysis)
Rscript scripts/03_kinship_betadosage.R

# Pi, Tajima's D (requires IUPAC VCFs from Step 2; conda env: analysis)
bash scripts/03_pixy_pi_tajimad.sh
```

Note: IQ-TREE2 was run directly on the `A1_all_taxa_phylo` relaxed PHYLIP output from `02_fix_phylip.py` using ModelFinder for model selection followed by tree inference with 1000 bootstrap replicates (`-B 1000`). A wrapper script is not included as IQ-TREE2 was run interactively; see manuscript Methods for the full command.

### Step 4 — Demographic history inference (GADMA2)

```bash
# Step 4a: Generate the joint SFS (conda env: gadma2)
# Review easySFS preview output before running the generate step
bash scripts/04_easysfs.sh

# Step 4b: Run GADMA2 optimization sets (conda env: gadma2)
# See gadma_params/README_gadma.md for run order and rationale
gadma --params gadma_params/runA_Si_outgroup_primary.txt
gadma --params gadma_params/runB_Sk_outgroup_topology_test.txt
# After confirming Si-outgroup topology from AIC comparison (Run A vs Run B):
gadma --params gadma_params/runC_Si_outgroup_convergence_mixed.txt

# Step 4c: Rescale best-fit parameters to demographic units (conda env: gadma2)
python3 scripts/04_gadma2_rescale_fix.py
```

Full details of model parameterization, topology testing, convergence assessment, and parameter rescaling are described in the manuscript Methods and Supplementary Methods.

---

## Key Analytical Decisions

**One assembly, multiple populations runs.** A single `denovo_map.pl` run builds the shared locus catalog across all 37 samples. All downstream filtering is applied at the `populations` stage, which re-reads from the same gstacks catalog. This ensures locus definitions are consistent across datasets and avoids catalog inconsistencies between analyses.

**Paralog pre-filtering.** Initial He/Ho/F~IS~ results showed strongly negative F~IS~ values inconsistent with biology, diagnosed as paralogous loci inflating observed heterozygosity. `03_prefilter_vcf_popgen.sh` applies a maximum observed heterozygosity filter (≤0.60 per site) as the primary paralog filter, followed by call-rate and HWE filters for within-taxon datasets. The HWE filter is deliberately skipped for the all-taxa Dataset A to avoid removing genuine F~ST~-informative SNPs due to the Wahlund effect.

**Random SNP selection.** `--write-random-snp` is used in preference to `--write-single-snp` for PCA and ADMIXTURE datasets. Selecting the first SNP per locus introduces ascertainment bias toward the 5′ end of RAD loci; random selection samples more evenly across loci.

**Individual-level popmaps for PHYLIP outputs.** `populations --phylip-var-all` with a group-level popmap outputs one consensus sequence per population group rather than per individual. All PHYLIP outputs intended for IUPAC VCF conversion (pixy input) or phylogenetic analysis use individual-level popmaps to ensure one sequence per sample.

**GADMA2 model structure.** The `[1,1,1]` structure parameter specifies one size-change epoch per inter-split time period, representing the most parsimonious parameterization for a three-population two-split model. Migration was treated as a free parameter in all periods; consistently near-zero migration values across clean optimization runs are interpreted as supporting allopatric divergence among island populations rather than as optimization failure. See `gadma_params/README_gadma.md` for details on the three optimization sets.

---

## Citation

If you use these scripts or the associated data, please cite:

> Publication in preparation
