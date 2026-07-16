#!/usr/bin/env bash
# =============================================================================
# 03_prefilter_vcf_popgen.sh
# Pre-filters all-SNP VCFs to remove paralogous loci before He/Ho/FIS/FST.
#
# FILTERS AND RATIONALE:
#
#   Call-rate filter (F_MISSING):
#     Within-taxon datasets (A2-A4) use --max-missing 0.80 (>=80% genotyped).
#     The populations runs used -p 1 -r 0.70, so 90% was over-aggressive and
#     removed too many otherwise usable loci. 80% is consistent with the
#     call-rate filter applied in the R script (LOCUS_CALLRATE_MIN=0.80) and
#     avoids redundant over-filtering before the het step does its work.
#     All-taxa dataset (A1) also uses 0.80 for consistency.
#
#   Max observed heterozygosity <= 0.60 (python from GT fields):
#     Primary paralog filter. Loci where >60% of called individuals are
#     heterozygous are almost certainly collapsed paralogs. Computed directly
#     from GT fields to avoid bcftools +fill-tags version dependencies.
#
#   HWE filter (A2/A3/A4 only — NOT applied to A1):
#     Within-taxon datasets: removes loci with HWE p < 0.001. Paralogs cause
#     near-universal heterozygosity and extreme HWE deviation.
#     A1 all-taxa: HWE filter is SKIPPED. Mixing three differentiated species
#     in a single HWE test causes the Wahlund effect — heterozygote deficiency
#     at differentiated loci due to population structure, not paralogs. Applying
#     HWE to the all-taxa dataset would remove genuine FST-informative SNPs.
#     The max-het filter alone is sufficient for A1 paralog removal.
#
# Requirements: bcftools (any recent version), bgzip, tabix, python3
# Run from project root.
# =============================================================================
set -euo pipefail

POPRUN_BASE="02_datasets/populations_runs"
LOGDIR="02_datasets/logs"
mkdir -p "${LOGDIR}"

# Thresholds
MAX_F_MISSING=0.20    # >=80% samples genotyped per locus (F_MISSING <= 0.20)
MAX_SITE_HET=0.60     # max observed heterozygosity per site
HWE_PVAL=0.001        # HWE filter applied to within-taxon datasets only

# Dataset definitions: "tag|apply_hwe"
# apply_hwe=1 -> run HWE filter; apply_hwe=0 -> skip (Wahlund effect)
DATASETS=(
    "A1_all_taxa_noMAF_allSNP|0"
    "A2_Ss_F1F2_noMAF_allSNP|1"
    "A3_Si_wildstock_noMAF_allSNP|1"
    "A4_Sk_F1F2_noMAF_allSNP|1"
    # Individual-level popmap runs for kinship analyses
    # HWE filter skipped: each sample is its own population (n=1),
    # making HWE undefined and p-values meaningless
    "A2_Ss_individual_noMAF_allSNP|0"
)

# Inline Python het filter: reads VCF body from stdin,
# writes CHROM<tab>POS of sites with obs_het <= threshold to outfile
cat > /tmp/het_filter.py << 'PYEOF'
import sys, re
max_het  = float(sys.argv[1])
out_path = sys.argv[2]
het_re   = re.compile(r'^0[/|]1$|^1[/|]0$')
miss_re  = re.compile(r'^\.$|^\./\.$|^\.\|\.$')
passed   = []
for line in sys.stdin:
    if line.startswith('#'):
        continue
    fields = line.rstrip('\n').split('\t')
    chrom, pos = fields[0], fields[1]
    gts    = [f.split(':')[0] for f in fields[9:]]
    called = [g for g in gts if not miss_re.match(g)]
    if not called:
        continue
    n_het   = sum(1 for g in called if het_re.match(g))
    obs_het = n_het / len(called)
    if obs_het <= max_het:
        passed.append(chrom + '\t' + pos + '\n')
with open(out_path, 'w') as fh:
    fh.writelines(passed)
print(f"Sites passing het filter: {len(passed)}", file=sys.stderr)
PYEOF

filter_vcf() {
    local tag="$1"
    local apply_hwe="$2"
    local invcf="${POPRUN_BASE}/${tag}/populations.snps.fixed.vcf.gz"
    local outvcf="${POPRUN_BASE}/${tag}/populations.snps.filtered.vcf.gz"
    local logfile="${LOGDIR}/prefilter_${tag}.log"
    local tmpdir
    tmpdir="$(mktemp -d)"
    trap 'rm -rf "${tmpdir}"' RETURN

    echo "================================================================"
    echo "  TAG       : ${tag}"
    echo "  HWE filter: $([ "${apply_hwe}" -eq 1 ] && echo yes || echo 'no -- skipped (Wahlund effect or individual-level popmap)')"
    echo "  INPUT     : ${invcf}"

    if [[ ! -s "${invcf}" ]]; then
        echo "  SKIP -- input VCF not found"; return
    fi

    bcftools view -h "${invcf}" > /dev/null 2>&1 \
        || { echo "  ERROR: bcftools cannot read ${invcf}"; return; }

    local n_in
    n_in=$(bcftools view -H "${invcf}" | wc -l)
    echo "  Sites before filtering          : ${n_in}"

    # -- Filter 1: call-rate --------------------------------------------------
    local vcf_cr="${tmpdir}/step1_callrate.vcf.gz"
    bcftools view \
        -i "F_MISSING <= ${MAX_F_MISSING}" \
        -Oz -o "${vcf_cr}" "${invcf}"
    tabix -f -p vcf "${vcf_cr}"
    local n_cr
    n_cr=$(bcftools view -H "${vcf_cr}" | wc -l)
    echo "  After call-rate (>=80%)         : ${n_cr}  (removed $(( n_in - n_cr )))"

    # -- Filter 2: max observed heterozygosity --------------------------------
    local het_sites="${tmpdir}/het_pass_sites.txt"
    bcftools view -H "${vcf_cr}" \
        | python3 /tmp/het_filter.py "${MAX_SITE_HET}" "${het_sites}" 2>/dev/null

    local vcf_het="${tmpdir}/step2_het.vcf.gz"
    bcftools view \
        -T "${het_sites}" \
        -Oz -o "${vcf_het}" "${vcf_cr}"
    tabix -f -p vcf "${vcf_het}"
    local n_het
    n_het=$(bcftools view -H "${vcf_het}" | wc -l)
    echo "  After max-het (<=0.60)          : ${n_het}  (removed $(( n_cr - n_het )))"

    # -- Filter 3: HWE (within-taxon datasets only) ---------------------------
    local n_out
    if [[ "${apply_hwe}" -eq 1 ]]; then
        local vcf_hwe_tags="${tmpdir}/step3_hwe_tags.vcf.gz"
        bcftools +fill-tags \
            "${vcf_het}" \
            -Oz -o "${vcf_hwe_tags}" \
            -- -t HWE
        tabix -f -p vcf "${vcf_hwe_tags}"

        bcftools view \
            -i "HWE >= ${HWE_PVAL}" \
            -Oz -o "${outvcf}" "${vcf_hwe_tags}"
        tabix -f -p vcf "${outvcf}"
        n_out=$(bcftools view -H "${outvcf}" | wc -l)
        echo "  After HWE (p>=0.001)            : ${n_out}  (removed $(( n_het - n_out )))"
    else
        # No HWE filter for all-taxa dataset -- copy het-filtered VCF as final output
        cp "${vcf_het}" "${outvcf}"
        cp "${vcf_het}.tbi" "${outvcf}.tbi"
        n_out="${n_het}"
        echo "  HWE filter                      : skipped (Wahlund effect)"
    fi

    # -- Summary --------------------------------------------------------------
    local n_removed=$(( n_in - n_out ))
    local pct_removed
    pct_removed=$(awk "BEGIN {printf \"%.1f\", ${n_removed}/${n_in}*100}")
    echo ""
    echo "  Total removed : ${n_removed} / ${n_in} (${pct_removed}%)"
    echo "  Final sites   : ${n_out}"
    echo "  Output        : ${outvcf}"

    {
        echo "Tag: ${tag}"
        echo "HWE filter applied: $([ "${apply_hwe}" -eq 1 ] && echo yes || echo no)"
        echo "Input sites           : ${n_in}"
        echo "After call-rate (80%) : ${n_cr}  (removed $(( n_in  - n_cr  )))"
        echo "After max-het (0.60)  : ${n_het}  (removed $(( n_cr  - n_het )))"
        if [[ "${apply_hwe}" -eq 1 ]]; then
        echo "After HWE (0.001)     : ${n_out}  (removed $(( n_het - n_out )))"
        else
        echo "After HWE             : skipped"
        fi
        echo "Total removed         : ${n_removed} (${pct_removed}%)"
        echo "Final sites           : ${n_out}"
    } > "${logfile}"
    echo "  Log           : ${logfile}"
    echo ""
}

echo "Thresholds:"
echo "  Call-rate    : F_MISSING <= ${MAX_F_MISSING}  (>=80% samples genotyped)"
echo "  Max obs. het : <= ${MAX_SITE_HET}"
echo "  HWE p-value  : >= ${HWE_PVAL}  (within-taxon datasets only)"
echo ""

for entry in "${DATASETS[@]}"; do
    tag="${entry%%|*}"
    apply_hwe="${entry##*|}"
    filter_vcf "${tag}" "${apply_hwe}"
done

rm -f /tmp/het_filter.py

echo "================================================================"
echo "Pre-filtering complete."
echo ""
echo "Expected final site counts (approximate):"
echo "  A1 all-taxa  : ~6000  (call-rate + max-het only)"
echo "  A2 Ss F1/F2  : ~1000-1500"
echo "  A3 Si        : ~2500-3000"
echo "  A4 Sk        : ~3000-3500"
echo ""
echo "Next: Rscript scripts/03_popgen_diversity_fst.R"
