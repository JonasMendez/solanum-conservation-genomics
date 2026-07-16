#!/usr/bin/env bash
# =============================================================================
# 02_prep_vcfs.sh
# Post-processes all populations VCF outputs:
#   1. Fixes Stacks VCF headers (adds proper ##contig lines)
#   2. bgzips and tabix-indexes every VCF
#
# pixy-compatible VCFs (with invariant sites) are handled separately:
#   populations outputs --phylip-var-all for the four pixy-target runs,
#   and those PHYLIP files are converted to IUPAC-coded VCFs (with invariant
#   sites) using the external phylip-to-vcf converter script before running
#   pixy. See 02_populations_runs.sh for the PHYLIP output paths.
#
# Requirements (conda env: vcftools):
#   bcftools, bgzip, tabix, htslib
#
# Run from project root.
# =============================================================================
set -euo pipefail

POPBASE="02_datasets/populations_runs"
LOGDIR="02_datasets/logs"

mkdir -p "${LOGDIR}"

# ── helper: fix header + bgzip + tabix one VCF ────────────────────────────
fix_and_index_vcf() {
    local in_vcf="$1"    # raw .vcf from populations
    local out_gz="$2"    # final .vcf.gz

    mkdir -p "$(dirname "${out_gz}")"
    local tmpdir
    tmpdir="$(mktemp -d)"
    trap 'rm -rf "${tmpdir}"' RETURN

    # Extract contigs actually present in the body
    awk '!/^#/ {print $1}' "${in_vcf}" | sort -u \
        | awk '{print "##contig=<ID=" $1 ">"}' > "${tmpdir}/contig_lines.txt"

    # Rebuild header: meta lines (minus any old ##contig=) + new contigs + #CHROM line
    grep '^##' "${in_vcf}" | grep -v '^##contig=' > "${tmpdir}/meta.txt" || true
    grep '^#CHROM' "${in_vcf}" > "${tmpdir}/chrom.txt"

    cat "${tmpdir}/meta.txt" \
        "${tmpdir}/contig_lines.txt" \
        "${tmpdir}/chrom.txt" > "${tmpdir}/header_fixed.txt"

    # Write fixed VCF then compress + index
    {
        cat "${tmpdir}/header_fixed.txt"
        awk '!/^#/' "${in_vcf}"
    } | bgzip -c > "${out_gz}"

    tabix -f -p vcf "${out_gz}"
    echo "    Indexed: ${out_gz}"
}

# ══════════════════════════════════════════════════════════════════════════════
# Fix headers + bgzip + tabix for all populations VCF outputs
# ══════════════════════════════════════════════════════════════════════════════

echo "════════════════════════════════════════════════════════════════"
echo "Fixing VCF headers and indexing all populations outputs"
echo "════════════════════════════════════════════════════════════════"

TAGS=(
    "A1_all_taxa_MAF_randomSNP"
    "A1_all_taxa_noMAF_randomSNP"
    "A1_all_taxa_noMAF_allSNP"
    "A2_Ss_F1F2_noMAF_randomSNP"
    "A2_Ss_F1F2_noMAF_allSNP"
    "A2_Ss_individual_noMAF_allSNP"
    "A3_Si_wildstock_noMAF_randomSNP"
    "A3_Si_wildstock_noMAF_allSNP"
    "A3_Si_individual_noMAF_allSNP"
    "A4_Sk_F1F2_noMAF_randomSNP"
    "A4_Sk_F1F2_noMAF_allSNP"
    "A4_Sk_individual_noMAF_allSNP"
)

for tag in "${TAGS[@]}"; do
    raw_vcf="${POPBASE}/${tag}/populations.snps.vcf"
    out_gz="${POPBASE}/${tag}/populations.snps.fixed.vcf.gz"

    if [[ ! -s "${raw_vcf}" ]]; then
        echo "  SKIP (no VCF): ${tag}"
        continue
    fi

    echo "  Processing: ${tag}"
    fix_and_index_vcf "${raw_vcf}" "${out_gz}"
done

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "VCF prep complete."
echo ""
echo "Fixed + indexed VCFs in each: 02_datasets/populations_runs/<tag>/"
echo ""
echo "Next:"
echo "  1. Run your phylip-to-vcf IUPAC converter on the four PHYLIP files:"
echo "       ${POPBASE}/A1_all_taxa_noMAF_randomSNP/populations.var.phylip"
echo "       ${POPBASE}/A2_Ss_F1F2_noMAF_randomSNP/populations.var.phylip"
echo "       ${POPBASE}/A3_Si_wildstock_noMAF_randomSNP/populations.var.phylip"
echo "       ${POPBASE}/A4_Sk_F1F2_noMAF_randomSNP/populations.var.phylip"
echo "  2. bash scripts/02_missingness_report.sh"
