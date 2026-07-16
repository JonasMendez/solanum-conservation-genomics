#!/usr/bin/env python3
# =============================================================================
# 02_fix_phylip.py
# De-interleaves Stacks PHYLIP output and fixes truncated/colliding sample
# labels for all PHYLIP-producing populations runs.
#
# Stacks outputs interleaved PHYLIP with labels truncated to 10 characters,
# which causes collisions when sample IDs share a prefix. This script:
#   1. Parses the interleaved PHYLIP by POSITION (not by label name),
#      matching rows to sample IDs from the corresponding popmap
#   2. Writes a relaxed full-ID PHYLIP (space-delimited, full sample names)
#      suitable for IQ-TREE and the IUPAC VCF converter
#   3. Writes a 10-char short-ID PHYLIP + bidirectional mapping TSVs
#      for tools that require strict PHYLIP format
#
# Usage:
#   python3 scripts/02_fix_phylip.py
#
# Inputs (from 02_populations_runs.sh outputs):
#   02_datasets/populations_runs/A1_all_taxa_indiv_phylip/populations.all.phylip
#   02_datasets/populations_runs/A2_Ss_indiv_phylip/populations.all.phylip
#   02_datasets/populations_runs/A3_Si_indiv_phylip/populations.all.phylip
#   02_datasets/populations_runs/A4_Sk_indiv_phylip/populations.all.phylip
#   02_datasets/populations_runs/A1_all_taxa_phylo/populations.all.phylip  (IQ-TREE)
#
# Outputs written to: 02_datasets/phylip_fixed/<tag>/
#   <tag>_relaxed_full.phy      relaxed PHYLIP with full sample IDs
#   <tag>_short10.phy           strict 10-char PHYLIP (short IDs)
#   <tag>_short_to_full.tsv     short_id -> full_sample_id mapping
#   <tag>_full_to_short.tsv     full_sample_id -> short_id mapping
#
# Run from project root.
# =============================================================================

import pathlib
import sys

POPMAP_DIR  = pathlib.Path("metadata/popmaps")
POPRUN_DIR  = pathlib.Path("02_datasets/populations_runs")
OUTBASE     = pathlib.Path("02_datasets/phylip_fixed")

# ── datasets to process ───────────────────────────────────────────────────
# (tag, popmap_file)
# The popmap determines the canonical sample order — must match the order
# populations used when generating the PHYLIP.
DATASETS = [
    # pixy-source PHYLIPs (individual-level, one seq per sample)
    ("A1_all_taxa_indiv_phylip", "popmap_all_taxa_individual.tsv"),
    ("A2_Ss_indiv_phylip",       "popmap_Ss_individual.tsv"),
    ("A3_Si_indiv_phylip",       "popmap_Si_individual.tsv"),
    ("A4_Sk_indiv_phylip",       "popmap_Sk_individual.tsv"),
    # IQ-TREE phylo run (also individual-level)
    ("A1_all_taxa_phylo",        "popmap_all_taxa_individual.tsv"),
]

# ── helpers ───────────────────────────────────────────────────────────────

def read_popmap_samples(popmap_path):
    """Return ordered list of sample_ids from column 1 of a popmap TSV."""
    samples = []
    with popmap_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            samples.append(line.split("\t")[0])
    return samples


def parse_interleaved_phylip(phy_path, expected_ntax):
    """
    Parse a Stacks interleaved PHYLIP file.
    Returns list of sequences in taxon order (positional, ignoring labels).
    Raises SystemExit on any format inconsistency.
    """
    lines = phy_path.read_text().splitlines()
    if not lines:
        raise SystemExit(f"Empty PHYLIP: {phy_path}")

    hdr = lines[0].split()
    if len(hdr) != 2:
        raise SystemExit(f"Bad PHYLIP header in {phy_path}: {lines[0]!r}")
    ntax, nchar = int(hdr[0]), int(hdr[1])

    if ntax != expected_ntax:
        raise SystemExit(
            f"{phy_path}: PHYLIP ntax={ntax} but popmap has {expected_ntax} samples.\n"
            f"Make sure the popmap passed to populations matches this PHYLIP."
        )

    seqs = [""] * ntax

    i = 1
    # ── first block: label + sequence chunk ───────────────────────────────
    taxa_read = 0
    while taxa_read < ntax and i < len(lines):
        line = lines[i].rstrip("\n")
        i += 1
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        # first token is the (potentially truncated) label — ignored
        chunk = "".join(parts[1:])
        seqs[taxa_read] += chunk
        taxa_read += 1

    if taxa_read != ntax:
        raise SystemExit(
            f"{phy_path}: first block yielded {taxa_read} taxa, expected {ntax}"
        )

    # ── subsequent blocks: sequence chunks only (no labels) ───────────────
    while i < len(lines) and any(len(s) < nchar for s in seqs):
        if not lines[i].strip():
            i += 1
            continue
        for t in range(ntax):
            if i >= len(lines):
                break
            line = lines[i].strip()
            i += 1
            if not line:
                continue
            seqs[t] += "".join(line.split())

    # ── validate ──────────────────────────────────────────────────────────
    bad = [(idx, len(s)) for idx, s in enumerate(seqs) if len(s) != nchar]
    if bad:
        details = "\n".join(
            f"  taxon index {idx+1}: length {L} (expected {nchar})"
            for idx, L in bad[:20]
        )
        raise SystemExit(f"PHYLIP parse length mismatch in {phy_path}:\n{details}")

    return seqs, nchar


def write_outputs(tag, samples, seqs, nchar, outdir):
    """Write relaxed full-ID PHYLIP, short-ID PHYLIP, and mapping TSVs."""
    outdir.mkdir(parents=True, exist_ok=True)
    ntax = len(samples)

    # Short IDs: S000000001 ... (10 chars, collision-proof)
    short_ids     = [f"S{i:09d}" for i in range(1, ntax + 1)]
    short_to_full = dict(zip(short_ids, samples))
    full_to_short = dict(zip(samples, short_ids))

    # Mapping TSVs
    (outdir / f"{tag}_short_to_full.tsv").write_text(
        "\n".join(f"{s}\t{short_to_full[s]}" for s in short_ids) + "\n"
    )
    (outdir / f"{tag}_full_to_short.tsv").write_text(
        "\n".join(f"{full}\t{sid}" for full, sid in full_to_short.items()) + "\n"
    )

    # Relaxed full-ID PHYLIP (IQ-TREE, IUPAC VCF converter)
    relaxed = outdir / f"{tag}_relaxed_full.phy"
    with relaxed.open("w") as out:
        out.write(f"{ntax} {nchar}\n")
        for full, seq in zip(samples, seqs):
            out.write(f"{full} {seq}\n")

    # Strict 10-char PHYLIP
    short = outdir / f"{tag}_short10.phy"
    with short.open("w") as out:
        out.write(f"{ntax} {nchar}\n")
        for full, seq in zip(samples, seqs):
            out.write(f"{full_to_short[full]} {seq}\n")

    print(f"  → {relaxed}")
    print(f"  → {short}")
    print(f"  → {outdir / f'{tag}_short_to_full.tsv'}")


# ── main ──────────────────────────────────────────────────────────────────

def main():
    errors = []

    for tag, popmap_file in DATASETS:
        print(f"\n{'━'*64}")
        print(f"  TAG   : {tag}")
        print(f"  POPMAP: {popmap_file}")

        popmap_path = POPMAP_DIR / popmap_file
        phy_path    = POPRUN_DIR / tag / "populations.all.phylip"
        outdir      = OUTBASE / tag

        # ── validate inputs ───────────────────────────────────────────────
        if not popmap_path.exists():
            msg = f"  ERROR: popmap not found: {popmap_path}"
            print(msg); errors.append(msg); continue

        if not phy_path.exists():
            msg = f"  ERROR: PHYLIP not found: {phy_path}"
            print(msg); errors.append(msg); continue

        samples = read_popmap_samples(popmap_path)
        print(f"  Samples in popmap : {len(samples)}")

        # ── parse + write ─────────────────────────────────────────────────
        try:
            seqs, nchar = parse_interleaved_phylip(phy_path, len(samples))
        except SystemExit as e:
            msg = f"  ERROR: {e}"
            print(msg); errors.append(msg); continue

        print(f"  Alignment length  : {nchar} bp")
        write_outputs(tag, samples, seqs, nchar, outdir)

    print(f"\n{'━'*64}")
    if errors:
        print(f"Completed with {len(errors)} error(s):")
        for e in errors:
            print(f"  {e}")
        sys.exit(1)
    else:
        print("All PHYLIP files processed successfully.")
        print(f"Fixed PHYLIPs in: {OUTBASE}/")
        print()
        print("Next steps:")
        print("  - Pass *_relaxed_full.phy files for A1-A4 indiv_phylip to your")
        print("    IUPAC VCF converter to generate pixy-compatible VCFs")
        print("  - Pass A1_all_taxa_phylo *_relaxed_full.phy to IQ-TREE, or")
        print("    use the --fasta-loci output with 03b_concat script for")
        print("    supermatrix concatenation before IQ-TREE")
        print("  - Then: bash scripts/02_missingness_report.sh")


if __name__ == "__main__":
    main()
