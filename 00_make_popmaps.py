#!/usr/bin/env python3
# =============================================================================
# 00_make_popmaps.py
# Generates all popmap and keeplist files needed for every populations run
# from the master metadata CSV.
#
# Corrections applied here (do not need to edit the CSV):
#   - Ss_PAKA_10: reassigned to F1 / PAKA (was erroneously coded F2/PAHY)
#   - Si_HAPAMA14_04: generation set to "stock" (was "unknown")
#
# Usage:
#   python3 scripts/00_make_popmaps.py
#
# Outputs written to metadata/popmaps/ and metadata/keeplists/
# =============================================================================

import csv
import pathlib

# ── paths ──────────────────────────────────────────────────────────────────
META_CSV  = pathlib.Path("metadata/Solanum_Metadata.csv")
POPMAP_DIR = pathlib.Path("metadata/popmaps")
KEEP_DIR   = pathlib.Path("metadata/keeplists")
POPMAP_DIR.mkdir(parents=True, exist_ok=True)
KEEP_DIR.mkdir(parents=True, exist_ok=True)

# ── load metadata with corrections ─────────────────────────────────────────
samples = []
with META_CSV.open() as f:
    reader = csv.DictReader(f)
    for row in reader:
        sid        = row["sample_id"].strip()
        species    = row["species"].strip()
        abbrev     = row["species_abbrev"].strip()
        generation = row["generation"].strip()
        pop_code   = row["population_code"].strip()
        location   = row["location"].strip()
        island     = row["island"].strip()
        source     = row["source"].strip()

        # ── corrections ────────────────────────────────────────────────────
        # Ss_PAKA_10 was erroneously coded F2/PAHY — correct to F1/PAKA
        if sid == "Ss_PAKA_10":
            generation = "F1"
            pop_code   = "PAKA"
            location   = "Palikea"

        # Si_HAPAMA14_04 has generation "unknown" — treat as stock
        if sid == "Si_HAPAMA14_04":
            generation = "stock"

        samples.append({
            "sample_id":  sid,
            "species":    species,
            "abbrev":     abbrev,
            "generation": generation,
            "pop_code":   pop_code,
            "location":   location,
            "island":     island,
            "source":     source,
        })

print(f"Loaded {len(samples)} samples from {META_CSV}")

# ── helper ─────────────────────────────────────────────────────────────────
def write_popmap(path, rows):
    """rows: list of (sample_id, population_label)"""
    path = pathlib.Path(path)
    path.write_text("\n".join(f"{s}\t{p}" for s, p in rows) + "\n")
    print(f"  Wrote {path}  ({len(rows)} samples)")

def write_keeplist(path, sample_ids):
    path = pathlib.Path(path)
    path.write_text("\n".join(sample_ids) + "\n")
    print(f"  Wrote {path}  ({len(sample_ids)} samples)")

# ── convenience filters ────────────────────────────────────────────────────
def by_abbrev(*abbrevs):
    return [s for s in samples if s["abbrev"] in abbrevs]

def by_generation(subset, *gens):
    return [s for s in subset if s["generation"] in gens]

# ── sample subsets ─────────────────────────────────────────────────────────
all_ss = by_abbrev("Ss")
all_sk = by_abbrev("Sk")
all_si = by_abbrev("Si")
all_taxa = samples  # all 36 samples

ss_F1 = by_generation(all_ss, "F1")   # PAKA (Palikea) — 10 samples
ss_F2 = by_generation(all_ss, "F2")   # PAHY (Pahole)  —  5 samples

sk_F1 = by_generation(all_sk, "F1")   # NUAA1, MOHA1, LAEA5 — 3 samples
sk_F2 = by_generation(all_sk, "F2")   # KOKA1/3, KOKB2, MAKA326/324 — 5 samples

si_stock = by_generation(all_si, "stock")  # 9 samples (incl. HAPAMA14)
si_wild  = by_generation(all_si, "wild")   # 5 samples

# ── sanity-check counts ────────────────────────────────────────────────────
print("\n── Sample counts after corrections ──────────────────────────────────")
print(f"  Ss total : {len(all_ss):3d}  (F1/PAKA={len(ss_F1)}, F2/PAHY={len(ss_F2)})")
print(f"  Sk total : {len(all_sk):3d}  (F1={len(sk_F1)}, F2={len(sk_F2)})")
print(f"  Si total : {len(all_si):3d}  (stock={len(si_stock)}, wild={len(si_wild)})")
print(f"  ALL total: {len(all_taxa):3d}")
print()

# ══════════════════════════════════════════════════════════════════════════════
# POPMAPS
# ══════════════════════════════════════════════════════════════════════════════

print("── Writing popmaps ───────────────────────────────────────────────────")

# --------------------------------------------------------------------------
# 1. ALL TAXA — species level (Ss / Sk / Si)
#    Used by: Analysis 1 (PCA/ADMIXTURE/pixy/private alleles), GADMA2 SFS
# --------------------------------------------------------------------------
write_popmap(
    POPMAP_DIR / "popmap_all_taxa_species.tsv",
    [(s["sample_id"], s["abbrev"]) for s in all_taxa]
)

# --------------------------------------------------------------------------
# 2. ALL TAXA — individual level (each sample = own population)
#    Used by: phylogeny (--fasta-loci --phylip-var-all), private alleles
# --------------------------------------------------------------------------
write_popmap(
    POPMAP_DIR / "popmap_all_taxa_individual.tsv",
    [(s["sample_id"], s["sample_id"]) for s in all_taxa]
)

# --------------------------------------------------------------------------
# 3. Ss only — F1 vs F2 (PAKA vs PAHY)
#    Used by: Analysis 2 (PCA/ADMIXTURE/He/Ho/FIS/Pi/private alleles)
# --------------------------------------------------------------------------
ss_f1f2 = ss_F1 + ss_F2
write_popmap(
    POPMAP_DIR / "popmap_Ss_F1F2.tsv",
    [(s["sample_id"], s["pop_code"]) for s in ss_f1f2]
)

# --------------------------------------------------------------------------
# 4. Ss only — individual level (for per-individual private alleles/kinship)
# --------------------------------------------------------------------------
write_popmap(
    POPMAP_DIR / "popmap_Ss_individual.tsv",
    [(s["sample_id"], s["sample_id"]) for s in ss_f1f2]
)

# --------------------------------------------------------------------------
# 5. Si only — wild vs stock
#    Used by: Analysis 3 (PCA/ADMIXTURE/pixy/private alleles)
# --------------------------------------------------------------------------
si_all = si_stock + si_wild
write_popmap(
    POPMAP_DIR / "popmap_Si_wildstock.tsv",
    [(s["sample_id"], s["generation"]) for s in si_all]
)

# --------------------------------------------------------------------------
# 6. Si only — individual level (for per-individual private alleles/kinship)
# --------------------------------------------------------------------------
write_popmap(
    POPMAP_DIR / "popmap_Si_individual.tsv",
    [(s["sample_id"], s["sample_id"]) for s in si_all]
)

# --------------------------------------------------------------------------
# 7. Sk only — F1 vs F2
#    Used by: Analysis 4 (PCA/ADMIXTURE/pixy/private alleles)
# --------------------------------------------------------------------------
sk_all = sk_F1 + sk_F2
write_popmap(
    POPMAP_DIR / "popmap_Sk_F1F2.tsv",
    [(s["sample_id"], s["generation"]) for s in sk_all]
)

# --------------------------------------------------------------------------
# 8. Sk only — individual level (for per-individual private alleles/kinship)
# --------------------------------------------------------------------------
write_popmap(
    POPMAP_DIR / "popmap_Sk_individual.tsv",
    [(s["sample_id"], s["sample_id"]) for s in sk_all]
)

# ══════════════════════════════════════════════════════════════════════════════
# KEEPLISTS (sample_id lists for bcftools subsetting)
# ══════════════════════════════════════════════════════════════════════════════

print("\n── Writing keeplists ─────────────────────────────────────────────────")

write_keeplist(KEEP_DIR / "keep_Ss_all.txt",       [s["sample_id"] for s in all_ss])
write_keeplist(KEEP_DIR / "keep_Ss_F1_PAKA.txt",   [s["sample_id"] for s in ss_F1])
write_keeplist(KEEP_DIR / "keep_Ss_F2_PAHY.txt",   [s["sample_id"] for s in ss_F2])
write_keeplist(KEEP_DIR / "keep_Sk_all.txt",        [s["sample_id"] for s in all_sk])
write_keeplist(KEEP_DIR / "keep_Sk_F1.txt",         [s["sample_id"] for s in sk_F1])
write_keeplist(KEEP_DIR / "keep_Sk_F2.txt",         [s["sample_id"] for s in sk_F2])
write_keeplist(KEEP_DIR / "keep_Si_all.txt",        [s["sample_id"] for s in all_si])
write_keeplist(KEEP_DIR / "keep_Si_stock.txt",      [s["sample_id"] for s in si_stock])
write_keeplist(KEEP_DIR / "keep_Si_wild.txt",       [s["sample_id"] for s in si_wild])
write_keeplist(KEEP_DIR / "keep_all_taxa.txt",      [s["sample_id"] for s in all_taxa])

print("\nDone. All popmaps and keeplists written.")
print("Review metadata/popmaps/ and metadata/keeplists/ before running populations.")
