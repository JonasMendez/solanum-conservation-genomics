#!/usr/bin/env python3
# =============================================================================
# 04_gadma2_rescale.py
# Converts GADMA2 output from coalescent-scaled units to demographic units
# (effective population size in individuals, time in generations and years).
#
# GADMA2 outputs parameters scaled relative to a reference population size
# (N_ref) derived from theta (theta = 4 * N_ref * mu * L, where L is the
# number of sites in the SFS). This script:
#   1. Uses hardcoded parameters from GADMA2 Run A Run 10 (best AIC = 2912.12)
#   2. Computes N_ref from theta given mu and L
#   3. Rescales all Ne, time, and migration values
#   4. Applies a range of mu and generation time values to propagate uncertainty
#
# USAGE:
#   python3 scripts/04_gadma2_rescale.py
#
# Run from project root after GADMA2 completes.
# Population order in SFS: [Si, Ss_Oahu, Ss_Kauai]
# =============================================================================

import sys
import csv
from pathlib import Path
from itertools import product

# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================

# Total number of sites (L) used to generate the SFS
L = 462971

# Mutation rate range (substitutions per site per generation)
MU_VALUES = [7e-9, 1e-8]

# Generation time range (years per generation)
GEN_TIME_VALUES = [2, 3, 5]

# Output directory
OUT_DIR = Path("03_analyses/gadma")

# =============================================================================
# ── HARDCODED PARAMETERS — Run A Run 10 (best AIC = 2912.12, theta = 2276.21)
# Population order: [Si, Ss_Oahu, Ss_Kauai]
# Parameter names follow GADMA2 model function variable order:
#   nu_1, nu_2         : ancestral sizes before first split (pre Si/SsSk split)
#   t1                 : time of Si / SsSk split
#   nu11, nu12         : sizes of Si and SsSk lineages in period 1
#   m1_12, m1_21       : migration Si->SsSk and SsSk->Si in period 1
#   nu12_1, nu12_2     : sizes at moment of Ss_Oahu / Ss_Kauai split
#   t2                 : time of Ss_Oahu / Ss_Kauai split
#   nu21, nu22, nu23   : final Ne of Si, Ss_Oahu, Ss_Kauai respectively
#   m2_*               : migration rates in period 2 among all three pops
# =============================================================================

THETA = 2276.21

# Ancestral period sizes (pre-first split)
nu_ancestral = {
    "nu_1_anc_Si":    10.514,
    "nu_2_anc_SsSk":  0.757,
}

# Period 1 sizes (post Si/SsSk split, pre Ss_Oahu/Ss_Kauai split)
nu_period1 = {
    "nu11_Si_p1":     0.139,
    "nu12_SsSk_p1":   0.223,
}

# Sizes at moment of second split
nu_split2 = {
    "nu12_1_SsOahu_at_split": 0.462,
    "nu12_2_SsKauai_at_split": 3.6,
}

# Final period sizes (present-day Ne) — most biologically meaningful
nu_final = {
    "nu21_Si_final":      0.177,
    "nu22_SsOahu_final":  0.020,
    "nu23_SsKauai_final": 0.032,
}

# Combine all Nu parameters in logical order
nu_params = {}
nu_params.update(nu_ancestral)
nu_params.update(nu_period1)
nu_params.update(nu_split2)
nu_params.update(nu_final)

# Divergence times (in units of 2*N_ref generations)
t_params = {
    "t1_Si_SsSk_split":        0.789,
    "t2_SsOahu_SsKauai_split": 0.229,
}

# Migration rates (in units of 1/(2*N_ref))
# Period 1: Si <-> ancestral SsSk
# m1_12 = Si -> SsSk; m1_21 = SsSk -> Si
# Period 2: all pairwise among Si, Ss_Oahu, Ss_Kauai
# m2_12 = Si->SsOahu; m2_13 = Si->SsKauai
# m2_21 = SsOahu->Si; m2_23 = SsOahu->SsKauai
# m2_31 = SsKauai->Si; m2_32 = SsKauai->SsOahu
m_params = {
    "m1_12_Si_to_SsSk":        2.401,
    "m1_21_SsSk_to_Si":        3.58e-5,
    "m2_12_Si_to_SsOahu":      0.286,
    "m2_13_Si_to_SsKauai":     0.079,
    "m2_21_SsOahu_to_Si":      0.134,
    "m2_23_SsOahu_to_SsKauai": 0.611,
    "m2_31_SsKauai_to_Si":     0.045,
    "m2_32_SsKauai_to_SsOahu": 0.280,
}

# =============================================================================
# ── RESCALING ─────────────────────────────────────────────────────────────────
# =============================================================================

OUT_DIR.mkdir(parents=True, exist_ok=True)
out_path = OUT_DIR / "gadma2_rescaled_parameters.csv"

# Build CSV header
header = ["mu", "gen_time_yr", "N_ref"]
for k in nu_params:
    header.append(f"{k}_Ne")
for k in t_params:
    header.append(f"{k}_gen")
    header.append(f"{k}_yr")
for k in m_params:
    header.append(f"{k}_migrants_per_gen")

rows = []

print("=" * 70)
print("GADMA2 PARAMETER RESCALING — Run A Run 10 (best AIC = 2912.12)")
print("=" * 70)
print(f"  Theta            : {THETA}")
print(f"  L (total sites)  : {L:,}")
print(f"  mu values        : {MU_VALUES}")
print(f"  gen time values  : {GEN_TIME_VALUES} years")
print()

for mu, gen_time in product(MU_VALUES, GEN_TIME_VALUES):

    # N_ref from theta: theta = 4 * N_ref * mu * L
    N_ref = THETA / (4 * mu * L)

    print(f"mu={mu:.2e}, gen_time={gen_time}yr  ->  N_ref = {N_ref:,.0f} individuals")
    print("-" * 60)

    row = {"mu": mu, "gen_time_yr": gen_time, "N_ref": round(N_ref)}

    # Rescale Ne: nu_i (fraction of N_ref) -> individuals
    print("  Population sizes:")
    for k, nu_val in nu_params.items():
        ne = nu_val * N_ref
        row[f"{k}_Ne"] = round(ne)
        print(f"    {k:<35} nu={nu_val:.4f}  ->  {ne:>10,.0f} ind")

    # Rescale time: T (units of 2*N_ref generations) -> generations and years
    print("  Divergence times:")
    for k, t_val in t_params.items():
        t_gen = t_val * 2 * N_ref
        t_yr = t_gen * gen_time
        row[f"{k}_gen"] = round(t_gen)
        row[f"{k}_yr"] = round(t_yr)
        print(f"    {k:<35} T={t_val:.4f}  ->  {t_gen:>10,.0f} gen  =  {t_yr:>10,.0f} yr")

    # Rescale migration: m (units of 1/(2*N_ref)) -> migrants per generation
    # Effective migrants per generation = m_scaled * N_ref
    print("  Migration (effective migrants per generation):")
    for k, m_val in m_params.items():
        m_eff = m_val * N_ref
        row[f"{k}_migrants_per_gen"] = round(m_eff, 4)
        print(f"    {k:<35} m={m_val:.5f}  ->  {m_eff:>10.2f} migrants/gen")

    rows.append(row)
    print()

# Write CSV
if rows:
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Rescaled parameters written to: {out_path}")
else:
    print("ERROR: No rows written.")
    sys.exit(1)

print()
print("=" * 70)
print("INTERPRETATION NOTES")
print("=" * 70)
print("  N_ref       : ancestral reference Ne used for coalescent scaling")
print("  nu=1.0      : Ne equals N_ref (no size change)")
print("  nu<1.0      : bottleneck / size reduction relative to N_ref")
print("  nu>1.0      : expansion relative to N_ref")
print("  T_gen       : time in generations (= T_coal * 2 * N_ref)")
print("  T_yr        : time in years (= T_gen * generation_time)")
print("  migrants/gen: effective number of migrants per generation")
print("                values near 0 support allopatric divergence")
print()
print("Population order in SFS and model: [Si, Ss_Oahu, Ss_Kauai]")
print("Final-period Ne (nu21, nu22, nu23) = present-day effective sizes")
print("Report ranges across mu and gen_time to reflect parameter uncertainty")
