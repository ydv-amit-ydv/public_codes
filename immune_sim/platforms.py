#!/usr/bin/env python3
"""
Vaccine platforms and adjuvants, parameterised from wet-lab literature.

The point of this module is that a vaccine is not just "antigen". How the
antigen gets made, how long it lasts, and which presentation pathway it enters
differ enormously between platforms, and those differences are measurable:

  * A recombinant PROTEIN is a bolus. It decays. It enters the endosome and is
    presented on MHC class II, so it drives CD4 help and antibody well, but it
    reaches MHC class I only by cross-presentation, which is inefficient.
    This is why protein + alum vaccines give good titers and poor CTL.

  * mRNA-LNP is a TEMPLATE. The mRNA is translated inside the cell for a day
    or two, so antigen is made endogenously and enters the MHC class I pathway
    directly - hence the strong CD8 responses seen with mRNA vaccines. The
    ionizable lipid is itself adjuvanting, so no separate adjuvant is needed.

  * circRNA has no free 5'/3' ends, so exonucleases cannot chew it. Expression
    lasts substantially longer than linear mRNA from the same dose.

  * saRNA carries a replicase, so the template amplifies before it decays;
    protective responses at 1/10 to 1/100 of the mRNA dose.

  * Viral vectors express strongly and are intrinsically adjuvanting, but
    anti-vector immunity blunts a homologous boost - the reason ChAdOx and
    Ad26 regimens use heterologous boosting.

NUMBERS AND THEIR STATUS
------------------------
`half_life_h`, `expression_days` and the dose ranges are order-of-magnitude
values taken from the published behaviour of each platform. Everything else
(`mhc1_routing`, `innate`, `th1_bias`, ...) is a RELATIVE scale on [0, 1] or a
multiplier, set so the simulated ranking reproduces documented qualitative
findings - not a measured constant. The benchmark suite in
`test_platforms.py` states exactly which published comparisons are being
reproduced. Treat these as a structured, falsifiable prior, not as data.
"""
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Platform:
    key: str
    label: str
    # --- antigen supply ---
    template: bool          # True = antigen is translated from a template
    half_life_h: float      # template half-life (h); for protein, antigen itself
    expression_days: float  # how long antigen is produced at a useful level
    amplifies: float        # replicon amplification rate (/day); 0 = none
    translation: float      # antigen made per unit template per day
    # --- presentation routing ---
    mhc1_routing: float     # [0,1] efficiency into MHC class I (endogenous)
    mhc2_routing: float     # [0,1] efficiency into MHC class II (exogenous)
    secreted: float         # [0,1] fraction reaching B cells as native antigen
    # --- innate activation ---
    innate: float           # intrinsic danger signal per dose unit
    innate_days: float      # how long that innate signal persists
    self_adjuvanting: bool  # carries its own adjuvant activity
    # --- practicalities that change a schedule ---
    anti_vector: float      # [0,1] loss of boost efficacy from anti-vector immunity
    reactogenicity: float   # [0,1] relative, for dose-limiting discussion
    typical_dose_ug: tuple  # (low, high) human dose range in micrograms
    notes: str = ""


PLATFORMS = {p.key: p for p in [
    Platform(
        "protein", "Recombinant protein subunit",
        template=False, half_life_h=36.0, expression_days=0.0, amplifies=0.0,
        translation=0.0,
        mhc1_routing=0.08,      # cross-presentation only
        mhc2_routing=1.00, secreted=1.00,
        innate=0.02, innate_days=1.0, self_adjuvanting=False,
        anti_vector=0.0, reactogenicity=0.20, typical_dose_ug=(5.0, 50.0),
        notes="Needs an exogenous adjuvant. Strong antibody, weak CTL."),

    Platform(
        "peptide", "Synthetic long peptide",
        template=False, half_life_h=8.0, expression_days=0.0, amplifies=0.0,
        translation=0.0,
        mhc1_routing=0.25, mhc2_routing=0.85, secreted=0.30,
        innate=0.01, innate_days=1.0, self_adjuvanting=False,
        anti_vector=0.0, reactogenicity=0.10, typical_dose_ug=(100.0, 1000.0),
        notes="Short-lived; depends entirely on the adjuvant."),

    Platform(
        "mrna_lnp", "mRNA in lipid nanoparticle",
        template=True, half_life_h=14.0, expression_days=3.0, amplifies=0.0,
        translation=380.0,
        mhc1_routing=0.95, mhc2_routing=0.80, secreted=0.55,
        innate=0.55, innate_days=2.5, self_adjuvanting=True,
        anti_vector=0.0, reactogenicity=0.70, typical_dose_ug=(30.0, 100.0),
        notes="Endogenous expression gives strong CD8. Ionizable lipid "
              "is intrinsically adjuvanting."),

    Platform(
        "circrna_lnp", "circRNA in lipid nanoparticle",
        template=True, half_life_h=48.0, expression_days=8.0, amplifies=0.0,
        translation=300.0,
        mhc1_routing=0.95, mhc2_routing=0.80, secreted=0.55,
        innate=0.40, innate_days=3.0, self_adjuvanting=True,
        anti_vector=0.0, reactogenicity=0.55, typical_dose_ug=(10.0, 50.0),
        notes="No free ends, so exonuclease-resistant: expression lasts far "
              "longer than linear mRNA. Pairs with circrna_opt in this repo."),

    Platform(
        "sarna_lnp", "Self-amplifying RNA in LNP",
        template=True, half_life_h=20.0, expression_days=10.0, amplifies=1.1,
        translation=260.0,
        mhc1_routing=0.95, mhc2_routing=0.80, secreted=0.55,
        innate=0.75, innate_days=4.0, self_adjuvanting=True,
        anti_vector=0.0, reactogenicity=0.80, typical_dose_ug=(0.5, 10.0),
        notes="Replicase amplifies the template: dose-sparing, but dsRNA "
              "intermediates raise innate activation and reactogenicity."),

    Platform(
        "adenovirus", "Adenoviral vector",
        template=True, half_life_h=30.0, expression_days=14.0, amplifies=0.0,
        translation=420.0,
        mhc1_routing=0.95, mhc2_routing=0.75, secreted=0.45,
        innate=0.65, innate_days=5.0, self_adjuvanting=True,
        anti_vector=0.65, reactogenicity=0.55, typical_dose_ug=(1.0, 10.0),
        notes="Strong, durable expression and strong CTL, but anti-vector "
              "immunity blunts a homologous boost - use heterologous."),

    Platform(
        "mva", "Modified vaccinia Ankara",
        template=True, half_life_h=24.0, expression_days=7.0, amplifies=0.0,
        translation=350.0,
        mhc1_routing=0.90, mhc2_routing=0.80, secreted=0.50,
        innate=0.70, innate_days=4.0, self_adjuvanting=True,
        anti_vector=0.50, reactogenicity=0.50, typical_dose_ug=(1.0, 10.0),
        notes="Common heterologous boost after an adenoviral prime."),

    Platform(
        "dna", "Plasmid DNA",
        template=True, half_life_h=72.0, expression_days=14.0, amplifies=0.0,
        translation=45.0,
        mhc1_routing=0.90, mhc2_routing=0.70, secreted=0.40,
        innate=0.30, innate_days=3.0, self_adjuvanting=True,
        anti_vector=0.0, reactogenicity=0.15, typical_dose_ug=(500.0, 4000.0),
        notes="Durable but weakly expressed in humans without "
              "electroporation; cGAS-STING sensing of cytosolic DNA."),

    Platform(
        "inactivated", "Inactivated whole pathogen",
        template=False, half_life_h=48.0, expression_days=0.0, amplifies=0.0,
        translation=0.0,
        mhc1_routing=0.12, mhc2_routing=1.00, secreted=1.00,
        innate=0.25, innate_days=2.0, self_adjuvanting=True,
        anti_vector=0.0, reactogenicity=0.30, typical_dose_ug=(3.0, 30.0),
        notes="Native epitope conformation and some intrinsic PAMPs; "
              "still largely an exogenous, class II route."),

    Platform(
        "live_attenuated", "Live attenuated",
        template=True, half_life_h=36.0, expression_days=21.0, amplifies=0.52,
        translation=300.0,
        mhc1_routing=1.00, mhc2_routing=0.90, secreted=0.70,
        innate=0.80, innate_days=8.0, self_adjuvanting=True,
        anti_vector=0.0, reactogenicity=0.45, typical_dose_ug=(0.1, 5.0),
        notes="Replicates in vivo: the strongest and most durable of all, "
              "with the safety trade-off that implies."),
]}


@dataclass(frozen=True)
class Adjuvant:
    key: str
    label: str
    magnitude: float        # danger signal per dose
    duration_days: float    # how long it keeps signalling
    th1_bias: float         # [0,1]; 0.5 balanced, >0.5 Th1/CTL, <0.5 Th2/antibody
    ctl_boost: float        # multiplier on CD8 priming
    ab_boost: float         # multiplier on B cell activation
    depot: float            # extra antigen retention (slows clearance)
    reactogenicity: float
    notes: str = ""


ADJUVANTS = {a.key: a for a in [
    Adjuvant("none", "No adjuvant", 0.0, 0.0, 0.50, 1.00, 1.00, 0.0, 0.00,
             "Protein alone is poorly immunogenic and can drive tolerance."),
    Adjuvant("alum", "Aluminium salt", 55.0, 4.0, 0.25, 0.85, 1.65, 0.55, 0.20,
             "Classic depot + NLRP3. Th2-biased: strong antibody, poor CTL."),
    Adjuvant("mf59", "MF59 squalene emulsion", 75.0, 3.0, 0.45, 1.05, 1.85, 0.25,
             0.30, "Strong antibody and broadened response; used in flu."),
    Adjuvant("as03", "AS03 squalene + tocopherol", 85.0, 3.0, 0.48, 1.10, 1.95,
             0.25, 0.40, "Potent antibody enhancement, dose-sparing."),
    Adjuvant("as01", "AS01 (MPL + QS-21)", 110.0, 3.5, 0.80, 1.90, 1.60, 0.20,
             0.60, "TLR4 + saponin. Strong Th1 and CD8; Shingrix-type."),
    Adjuvant("cpg1018", "CpG 1018 (TLR9)", 95.0, 3.0, 0.82, 1.70, 1.45, 0.15,
             0.40, "Th1-biased; licensed in HepB-CpG."),
    Adjuvant("polyic", "Poly(I:C) (TLR3)", 100.0, 2.0, 0.85, 1.80, 1.30, 0.10,
             0.55, "Strong type-I interferon, Th1 and cross-presentation."),
    Adjuvant("mpla", "MPLA (TLR4)", 80.0, 2.5, 0.72, 1.45, 1.40, 0.15, 0.35,
             "Detoxified LPS derivative."),
    Adjuvant("lps", "LPS / generic danger signal", 100.0, 2.0, 0.70, 1.40, 1.35,
             0.05, 0.90, "The C-ImmSim default. Too reactogenic for human use."),
    Adjuvant("saponin", "Matrix-M saponin", 95.0, 3.0, 0.68, 1.55, 1.70, 0.20,
             0.35, "Used with Novavax protein vaccine."),
    Adjuvant("intrinsic", "Platform-intrinsic only", 0.0, 0.0, 0.65, 1.15, 1.05,
             0.0, 0.0, "For RNA/vector platforms that adjuvant themselves."),
]}


DEFAULT_ADJUVANT = {
    "protein": "alum", "peptide": "polyic", "inactivated": "alum",
    "mrna_lnp": "intrinsic", "circrna_lnp": "intrinsic", "sarna_lnp": "intrinsic",
    "adenovirus": "intrinsic", "mva": "intrinsic", "dna": "intrinsic",
    "live_attenuated": "intrinsic",
}


def resolve(platform_key, adjuvant_key=None, dose_ug=None):
    """Return (Platform, Adjuvant, dose) with sensible defaults per platform."""
    if platform_key not in PLATFORMS:
        raise KeyError(f"unknown platform {platform_key!r}; "
                       f"choose from {sorted(PLATFORMS)}")
    pf = PLATFORMS[platform_key]
    ak = adjuvant_key or DEFAULT_ADJUVANT.get(platform_key, "none")
    if ak not in ADJUVANTS:
        raise KeyError(f"unknown adjuvant {ak!r}; choose from {sorted(ADJUVANTS)}")
    adj = ADJUVANTS[ak]
    if dose_ug is None:
        lo, hi = pf.typical_dose_ug
        dose_ug = (lo * hi) ** 0.5          # geometric mid of the human range
    return pf, adj, float(dose_ug)


def describe(platform_key, adjuvant_key=None, dose_ug=None):
    pf, adj, dose = resolve(platform_key, adjuvant_key, dose_ug)
    lines = [
        f"{pf.label}  +  {adj.label}",
        f"  dose              {dose:g} ug  (typical {pf.typical_dose_ug[0]:g}-"
        f"{pf.typical_dose_ug[1]:g})",
        f"  antigen supply    " + (
            f"template, t1/2 {pf.half_life_h:g} h, ~{pf.expression_days:g} d of expression"
            + (f", amplifying at {pf.amplifies:g}/d" if pf.amplifies else "")
            if pf.template else f"bolus, t1/2 {pf.half_life_h:g} h"),
        f"  routing           MHC-I {pf.mhc1_routing:.2f}   MHC-II {pf.mhc2_routing:.2f}"
        f"   secreted {pf.secreted:.2f}",
        f"  innate            platform {pf.innate:.2f} for {pf.innate_days:g} d"
        f"   + adjuvant {adj.magnitude:g} for {adj.duration_days:g} d",
        f"  bias              Th1 {adj.th1_bias:.2f}   CTL x{adj.ctl_boost:.2f}"
        f"   antibody x{adj.ab_boost:.2f}",
        f"  reactogenicity    platform {pf.reactogenicity:.2f}, "
        f"adjuvant {adj.reactogenicity:.2f}",
    ]
    if pf.anti_vector:
        lines.append(f"  anti-vector       {pf.anti_vector:.2f} loss on a "
                     f"homologous boost - prefer heterologous")
    if pf.notes:
        lines.append(f"  note              {pf.notes}")
    if adj.notes:
        lines.append(f"  adjuvant note     {adj.notes}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    keys = sys.argv[1:] or list(PLATFORMS)
    for k in keys:
        print(describe(k)); print()
