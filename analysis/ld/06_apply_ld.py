# Step 6 - apply the linkage-disequilibrium revision to the manuscript.
# Every replacement is asserted so a silent no-op cannot happen.

import io
import os
import re

P = (r"F:\longxiawenjian\gaoxueyaqiguanzhou\Manuscript_GutHeart_8genes"
     r"\Manuscript_GutHeart_8genes_v1.md")

t = io.open(P, encoding="utf-8").read()
done = []


def sub(old, new, why):
    global t
    n = t.count(old)
    assert n == 1, "expected 1 occurrence, found %d for: %r" % (n, old[:70])
    t = t.replace(old, new)
    done.append(why)


# ---------------------------------------------------------------- Methods
METHODS = """## Linkage disequilibrium reference panel, fine-mapping and conditional analysis

To remove the single-causal-variant assumption from the two loci that could be followed
up, a reference panel was built from the 1000 Genomes high-coverage phased panel in GRCh38
coordinates [15], restricted to its 633 European-ancestry samples. Alt-allele dosages were
retrieved with tabix for the same ±100 kb windows used throughout, and the pairwise
correlation matrix was computed directly from those dosages after excluding monomorphic
sites. Only variants present in the eQTL, outcome and reference sources jointly were
retained, giving 884 variants at 1p36 and 710 at 15q26.1. Alleles were compared variant by
variant and required no flips.

Fine-mapping used susieR [16] on the standardised z scores with this reference matrix, ten
single-effect components and residual variance fixed, with credible sets defined at 95%
coverage and a minimum absolute within-set correlation of 0.5. Colocalization was then
re-estimated with coloc.susie, which tests every pair of credible sets and therefore
returns a posterior for each pair instead of assuming one causal variant per region.

Conditional association was computed from the summary statistics and the reference matrix
under the joint model, in which the z score at variant i given a conditioning variant j is
(z_i − r_ij z_j) / sqrt(1 − r_ij²). The same expression was evaluated in both directions,
so that each trait was tested for survival of its own signal once the other lead had been
accounted for. Under one shared causal variant the two traits show the same conditional
pattern; when each trait retains a signal the other loses, the region contains two
distinct causal variants.

## cis-Mendelian randomization"""
sub("## cis-Mendelian randomization", METHODS,
    "Methods: new section on the LD panel, fine-mapping and conditional analysis")

# ------------------------------------------------------- Results: MTHFR
MTHFR_PAR = (
    "Because H3 posteriors can be inflated by combinatorial pairing, and because coloc.abf "
    "does not model linkage disequilibrium, we rebuilt the locus with a reference panel and "
    "re-examined it three ways. Fine-mapping returned a single credible set for expression, "
    "containing rs3818762 and rs1023252 (maximum posterior inclusion probability 0.664), so "
    "the *MTHFR* eQTL carries one signal and not two; hypertension, by contrast, resolved "
    "into six credible sets, none of which contained rs3818762. Conditional analysis then "
    "separated the two traits asymmetrically. The expression signal at rs3818762 survived "
    "conditioning on the hypertension lead (z = +7.01, P = 2.4 × 10⁻¹²), whereas the "
    "apparent expression signal at the hypertension lead did not (z = +1.16, P = 0.25); in "
    "the other direction the hypertension signal at rs12567136 survived conditioning on the "
    "expression lead almost undiminished (z = −12.20, P = 3.0 × 10⁻³⁴), whereas the "
    "hypertension signal at rs3818762 did not (z = −2.21, P = 0.027). Each trait therefore "
    "carries its own causal variant, correlated at r = 0.358 but not shared. Colocalization "
    "on credible sets agreed: all six credible-set pairs returned PPH3 between 0.856 and "
    "1.000 and PPH4 at or below 0.002 (Supplementary Table S5), so the negative rests on two "
    "distinct causal variants rather than on an artefact of a single-causal-variant model. "
    "The per-variant Wald ratios are consistent with that architecture. Their median drifted "
    "with distance from the hypertension lead: −0.139 (n = 147), −0.141 (n = 169) and "
    "−0.063 (n = 61) for variants within 20 kb, within 20–50 kb and beyond 50 kb, and the "
    "same monotonic drift was present when the analysis was restricted to transverse colon "
    "alone (−0.219, −0.190 and −0.135 for 19, 44 and 25 variants). One caution follows from "
    "the reference panel. rs12404124 and rs198397, 9.5 kb apart, carry cis-eQTL effects of "
    "opposite sign (β = −0.217 and +0.214), but they are in almost complete negative linkage "
    "disequilibrium (r = −0.992, r² = 0.984), and a single causal variant predicts exactly "
    "that reversal, the expected z at rs198397 being +6.66 against an observed +6.62. Sign "
    "reversal between two markers is therefore not by itself evidence of two causal "
    "variants, and we did not rely on it; the conclusion rests on the conditional and "
    "fine-mapping analyses above. Consistently with linkage contamination rather than a "
    "genuine effect, the *MTHFR* cis-MR estimate pointed in opposite directions across "
    "outcomes: apparently protective for hypertension (OR 0.869, 95% CI 0.835–0.903; "
    "P = 1.9 × 10⁻¹²) but risk-increasing for heart failure (OR 1.096, 95% CI 1.028–1.167; "
    "P = 4.7 × 10⁻³)."
)

lines = t.split("\n")
hit = [i for i, ln in enumerate(lines)
       if ln.startswith("Because H3 posteriors can be inflated by combinatorial pairing")]
assert len(hit) == 1, "MTHFR paragraph not found uniquely: %s" % hit
lines[hit[0]] = MTHFR_PAR
t = "\n".join(lines)
done.append("Results: MTHFR paragraph rewritten around fine-mapping and conditional analysis")

# --------------------------------------------------- Results: FES control
FES_OLD = ("The causal transcript at this blood pressure locus is therefore FES, and study "
           "designs that treat *FURIN* eQTLs as the causal route are blind to this signal "
           "within a ±100 kb window under a single-causal-variant model.")
FES_NEW = ("The causal transcript at this blood pressure locus is therefore FES, and study "
           "designs that treat *FURIN* eQTLs as the causal route are blind to this signal "
           "within a ±100 kb window under a single-causal-variant model. Fine-mapping on the "
           "reference panel reached the same assignment without that assumption: expression "
           "resolved into one credible set of eight variants and hypertension into a "
           "credible set of eight variants containing rs6227, and that pair returned "
           "PPH4 = 0.985 against PPH3 = 0.015, the strongest colocalization in the study. "
           "Conditional analysis confirmed a single shared signal rather than two: "
           "conditioning either lead on the other abolished both the *FES* eQTL association "
           "(residual |z| ≤ 0.80, P ≥ 0.42) and the hypertension association (residual "
           "|z| ≤ 1.11, P ≥ 0.27), the pattern the positive control was designed to "
           "produce.")
sub(FES_OLD, FES_NEW, "Results: FES control now reported with coloc.susie and conditioning")

# ------------------------------------------------------------ Limitations
LIM_OLD = ("Second, this analysis used no linkage disequilibrium reference panel, so "
           "Mendelian randomization relied on single-instrument Wald ratios without "
           "heterogeneity or pleiotropy diagnostics, and neither implicated locus could be "
           "followed with statistical fine-mapping or formal conditional analysis; the "
           "regulatory layer is conservative by design.")
LIM_NEW = ("Second, Mendelian randomization relied on single-instrument Wald ratios without "
           "heterogeneity or pleiotropy diagnostics, and the reference panel covered only "
           "the two loci followed up in detail, so fine-mapping and conditional analysis "
           "were possible for *MTHFR* and *FES* but were not extended across all 240 tests; "
           "the regulatory layer remains conservative by design.")
sub(LIM_OLD, LIM_NEW, "Limitations: LD-panel limitation replaced by the residual one")

# --------------------------------------------------------------- Abstract
AB_OLD = ("all four implicated *MTHFR*, whose expression and hypertension lead variants lay "
          "32.7 kb apart in different genes, with distance-dependent Wald ratio drift "
          "indicating two distinct causal variants.")
AB_NEW = ("all four implicated *MTHFR*, whose expression and hypertension lead variants lay "
          "32.7 kb apart in different genes and which fine-mapping against a 1000 Genomes "
          "reference panel resolved into two trait-specific causal variants (posterior for "
          "a shared variant ≤ 0.002).")
sub(AB_OLD, AB_NEW, "Abstract: fine-mapping result stated")

# --------------------------------------------------------- Discussion text
DIS_OLD = ("At 1p36, the *MTHFR* eQTL lead and the hypertension lead lie 32.7 kb apart in "
           "different nearest genes, and the per-variant Wald ratio drifts with distance and "
           "reverses sign within 10 kb, a signature that coincident MR with a single causal "
           "variant cannot produce.")
DIS_NEW = ("At 1p36, the *MTHFR* eQTL lead and the hypertension lead lie 32.7 kb apart in "
           "different nearest genes, and conditioning each trait on the other lead leaves "
           "exactly one signal standing, a signature that coincident MR with a single causal "
           "variant cannot produce.")
sub(DIS_OLD, DIS_NEW, "Discussion: sign-reversal argument replaced by the conditional result")

# ------------------------------------------------------- Figure 6b caption
FIG_OLD = ("A single shared causal\nvariant would give a horizontal band of constant ratio; "
           "instead the median ratio drifts with\ndistance from the blood pressure lead "
           "(−0.139, −0.141 and −0.063 for variants within 20 kb,\n20–50 kb and beyond "
           "50 kb) while the sign of individual ratios is not preserved: rs12404124 "
           "and\nrs198397, 9.5 kb apart, carry cis-eQTL effects of opposite sign at the same "
           "outcome signal. Both\nfeatures indicate at least two causal variants in partial "
           "linkage disequilibrium.")
FIG_NEW = ("A single shared causal\nvariant would give a horizontal band of constant ratio; "
           "instead the median ratio drifts with\ndistance from the blood pressure lead "
           "(−0.139, −0.141 and −0.063 for variants within 20 kb,\n20–50 kb and beyond "
           "50 kb), as expected when each trait carries its own causal variant. The\nsign of "
           "individual ratios is not conserved, but that feature is not diagnostic on its "
           "own:\nrs12404124 and rs198397 are in almost complete negative linkage "
           "disequilibrium (r = −0.992),\nso opposite signs at the two markers are what a "
           "single causal variant predicts.")
sub(FIG_OLD, FIG_NEW, "Figure 6b caption: sign reversal no longer presented as evidence")

# ------------------------------------- Supplementary table list and table
S4_ROW = ("| Supplementary Table S4 | STROBE-MR checklist, with the location of each item in "
          "this manuscript |")
S5_ROW = S4_ROW + ("\n| Supplementary Table S5 | Reference-panel fine-mapping, colocalization "
                   "on credible sets and conditional analysis for *MTHFR* and *FES* |")
sub(S4_ROW, S5_ROW, "Supplementary table list gains S5")

S5_BLOCK = """
## Supplementary Table S5 Reference-panel fine-mapping, colocalization and conditional analysis

| Locus | Trait | Credible sets | Lead variant | Max PIP |
|---|---|---|---|---|
| 1p36 (*MTHFR*) | Expression, transverse colon | 1 | rs1023252 | 0.664 |
| 1p36 (*MTHFR*) | Hypertension | 6 | rs12567136 (set 4) | 0.220 |
| 15q26.1 (*FES*) | Expression, transverse colon | 1 | rs4932373 | 0.252 |
| 15q26.1 (*FES*) | Hypertension | 3 | rs6227 (set 2) | 0.247 |

| Comparison | coloc.abf PPH3 | coloc.abf PPH4 | coloc.susie PPH3 | coloc.susie PPH4 |
|---|---|---|---|---|
| *MTHFR* expression vs hypertension | 1.000 | 6.5 x 10^-8 | 0.856-1.000 | <= 0.002 |
| *FES* expression vs hypertension | 0.008 | 0.992 | 0.015 | 0.985 |

| Conditional test | z | P |
|---|---|---|
| *MTHFR* expression at rs3818762 given rs12567136 | +7.01 | 2.4 x 10^-12 |
| *MTHFR* expression at rs12567136 given rs3818762 | +1.16 | 0.25 |
| Hypertension at rs12567136 given rs3818762 | -12.20 | 3.0 x 10^-34 |
| Hypertension at rs3818762 given rs12567136 | -2.21 | 0.027 |
| *FES* expression at rs6227 given rs1894400 | -0.08 | 0.94 |
| *FES* expression at rs1894400 given rs6227 | -0.80 | 0.42 |
| Hypertension at rs6227 given rs1894400 | +1.11 | 0.27 |
| Hypertension at rs1894400 given rs6227 | +0.65 | 0.52 |

Note. Reference panel: 1000 Genomes high-coverage phased panel [15], 633 European-ancestry
samples; 884 variants at 1p36 and 710 at 15q26.1 were shared by the eQTL, outcome and
reference sources. Correlation between the two leads was r = 0.358 at 1p36 and r = 0.990 at
15q26.1. Colocalization priors as in the main analysis (p1 = p2 = 1 x 10^-4,
p12 = 1 x 10^-5). The *MTHFR* range spans all six credible-set pairs.
"""

REF_ANCHOR = "\n# References"
assert t.count(REF_ANCHOR) == 1
t = t.replace(REF_ANCHOR, "\n" + S5_BLOCK + "\n# References")
done.append("new Supplementary Table S5 inserted before References")

# ---------------------------------------------------------- New references
lines = t.split("\n")
last_ref = [i for i, ln in enumerate(lines) if re.match(r"^14\. ", ln)]
assert len(last_ref) == 1, "reference 14 not found uniquely: %s" % last_ref
INS = ["15. 1000 Genomes Project Consortium. A global reference for human genetic variation. "
       "*Nature*. 2015;526(7571):68–74. doi:10.1038/nature15393",
       "16. Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable "
       "selection in regression, with application to genetic fine mapping. *J R Stat Soc "
       "Series B*. 2020;82(5):1273–1300. doi:10.1111/rssb.12388"]
lines[last_ref[0] + 1:last_ref[0] + 1] = INS
t = "\n".join(lines)
done.append("two method references appended as [15] and [16]")

io.open(P, "w", encoding="utf-8", newline="\n").write(t)
for i, d in enumerate(done, 1):
    print("  %d. %s" % (i, d))
print("written:", P, os.path.getsize(P), "bytes")
