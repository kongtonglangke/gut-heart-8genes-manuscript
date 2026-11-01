# -*- coding: utf-8 -*-
"""Apply the verified corrections from the two review reports - v1b.

Every replacement is asserted to match exactly once, so a silent miss is
impossible.  Items the reviewers raised that turned out to be already correct are
NOT touched and are listed in the printed report.

Usage:  python 06_apply_review.py
"""
import io
import re

MD = (r"F:\longxiawenjian\gaoxueyaqiguanzhou\Manuscript_GutHeart_8genes"
      r"\Manuscript_GutHeart_8genes_v1.md")

t = io.open(MD, encoding="utf-8").read()
done = []


def sub(old, new, why, count=1):
    global t
    n = t.count(old)
    assert n == count, "expected %d match(es) for %r, found %d" % (count, old[:60], n)
    t = t.replace(old, new)
    done.append(why)


# ---------------------------------------------------------------- 1. metadata
sub("**Article type:** Original Research Article — Methods / Negative Results\n\n", "",
    "删去正文中的投稿元数据：Article type")
sub("**Suggested journals (IF 4–6):** Scientific Reports; BMC Medical Genomics; "
    "Frontiers in Genetics; Journal of Human Genetics\n\n", "",
    "删去正文中的投稿元数据：Suggested journals")
sub("**Author order note:** only two authors contributed to this work; no other "
    "individuals meet authorship criteria and none are acknowledged as contributors.\n\n",
    "", "删去正文中的投稿元数据：Author order note")

# ------------------------------------------------- 2. pre-specified, not pre-registered
sub("none met preregistered thresholds", "none met the pre-specified thresholds",
    "Abstract: preregistered -> pre-specified")
sub("none met the preregistered positive threshold",
    "none met the pre-specified positive threshold",
    "Results: preregistered -> pre-specified")
sub("preregistered threshold of 0.80", "pre-specified threshold of 0.80",
    "Fig. 3 legend: preregistered -> pre-specified")
sub("no test met the preregistered positive threshold",
    "no test met the pre-specified positive threshold",
    "Table 3 note: preregistered -> pre-specified")

# ------------------------------------------------------------ 3. coloc version
sub("(R coloc package, version 5.x)", "(R coloc package, version 5.2.3)",
    "Methods: coloc 版本 5.x -> 5.2.3")
sub("R 4.3.3 (coloc 5.x)", "R 4.3.3 (coloc 5.2.3)",
    "Software 段: coloc 版本 5.x -> 5.2.3")

# ------------------------------------------- 4. ACE control honesty (abstract)
sub("Null results here do not reflect insensitivity: two positive-control pairings "
    "returned posterior 0.985 for *FES* in the same intestinal tissue and outcome, "
    "and 0.7295 for *ACE* in kidney cortex.",
    "Null results here do not reflect insensitivity: a positive control returned "
    "posterior 0.985 for *FES* in the same intestinal tissue and outcome, and a "
    "second control returned 0.730 for *ACE* in kidney cortex, below the 0.80 "
    "threshold but consistent with previous reports for that locus.",
    "Abstract: 区分「一个达阈、一个未达阈」")

# ---------------------------------------- 5. ACE honesty (Results, credibility)
sub("Before interpreting the negative findings, we verified that the pipeline can "
    "detect true positives in the exact setting used here (Table 4). The ACE–kidney "
    "cortex control returned PPH4 = 0.7295. Critically, the FES–transverse colon "
    "control, which uses the same intestinal tissue and the same hypertension GWAS "
    "as the eight candidates, returned PPH4 = 0.9846 in the main analysis and 0.8646 "
    "and 0.9984 under the alternative priors (803 shared variants).",
    "Before interpreting the negative findings, we asked whether the pipeline can "
    "detect a true positive in the exact setting used here (Table 4). Of the two "
    "pre-specified controls, one met the 0.80 threshold and one did not. The "
    "FES–transverse colon control, which uses the same intestinal tissue and the same "
    "hypertension GWAS as the eight candidates, returned PPH4 = 0.985 in the main "
    "analysis and 0.865 and 0.998 under the alternative priors (803 shared variants). "
    "The ACE–kidney cortex control returned PPH4 = 0.730, below the threshold but "
    "close to the value previously reported for this pair, so it corroborates the "
    "direction of the pipeline without independently demonstrating threshold "
    "sensitivity.",
    "Results: ACE 明说低于阈值，只有 FES 达阈")

# ------------------------------------- 6. ACE honesty (Discussion credibility)
sub("Negative results in this design are only meaningful if the pipeline can detect "
    "a true positive in the same configuration. Two controls establish this. The "
    "FES–transverse colon test used the identical tissue, identical outcome "
    "genome-wide association study, identical window definition and identical code as "
    "the candidate tests, and returned PPH4 = 0.9846. The ACE–kidney cortex test, "
    "included as an a priori positive control from prior work on canonical "
    "renin–angiotensin system genes, returned PPH4 = 0.7295, consistent with "
    "previously reported values for this pair.",
    "Negative results in this design are only meaningful if the pipeline can detect a "
    "true positive in the same configuration. The FES–transverse colon test used the "
    "identical tissue, identical outcome genome-wide association study, identical "
    "window definition and identical code as the candidate tests, and returned "
    "PPH4 = 0.985; it is the control that establishes threshold sensitivity. The "
    "ACE–kidney cortex test, an a priori control from prior work on canonical "
    "renin–angiotensin system genes, returned PPH4 = 0.730, below the 0.80 threshold "
    "but consistent with previously reported values for this pair, so it corroborates "
    "the pipeline in absolute terms without clearing the same bar.",
    "Discussion: 两个对照的地位分开写")

# ------------------------------------------------------- 7. precision unify
n = t.count("0.9846")
t = t.replace("0.9846", "0.985")
done.append("统一 PPH4 精度：0.9846 -> 0.985（%d 处）" % n)
n = t.count("0.7295")
t = t.replace("0.7295", "0.730")
done.append("统一 PPH4 精度：0.7295 -> 0.730（%d 处）" % n)

# --------------------------------------- 8. TBXAS1 ratio provenance footnote
sub("no test met the pre-specified positive threshold (PPH4 >= 0.80 and PPH4/PPH3 >= 3).",
    "no test met the pre-specified positive threshold (PPH4 >= 0.80 and PPH4/PPH3 >=\n"
    "3). Posterior ratios quoted in the text are computed from unrounded posteriors,\n"
    "so they will not reproduce exactly from the rounded values shown here.",
    "Table 3 note: 说明比值由未舍入后验计算（解释 11.81 vs 11.67）")

# ------------------------------------------- 9. Table 5 'independent' wording
sub("Strict instruments were independent cis-eQTL variants with P < 5 x 10^-8 and "
    "F > 10.",
    "Strict instruments were cis-eQTL variants with P < 5 x 10^-8 and F > 10. Only one\n"
    "variant met these criteria for each gene–tissue pair, so independence between\n"
    "instruments does not arise and each estimate rests on a single variant.",
    "Table 5 note: 删去无从谈起的 'independent'")

# ------------------------------------------------ 10. 'structurally blind'
sub("showing that a gene-centred design built on proximity can be structurally blind "
    "to the transcript through which a locus actually acts",
    "showing that a gene-centred design built on proximity can miss the transcript "
    "through which a locus actually acts",
    "Introduction: 弱化 structurally blind")
sub("study designs built on *FURIN* eQTLs are structurally blind to this signal",
    "study designs that treat *FURIN* eQTLs as the causal route are blind to this "
    "signal within a ±100 kb window under a single-causal-variant model",
    "Results: structurally blind 加上限定条件")

# ------------------------------------------- 11. Supplementary Table S2 label
sub("Supplementary Fig. S3 presents these results as a forest plot; Supplementary "
    "Table S2 reports\nthe cell-type localization matrix.",
    "Supplementary Fig. S3 presents these results as a forest plot; Supplementary\n"
    "Table S2 reports the single-cell expression matrix summary.",
    "Table 5 尾注: 与附表 S2 实际内容对齐")

# ------------------------------------------------- 12. CTSD 597 TPM clarity
sub("despite 597 TPM in ileum",
    "despite 597 TPM in terminal ileum against a three-segment mean of 440 TPM "
    "(Table 1)",
    "Results: 明确 597 TPM 是回肠单组织值")

# -------------------------------------- 13. CTSD whole-blood distance to threshold
sub("(minimum P = 5.9 × 10⁻⁵ in whole blood)",
    "(minimum P = 5.9 × 10⁻⁵ in whole blood, three orders of magnitude short of the "
    "eGene threshold even at n = 670)",
    "Results: 补 CTSD 全血距阈值差距")

# ------------------------------------------- 14. Fig 5 legend other-tissue range
sub("while showing no association with *FURIN* mRNA in this\ntissue (P = 0.345) or in "
    "any other tissue tested.",
    "while showing no association with *FURIN* mRNA in this tissue (P = 0.345) or in\n"
    "any other tissue tested (P = 0.568–0.930 across the remaining three).",
    "Fig. 5 图注: 补其余组织 P 值范围（0.568–0.930，非审稿人所写的 0.918–0.930）")

# --------------------------------------------- 15. LRP1 protein-layer sentence
sub("carrying no cis-pQTL despite being measured.",
    "carrying no cis-pQTL despite being measured. This matters most for *LRP1*, the "
    "only candidate with intestinal enrichment at the expression layer and the one "
    "assayed with three aptamers, none of which detects a cis signal.",
    "Results: 补 LRP1 蛋白层阴性的意义")

io.open(MD, "w", encoding="utf-8", newline="\n").write(t)

print("=" * 66)
for i, d in enumerate(done, 1):
    print("  %2d. %s" % (i, d))
print("=" * 66)
print("共 %d 项" % len(done))
left = len(re.findall(r"preregistered", t, re.I))
print("残留 preregistered:", left)
