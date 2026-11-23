# -*- coding: utf-8 -*-
"""Apply the verified corrections from review round 3 (v1d).

Each replacement is asserted, so a silent miss cannot pass.  Items the reviewers
raised that the manuscript already satisfies are left untouched.
"""
import io

MD = (r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Manuscript_GutHeart_8genes"
      r"\Manuscript_GutHeart_8genes_v1.md")
t = io.open(MD, encoding="utf-8").read()
done = []


def sub(old, new, why, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d got %d for %r" % (n, c, old[:60])
    t = t.replace(old, new)
    done.append(why)


# ------------------------------------------------- 1. Table 5 note duplication
sub("""Strict instruments were cis-eQTL variants with P < 5 x 10^-8 and F > 10. Only one
variant met these criteria for each gene–tissue pair, so independence between
instruments does not arise and each estimate rests on a single variant.
Estimates are odds ratios per one standard deviation higher tissue expression. Only one variant
met these criteria for each gene–tissue pair, so every estimate is a single-instrument Wald
ratio. This avoids the variance inflation that a mismatched linkage disequilibrium reference
panel would introduce in a multi-instrument model, but it forfeits both statistical power and
any test for heterogeneity or directional pleiotropy across instruments. Five of the eight genes
had no instrument meeting these criteria.""",
    """Strict instruments were cis-eQTL variants with P < 5 x 10^-8 and F > 10. Only one
variant met these criteria for each gene–tissue pair, so every estimate is a
single-instrument Wald ratio and independence between instruments does not arise. Estimates
are odds ratios per one standard deviation higher tissue expression. A single-instrument
design avoids the variance inflation that a mismatched linkage disequilibrium reference panel
would introduce in a multi-instrument model, but it forfeits both statistical power and any
test for heterogeneity or directional pleiotropy across instruments. Five of the eight genes
had no instrument meeting these criteria.""",
    "Table 5 表注：删去重复的「只有一个工具」句")

# ------------------------------------------- 2. Methods: split control roles
sub("we pre-specified two positive controls processed with identical code: ACE eQTL in "
    "kidney cortex versus hypertension (expected PPH4 ≈ 0.73 from prior work on this locus), "
    "and FES eQTL in transverse colon versus hypertension, which uses exactly the same tissue "
    "and outcome as the eight candidates.",
    "we pre-specified two controls processed with identical code, chosen for different "
    "purposes. The FES eQTL in transverse colon versus hypertension is the threshold-calibrating "
    "control, because it uses exactly the same tissue and outcome as the eight candidates. The "
    "ACE eQTL in kidney cortex versus hypertension is a literature-replicating control, included "
    "because prior work on this locus reported PPH4 near 0.73; it was not expected to cross the "
    "0.80 threshold and therefore tests reproducibility rather than threshold sensitivity.",
    "Methods：两个对照的角色分开（阈值校准 vs 文献复现）")

# ------------------------------------------------- 3. Table 4 two FES P values
sub("The lead\nvariant rs6227 is simultaneously a *FES* cis-eQTL in transverse colon and the top "
    "hypertension\nvariant in the window.",
    "The window-minimum *FES* eQTL P (5.21 × 10^-13) belongs to a variant other than rs6227; "
    "rs6227 itself reaches P = 4.2 × 10^-9 for *FES* in the same tissue. The lead variant is "
    "simultaneously a *FES* cis-eQTL in transverse colon and the top hypertension variant in "
    "the window.",
    "Table 4 表注：说明 5.21e-13 与 rs6227 的 4.2e-9 是两个不同变异")

# ------------------------------------- 4. Limitations: merge the LD limitation
sub("Second, no linkage disequilibrium reference panel was available, so Mendelian "
    "randomization relied on single-instrument Wald ratios without heterogeneity or pleiotropy "
    "diagnostics. Third, power was limited in kidney cortex and left ventricle, and their null "
    "results cannot be read biologically. Fourth, human bulk and single-cell atlases do not "
    "capture cell-type-specific eQTL or conditional\nregulation, and the absence of an LD "
    "reference panel meant that neither implicated locus could\nbe followed with statistical "
    "fine-mapping or formal conditional analysis, so our regulatory\nlayer is conservative by "
    "design. Fifth, all analyses were undertaken by a single analyst with no blinded duplicate "
    "screening, which we declare explicitly.",
    "Second, this analysis used no linkage disequilibrium reference panel, so Mendelian "
    "randomization relied on single-instrument Wald ratios without heterogeneity or pleiotropy "
    "diagnostics, and neither implicated locus could be followed with statistical fine-mapping "
    "or formal conditional analysis; the regulatory layer is conservative by design. Third, "
    "power was limited in kidney cortex and left ventricle, and their null results cannot be "
    "read biologically. Fourth, human bulk and single-cell atlases do not capture "
    "cell-type-specific eQTL or conditional regulation. Fifth, all analyses were undertaken by "
    "a single analyst with no blinded duplicate screening, which we declare explicitly.",
    "Limitations：合并第 2、4 条（同一 LD 面板根因），六条变五条")
sub("Six limitations bound the interpretation.", "Five limitations bound the interpretation.",
    "Limitations：计数改为五")

# ---------------------------------------- 5. Abstract wording: control gene FES
sub("a control in the same tissue and outcome returned 0.985 for *FES*",
    "a control gene in the same tissue and outcome (*FES*) returned 0.985",
    "摘要：点明 FES 是基因名")

# ------------------------- 6. Abstract: restrict to the array actually checked
sub("whereas *MTHFR* and *FES* are absent from plasma proteomic panels",
    "whereas *MTHFR* and *FES* are absent from the SomaScan v4 panel",
    "摘要：把「血浆蛋白组面板」收窄为实测的 SomaScan v4")

# --------------------------------------------- 7. CTSD: give the fold change
sub("three orders of magnitude short of the eGene threshold even at n = 670",
    "about 1,000-fold above the eGene threshold even at n = 670",
    "Results：CTSD 全血改写成约 1,000 倍（更稳）")

# ------------------------------- 8. Supplementary S4 row wording (Table 4 use)
sub("| Pleiotropy and sensitivity analyses | Table 4 (three priors); Limitations "
    "(single-instrument design precludes heterogeneity tests) |",
    "| Pleiotropy and sensitivity analyses | Table 4 (prior sensitivity); pleiotropy "
    "diagnostics precluded by the single-instrument design (Limitations) |",
    "附表 S4：多效性一行措辞纠正（Table 4 是先验敏感性，不是多效性检验）")

# ------------------------- 9. Data availability: no prospective registration
sub("deCODE plasma pQTL supplementary tables: https://www.decode.com/summarydata/. No new data "
    "were generated.",
    "deCODE plasma pQTL supplementary tables: https://www.decode.com/summarydata/. No new data "
    "were generated. All analytic thresholds were pre-specified by the authors before analysis; "
    "the study was not prospectively registered.",
    "数据可用性：明说阈值由作者事前设定、未前瞻性注册")

# ---------------- 10. Introduction: independent functional literature for FURIN
sub("in which two of them had received endothelial functional validation [1].",
    "in which two of them had received endothelial functional validation [1]. Independent work "
    "corroborates an endothelial role for one of the two: *FURIN* expression in vascular "
    "endothelial cells is modulated by a coronary artery disease-associated variant [13], and "
    "*FURIN* cleaves thrombospondin-1 to drive endothelial cell death and efferocytosis [14]. An "
    "endothelial function for the candidate is therefore well supported even where its "
    "cis-regulation is not.",
    "Introduction：补两篇独立的 FURIN 内皮功能文献（待重排编号）")

io.open(MD, "w", encoding="utf-8", newline="\n").write(t)
for i, d in enumerate(done, 1):
    print("  %2d. %s" % (i, d))
print("共 %d 项" % len(done))
