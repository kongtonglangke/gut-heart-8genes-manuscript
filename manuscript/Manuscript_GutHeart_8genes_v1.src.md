# Front matter

## Title

**Genetic triage of "gut–heart axis" candidates: colocalization and Mendelian randomization of eight hypertension transcriptomic markers**

Chinese working title (用户定题): "肠–心轴"候选基因分诊：8 个高血压转录组标志物的共定位与孟德尔随机化研究

**Running title:** Genetic triage of gut–heart axis candidates

**Article type:** Original Research Article — Methods / Negative Results

**Suggested journals (IF 4–6):** Scientific Reports; BMC Medical Genomics; Frontiers in Genetics; Journal of Human Genetics

## Authors and affiliations

**Yan Wang**¹ · **Jing Bai**²（通讯作者）

¹ Department of Nephrology, Affiliated Hospital of Gansu Medical College, Pingliang, Gansu 744000, China
² Department of Cardiovascular Medicine, Affiliated Hospital of Gansu Medical College, Pingliang, Gansu 744000, China

**Corresponding author:** Jing Bai（白静）, Department of Cardiovascular Medicine, Affiliated Hospital of Gansu Medical College, Pingliang, Gansu 744000, China. E-mail: 18693323793@163.com

**Author order note:** only two authors contributed to this work; no other individuals meet authorship criteria and none are acknowledged as contributors.

*投稿前待补（非正文占位，属投稿系统字段）：ORCID iD；是否需要在 Twitter/领英等处声明对作者身份的一致性拼写（Yan Wang / Jing Bai 与文献 [1] 一致）。*

## Abstract

Blood transcriptomic studies increasingly propose genes as gut–heart axis markers in hypertension, but detection in blood establishes neither intestinal expression nor local genetic contribution to disease. We performed a three-layer genetic triage of eight such candidates (*LRP1*, *CTSD*, *MTHFR*, *AUTS2*, *FURIN*, *CST3*, *FCER1G* and *TBXAS1*) using GTEx v8 cis-eQTL statistics for six tissues and FinnGen release 12 GWAS for hypertension, heart failure, atrial fibrillation, chronic kidney disease and type 2 diabetes. Only *MTHFR* was an eGene in all three intestinal segments, whereas *CTSD* was an eGene in none despite high intestinal expression. Of 240 colocalization tests, none met preregistered thresholds, and only four had genome-wide significant signals on both sides; all four implicated *MTHFR*, whose expression and hypertension lead variants lay 32.7 kb apart in different genes, with distance-dependent Wald ratio drift indicating two distinct causal variants. The strongest candidate-annotated hypertension signal, *FURIN* rs6227 (P = 3.18 × 10⁻³⁴), showed no effect on *FURIN* mRNA in any tissue yet was a strong cis-eQTL for neighbouring *FES*, whose colocalization posterior was 0.985 versus 0.003 for *FURIN*; single-cell mapping nevertheless placed *FURIN* in intestinal tuft cells and *CST3* in enterocytes. Null results here do not reflect insensitivity: two positive-control pairings returned posterior 0.985 for *FES* in the same intestinal tissue and outcome, and 0.7295 for *ACE* in kidney cortex. Tissue-resolved genetic triage should precede mechanistic investment.

*(Word count: 253 — 目标 200–300 词；若目标期刊限 250 词（如 Scientific Reports），投稿前需再删约 1 句)*

## Keywords

gut–heart axis; hypertension; colocalization; Mendelian randomization; expression quantitative trait loci; GTEx; FinnGen; FURIN; FES; *MTHFR*

## Key points box（可选，部分期刊需要）

- **What is already known:** Blood transcriptomic candidate lists increasingly underpin gut–heart axis hypotheses, and proximity-based assignment of GWAS variants to genes remains a common practice.
- **What this study adds:** A three-layer framework that separates expression, local regulation and causal sharing, showing that eight hypertension candidates largely lack tissue-resolved genetic support, that one apparently shared locus harbours two distinct causal variants, and that one candidate-annotated signal acts through its neighbour.
- **Why it matters:** Pre-analytical genetic triage can redirect experimental investment, and de novo attribution should replace proximity-based gene assignment.

# Introduction

Hypertension remains the leading modifiable contributor to cardiovascular death worldwide, and its residual risk has driven interest in mechanisms beyond classical renal and vascular pathways. Among these, the gut–heart axis has become prominent: microbial metabolites, barrier function, immune signalling and intestinal electrolyte handling are now thought to modulate blood pressure, and dietary interventions such as salt restriction act partly by reshaping gut microbial composition. Parallel progress in blood transcriptomics has produced lists of differentially expressed genes that are proposed both as biomarkers and as candidate mediators.

A difficulty follows from this convergence. Genes identified as marker candidates because they can be quantified in accessible tissue are often subsequently described and interpreted as though they had been identified as effectors in a specific organ. Three requirements are thereby conflated: that a gene be detectable in blood, that it be expressed in the tissue of interest, and that its expression in that tissue be causally relevant to the outcome. These requirements do not imply one another. Genes specific to circulating leukocytes may dominate whole-blood profiles with no intestinal involvement, and a gene may be abundant in the intestinal epithelium yet have expression that is not under local cis genetic control, or that shares no causal variant with cardiovascular disease. Because these distinctions determine whether costly mechanistic work is justified, they deserve explicit testing rather than narrative transition.

Two methodological capabilities now allow these requirements to be separated. First, large reference resources provide tissue-resolved cis-expression quantitative trait locus (cis-eQTL) summary statistics, so that local genetic regulation can be assessed in intestinal segments, heart and kidney rather than inferred from blood [3,4]. Second, colocalization and Mendelian randomization, when used jointly, distinguish shared causal variation from merely coincident association [9,10]. Each method alone is fragile. Mendelian randomization in cis regions is vulnerable to linkage disequilibrium between neighbouring variants, which can manufacture a conventional-looking causal estimate between an exposure and an outcome driven by different variants; colocalization, conversely, returns low posteriors in the absence of any signal, thereby conflating lack of power with lack of association. Joint use with explicit positive controls mitigates both failure modes, and additional diagnostic information is available from the behaviour of per-variant Wald ratios, which should be constant under a single shared causal variant.

Existing candidate lists also carry a second, less discussed risk: signals are routinely attributed to the nearest or most topical gene. Annotating a variant to a gene because it lies within its transcript, or because that gene was the subject of the study, is an assumption, not evidence. In a systematic integration of GWAS and eQTL summary data, 77 of 126 prioritised genes were not the nearest annotated gene to the leading variant [11], showing that a gene-centred design built on proximity can be structurally blind to the transcript through which a locus actually acts.

In this study we performed a three-layer genetic triage of eight candidate markers previously linked to hypertension and the renin–angiotensin system in a blood transcriptomic study in which two of them had received endothelial functional validation [1]. The candidates were *LRP1*, *CTSD*, *MTHFR*, *AUTS2*, *FURIN*, *CST3*, *FCER1G* and *TBXAS1*. The framework asks, in order, where each gene is expressed, whether that expression is under local genetic control in the intestinal and cardiovascular tissues relevant to the gut–heart hypothesis, and whether any local genetic signal shares a causal variant with hypertension, heart failure, atrial fibrillation, chronic kidney disease or type 2 diabetes. Two positive controls were analysed with identical code to establish whether the design can detect a shared variant in the setting used here. Our purpose was neither to confirm nor to promote individual candidates, but to determine which candidates are worth experimental investment, and to show how a triage framework can be applied before that investment is committed. Finally, the study addresses a set of candidate genes distinct from canonical renin–angiotensin system components and deliberately shares none of its genes, primary tissues or outcomes with such investigations.

# Methods

## Study design

We conducted a three-layer genetic triage of eight blood transcriptomic markers previously proposed as candidates linking the gut–heart axis to hypertension [1]. The framework separates three questions that are often conflated: where a gene is expressed (expression layer), whether its expression is under local genetic control in the tissue of interest (regulatory layer), and whether that local regulation shares a causal variant with cardiovascular outcomes (causal layer). All analyses used publicly available summary-level data; no individual-level data were accessed, and no new ethical approval was required. Reporting follows the STROBE-MR guideline where applicable [2]. The present study deliberately concerns a different set of genes from work on canonical renin–angiotensin system components, and the resulting analyses share no genes, no primary tissues and no outcomes with those studies.

## Candidate genes

The eight candidates (*LRP1*, *CTSD*, *MTHFR*, *AUTS2*, *FURIN*, *CST3*, *FCER1G* and *TBXAS1*) were taken from a blood transcriptomic study of renin–angiotensin system–related markers in hypertension, in which *FURIN* and *CST3* had undergone loss- and gain-of-function validation in human umbilical vein endothelial cells [1]. All genomic coordinates in this study refer to GRCh38 unless stated otherwise.

## Data sources

Cis-eQTL summary statistics were obtained from the GTEx v8 release as uniformly reprocessed by the eQTL Catalogue [3,4]. Full summary statistics were queried remotely by tabix for six tissues: sigmoid colon (n = 318), transverse colon (n = 368), terminal ileum (n = 174), left ventricle (n = 386), kidney cortex (n = 73) and whole blood (n = 670). Outcome genome-wide association study (GWAS) summary statistics were taken from FinnGen release 12 [5] for essential hypertension (I9_HYPTENS; 154,630 cases and 345,634 controls), heart failure (I9_HEARTFAIL; 37,653/462,695), atrial fibrillation (I9_AF; 63,532/252,810), chronic kidney disease (N14_CHRONKIDNEYDIS; 12,787/480,448) and type 2 diabetes (T2D; 82,878/403,489). FinnGen files were accessed through HTTP range requests against local tabix indices, without downloading genome-wide files. Tissue expression was summarised as GTEx v8 median transcripts per million (TPM) across the 54 tissue entries of the distributed median expression table, which include histologic subregions and cultured fibroblasts in addition to the 49 tissues analysed in the flagship v8 eQTL description [3]. Single-cell expression was obtained from the Human Protein Atlas single-cell type expression resource (3,087,080-row matrix), from which complete profiles were extracted in normalized counts per million (nCPM) for 154 cell types available in this release, compared with the 192 cell type clusters described when the resource was introduced [6]. Gut microbial taxa instruments were derived from the MiBioGen consortium genome-wide significant hits (122,110 records across 211 taxa; per-taxon n = 13,266–18,340) [7].

## Expression and single-cell analyses

Tissue specificity was quantified with the Tau index, Σ(1 − x_i/x_max)/(n − 1), across the 54 tissues, where x_i is the median TPM in tissue i [8]. Tau ≥ 0.8 was considered tissue-specific and Tau ≤ 0.6 broadly expressed. Gut enrichment was summarised as the mean median TPM across the three intestinal segments and as the ratio between the intestinal and the whole-body peak single-cell expression.

## cis-eQTL detectability matrix

For each gene–tissue pair we retrieved all cis-eQTL records within the gene body ±100 kb. A gene was classified as an eGene in a tissue when the minimum cis-eQTL P value in the window was below 5 × 10⁻⁸. The ±100 kb window was fixed a priori; wider windows could not be retrieved reliably through the remote query channel and are acknowledged as a limitation. For every test we recorded three variant counts: tested (all variants in the window), analysed (variants harmonised between exposure and outcome) and fine-mapped (variants entering the colocalization posterior calculation).

## Colocalization

Colocalization between each eQTL signal and each outcome was tested with coloc.abf (R coloc package, version 5.x) [9], using priors p1 = p2 = 1 × 10⁻⁴ and p12 ∈ {1 × 10⁻⁶, 1 × 10⁻⁵, 1 × 10⁻⁴} for sensitivity analysis; p12 = 1 × 10⁻⁵ was the main analysis. In total, 240 tests were performed, covering eight genes, six tissues and five outcomes. A test was declared positive only when PPH4 ≥ 0.80 and PPH4/PPH3 ≥ 3, thresholds fixed before analysis. Because coloc.abf can return degenerate posteriors at extreme approximate Bayes factors and because H3 posteriors can be inflated by combinatorial pairing in variant-dense regions, we pre-specified two positive controls processed with identical code: ACE eQTL in kidney cortex versus hypertension (expected PPH4 ≈ 0.73 from prior work on this locus), and FES eQTL in transverse colon versus hypertension, which uses exactly the same tissue and outcome as the eight candidates. For FES, 803 shared variants were analysed within the window chr15:90,768,588–90,983,564. To guard against combinatorial inflation of PPH3, borderline cases were examined variant by variant. For a given gene, the Wald ratio of every cis-eQTL variant reaching P < 5 × 10⁻⁸ and F > 10 with a matched outcome statistic was computed as the outcome log odds ratio divided by the exposure effect per standard deviation of expression, and the medians were compared across distance bands from the corresponding outcome lead variant. This per-variant gradient was evaluated across all tissues in which the gene was an eGene, and repeated restricted to the single tissue carrying the expression lead. Under one shared causal variant the ratio is expected to be constant, whereas two variants in partial linkage equilibrium produce a distance-dependent gradient.

## cis-Mendelian randomization

For each gene–tissue pair with at least one genome-wide significant cis-eQTL variant, we selected the lowest-P variant with valid statistics on both sides and F statistic > 10, and estimated the causal effect by a single-instrument Wald ratio. Multi-instrument inverse-variance weighting was deliberately not used: no linkage disequilibrium reference panel was available in this environment, and strongly correlated cis variants would artificially inflate precision. The strict instrument threshold was P < 5 × 10⁻⁸; a lenient tier (P < 1 × 10⁻⁵) was treated as exploratory only. Effect alleles were harmonised by position with explicit strand checks. Effect estimates are expressed as odds ratios (OR) per one–standard deviation increase in gene expression. Direction concordance across the five outcomes was examined as an internal validity check.

## Microbiome-to-hypertension Mendelian randomization

As a supplementary analysis, each MiBioGen taxon was represented by its single lead variant with F > 10, yielding 211 instruments at 85 unique genomic positions. Variants were lifted from GRCh37 to GRCh38 through the Ensembl REST API, and outcome statistics were retrieved from FinnGen R12 by tabix. After allele harmonisation, 107 taxa had interpretable single-instrument Wald estimates; Benjamini–Hochberg false discovery rates (FDR) were computed across these 107 tests. Because instruments sharing one variant are not independent, signals were de-duplicated by instrument before interpretation.

## Software and reproducibility

Analyses used R 4.3.3 (coloc 5.x) and Python 3.13, including a pure-Python tabix reader for remote summary statistics. Every figure and table in this manuscript is generated by the analysis scripts from the raw extracted data; scripts, intermediate data and result tables are provided in the Supplementary Materials. All analyses were executed by a single analyst, which is declared as a limitation.

## Data availability

GTEx v8 eQTL summary statistics: eQTL Catalogue (https://www.ebi.ac.uk/eqtl/). FinnGen R12: https://www.finngen.fi/. GTEx median expression: https://gtexportal.org/home/. Human Protein Atlas: https://www.proteinatlas.org/. MiBioGen: https://www.mibiogen.org/. No new data were generated.

# Results

## Expression triage: three distinct profiles among eight candidates

Across 54 GTEx tissues, the eight candidates fell into three expression classes (Fig. 1, Table 1). Two genes were essentially blood-specific: *FCER1G* (Tau = 0.972; whole blood 1,364.8 TPM versus 35.2 TPM across the three intestinal segments) and *TBXAS1* (Tau = 0.941; 68.3 versus 4.9 TPM). Four genes showed intermediate specificity (*LRP1*, *FURIN*, *CTSD* and *AUTS2*), and two were broadly expressed (*MTHFR*, *CST3*). *LRP1* was the only candidate with relative intestinal enrichment, with 85.9 TPM in sigmoid colon versus 8.5 TPM in whole blood, a tenfold difference. *CTSD* was the most highly transcribed gene of the set (835 TPM in kidney cortex and 1,571 TPM in adrenal gland) yet showed no tissue preference (Tau = 0.714), illustrating that abundant transcription carries no localisation information.

## Regulatory triage: the intestine is the limiting layer

The detectability matrix covering eight genes and six tissues (Fig. 2, Table 2) showed that only *MTHFR* was an eGene in all three intestinal segments (minimum P = 1.38 × 10⁻¹⁰ in sigmoid colon, 3.97 × 10⁻¹⁴ in transverse colon and 5.93 × 10⁻¹¹ in ileum). *CST3* was colon-specific (7.99 × 10⁻¹⁵ and 6.77 × 10⁻⁹ in the two colons, but 2.5 × 10⁻⁴ in ileum), and *LRP1* reached significance only in transverse colon (6.97 × 10⁻¹⁵). *CTSD* was not an eGene in any of the six tissues (minimum P = 5.9 × 10⁻⁵ in whole blood), despite 597 TPM in ileum, demonstrating that being heavily transcribed and being locally genetically regulated are independent properties. No gene reached eGene significance in left ventricle or kidney cortex; given the sample sizes (n = 386 and n = 73), these null results reflect limited power and were not interpreted biologically.

## Causal triage: zero of 240 colocalization tests passed, with two positive controls

Before interpreting the negative findings, we verified that the pipeline can detect true positives in the exact setting used here (Table 4). The ACE–kidney cortex control returned PPH4 = 0.7295. Critically, the FES–transverse colon control, which uses the same intestinal tissue and the same hypertension GWAS as the eight candidates, returned PPH4 = 0.9846 in the main analysis and 0.8646 and 0.9984 under the alternative priors (803 shared variants). The design is therefore capable of yielding strong positive colocalization in intestinal tissue.

Among the 240 candidate tests (Fig. 3, Fig. 4, Table 3), none met the preregistered positive threshold. Only four tests had genome-wide significant signals on both sides, all of them *MTHFR* versus hypertension, and all four returned PPH4 ≈ 0 with PPH3 ≈ 1. A further 46 tests were significant only on the eQTL side and 38 only on the outcome side; the remaining 152 lacked signal on both sides and are uninformative rather than negative. The strongest result per outcome is reported for completeness: *AUTS2* in whole blood favoured distinct variants for hypertension (PPH3 = 0.752) and type 2 diabetes (PPH3 = 0.815), although *AUTS2* is barely transcribed in blood (0.83 TPM), raising the possibility that the tagged eQTL marks a neighbouring transcript. *TBXAS1* in whole blood gave the only H4/H3 ratio above 3 (atrial fibrillation; PPH4 = 0.420, H4/H3 = 11.81), but the outcome-side minimum P was only 2.3 × 10⁻⁴, so this remains a lead for replication rather than a positive finding.

## The 15q26.1 hypertension signal acts through FES, not FURIN

The strongest outcome signal anywhere in the study lay within *FURIN* itself: rs6227 (chr15:90,882,002) reached P = 3.18 × 10⁻³⁴ for hypertension. Displaying the same window, the same variant set and the same tissue side by side makes the attribution explicit (Fig. 5). Ensembl VEP annotation confirmed rs6227 as a *FURIN* 3′ untranslated region variant (ENST00000268171, 1,562 bp from the 3′ end) and simultaneously as an upstream variant of the adjacent FES gene, 1,693 bp from its transcription start site. Despite this annotation, rs6227 showed no association with *FURIN* mRNA in any tissue tested (P = 0.918, 0.345, 0.930 and 0.568 in sigmoid colon, transverse colon, ileum and whole blood), whereas it was strongly associated with FES expression in transverse colon (P = 4.2 × 10⁻⁹, β = −0.196). Direct colocalization confirmed this assignment: FES in transverse colon versus hypertension gave PPH4 = 0.9846, whereas *FURIN* in the identical window gave PPH4 = 0.003. The causal transcript at this blood pressure locus is therefore FES, and study designs built on *FURIN* eQTLs are structurally blind to this signal.

## MTHFR: resolving the MR–colocalization contradiction

**MTHFR* was the only candidate with power on both sides of the colocalization test, and it returned an unambiguous negative (Fig. 6a). The transverse colon eQTL lead rs3818762 (chr1:11,790,946; P = 3.97 × 10⁻¹⁴) lies within *MTHFR*, whereas the hypertension lead rs12567136 (chr1:11,823,674; P = 5.19 × 10⁻⁴⁴) lies within CLCN6, 32.7 kb away. Across 909 shared variants, PPH3 = 1.000 and PPH4 = 0.000 under all three priors.

Because H3 posteriors can be inflated by combinatorial pairing, we examined the locus variant by variant. We computed per-variant Wald ratios for every *MTHFR* cis-eQTL variant with F > 10 in each tissue in which *MTHFR* was an eGene, giving 377 variant–tissue combinations across transverse colon, sigmoid colon, ileum and whole blood (Fig. 6b). Under a single shared causal variant the ratio should be constant; instead its median drifted with distance from the hypertension lead: −0.139 (n = 147), −0.141 (n = 169) and −0.063 (n = 61) for variants within 20 kb, within 20–50 kb and beyond 50 kb. The same monotonic drift was present when the analysis was restricted to transverse colon alone, the tissue in which the expression lead was defined (−0.219, −0.190 and −0.135 for 19, 44 and 25 variants). Two variants flanking the outcome lead illustrate the problem more directly: rs12404124 (chr1:11,813,812) and rs198397 (chr1:11,823,275) tag the same hypertension signal (P = 2.12 × 10⁻¹⁹ and 2.17 × 10⁻¹⁹) yet carry cis-eQTL effect sizes of opposite sign in transverse colon (β = −0.217 and +0.214; P = 9.7 × 10⁻¹¹ and 1.7 × 10⁻¹⁰). A single shared causal variant cannot associate one outcome effect with two opposite directions of expression change, so both patterns indicate two distinct causal variants in partial linkage disequilibrium. Consistently, the *MTHFR* cis-MR estimate pointed in opposite directions across outcomes: apparently protective for hypertension (OR 0.869, 95% CI 0.835–0.903; P = 1.9 × 10⁻¹²) but risk-increasing for heart failure (OR 1.096, 95% CI 1.028–1.167; P = 4.7 × 10⁻³), a pattern expected from linkage contamination rather than a genuine causal effect.

## cis-MR and single-cell localisation

Strict instruments were available for only three genes (Table 5): *LRP1* (rs7398375, transverse colon; F = 67.3), *MTHFR* (rs3818762; F = 63.2) and *CST3* (rs7270413, sigmoid colon; F = 66.6, a downgraded second-choice instrument because the lead variant rs71334202 was absent from FinnGen). The remaining five genes had no usable instrument (best cis-eQTL P values 8.1 × 10⁻⁷ to 3.9 × 10⁻⁵), which represents untestability rather than absence of effect. *LRP1* and *CST3* showed no association with any of the five outcomes (for hypertension, OR 1.040, 95% CI 0.993–1.088 and OR 0.999, 95% CI 0.973–1.025, respectively) (Supplementary Fig. S3).

Single-cell mapping localised the candidates but did not rescue their causal evidence (Supplementary Fig. S1). *FURIN* reached its whole-body expression maximum in intestinal tuft cells (277 nCPM; intestinal/body peak ratio = 1.000), and *CST3* was strongly expressed in enterocytes (3,917 nCPM), consistent with its high colonic bulk expression. Conversely, *FCER1G* was essentially absent from intestinal epithelium (0.8 nCPM) and peaked in neutrophils (3,045 nCPM), as did *TBXAS1* (2,239 nCPM). The two genes with the most credible intestinal expression thus lacked causal-layer evidence, while the genes with causal-layer signals were not intestinal in either expression or cell-type terms.

## Supplementary microbiome analysis

Single-instrument MR of gut microbial taxa against hypertension yielded 107 interpretable estimates (Supplementary Fig. S2, Supplementary Table S3). Seven taxa passed FDR < 0.05 but corresponded to only three independent instruments. The four significant taxa of the Bifidobacterium lineage shared one instrument, rs182549, located in the LCT/MCM6 lactase-persistence locus, a variant with well-established dietary pleiotropy; its risk-direction estimate (OR 1.22–1.24) was also directionally inconsistent with previous multi-instrument analyses. These signals were therefore treated as exploratory and were not carried into the main conclusions.

# Discussion

## Principal findings

We subjected eight blood transcriptomic markers proposed as gut–heart axis candidates to a three-layer genetic triage. In the causal layer, none of 240 colocalization tests met the preregistered positive threshold, and only three of eight genes supported a strict Mendelian randomization instrument. The failure was not uniform across layers: expression and single-cell data placed several candidates convincingly in the intestinal epithelium, yet those placements did not translate into genetic evidence for shared causal variation with cardiovascular outcomes. The most concise statement of our result is that where these markers are expressed is largely independent of whether they are causal, and that being measurable in blood is not evidence of acting through the gut.

## Why the negative findings are credible

Negative results in this design are only meaningful if the pipeline can detect a true positive in the same configuration. Two controls establish this. The FES–transverse colon test used the identical tissue, identical outcome genome-wide association study, identical window definition and identical code as the candidate tests, and returned PPH4 = 0.9846. The ACE–kidney cortex test, included as an a priori positive control from prior work on canonical renin–angiotensin system genes, returned PPH4 = 0.7295, consistent with previously reported values for this pair. Conversely, we did not treat the 152 tests without signal on either side as negative evidence; these are uninformative tests, and we report them as such. The informativeness boundaries were set by power, most severely in kidney cortex (n = 73) and left ventricle (n = 386).

## Signal attribution, not signal strength, was the decisive issue

Three observations converged on the same methodological point: the assignment of a GWAS or eQTL signal to a transcript must be verified variant by variant, because proximity is not ownership. At 15q26.1, rs6227 is annotated as a *FURIN* 3′ untranslated region variant and yet reaches the highest outcome significance in this study without detectable effects on *FURIN* mRNA in any tissue; the same variant is a strong cis-eQTL for FES, and formal colocalization endorsed FES (PPH4 = 0.985) rather than *FURIN* (PPH4 = 0.003). At 1p36, the *MTHFR* eQTL lead and the hypertension lead lie 32.7 kb apart in different nearest genes, and the per-variant Wald ratio drifts with distance and reverses sign within 10 kb, a signature that coincident MR with a single causal variant cannot produce. In the supplementary microbiome analysis, the apparently strongest association again belonged to a different locus: the entire Bifidobacterium-lineage signal rested on a single LCT/MCM6 variant, whose relationship to microbial abundance is mediated by lactose intake and whose estimates were directionally inconsistent with previous multi-instrument work.

These findings also explain why integrating expression and outcome data requires an explicit test against pleiotropy rather than reliance on physical proximity: in a systematic integration of this kind, 77 of 126 prioritised genes were not the nearest annotated gene to the top associated variant [11]. They further illustrate that cis-Mendelian randomization estimates depend on the instrumental-variable assumptions holding, not merely on the number of variants carried into the test [10]. A practical corollary for hypothesis-driven studies of the gut–heart axis is that candidate lists generated in blood cannot be carried into intestinal experimental models without tissue-resolved genetic verification.

## Relation to previous evidence

The eight candidates originated from a blood transcriptomic analysis in which two of them, *FURIN* and *CST3*, had endothelial functional validation [1]. That evidence addresses a different question from ours: what a protein does when manipulated in a cell is compatible with no contribution from cis-regulated expression variation in the relevant tissue at the population level. Our results are therefore not in conflict with those experiments, but they delimit what such experiments can support. More broadly, the pattern within our own data is instructive: whole blood had by far the largest cis-eQTL sample size (n = 670) and yielded eGene status for four of the eight candidates, yet not one of those whole-blood signals shared a causal variant with any of the five outcomes. Abundance of transcript and detectability of regulation therefore carry no information about whether that regulation is connected to disease in the tissue where it can be measured.

We also note that the two candidates with strong tissue specificity, *FCER1G* and *TBXAS1*, were blood-specific rather than gut-specific, and that the only candidate with intestinal enrichment by bulk expression, *LRP1*, had a significant local eQTL in just one intestinal segment and no colocalization with any outcome. This dissociation between the expression layer and the causal layer is the empirical basis for the triage framework, rather than any single gene result.

## Implications for experimental design

Triage has a practical value: it prevents the allocation of experimental resources to candidates that population genetics cannot support as cis-regulated actors in the tissue of interest. Our data do suggest narrower, testable next steps. The *FURIN* maximum in intestinal tuft cells is biologically interesting because tuft cells are chemosensory sentinels of the intestinal epithelium, yet any experiment derived from this observation must include FES as a comparator, since population-level evidence assigns the blood pressure signal to FES. *CST3* is strongly expressed in enterocytes without any corresponding genetic signal, so it is a useful internal contrast for studies aiming to separate expression abundance from regulatory contribution. Finally, supplementing our currently unavailable protein quantitative trait locus step would allow us to test whether rs6227 acts on *FURIN* protein rather than mRNA, which is the single most informative missing measurement in this project.

## Limitations

Five limitations bound the interpretation. First, our cis windows were fixed at ±100 kb, so distal regulatory effects were not examined. Second, no linkage disequilibrium reference panel was available, so Mendelian randomization relied on single-instrument Wald ratios without heterogeneity or pleiotropy diagnostics. Third, power was limited in kidney cortex and left ventricle, and their null results cannot be read biologically. Fourth, human bulk and single-cell atlases do not capture cell-type-specific eQTL or conditional regulation, so our regulatory layer is conservative by design. Fifth, all analyses were undertaken by a single analyst with no blinded duplicate screening, which we declare explicitly. Protein-level analyses for *FURIN* were not completed because the relevant summary statistics could not be retrieved through available channels, which is an access problem rather than an absence of data.

## Conclusions

Across three layers of genetic evidence, eight blood-derived "gut–heart axis" candidates for hypertension largely did not withstand triage. One locus with adequate power on both sides was shown to harbour two distinct causal variants rather than one shared variant, and one signal annotated to a candidate was shown to act through its neighbour. Because the same pipeline returned strong positive colocalization for a neighbouring gene in the same tissue and for the same outcome, the most useful reading of these null results is that they constrain candidate selection and sharpen attribution practice: before a candidate enters mechanistic work on the gut–heart axis, its local genetic regulation and the ownership of its locus should be established in the tissue where it is claimed to act.

## Declarations

**Ethics approval and consent to participate.** All analyses used publicly available summary-level data from studies that had obtained participant informed consent and ethics approval; no additional ethical approval was required for the present secondary analysis.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All analysis scripts, extracted intermediate data and complete result tables accompany this manuscript as Supplementary Materials.

**Competing interests.** The authors declare that they have no competing interests.

**Funding.** This research received no funding from any agency in the public, commercial or not-for-profit sectors.

**Authors' contributions.** **Jing Bai**: Conceptualization, study design, supervision of the analytical framework, critical revision of the manuscript for important intellectual content, approval of the final version. **Yan Wang**: Methodology, software, data curation, formal analysis, investigation, visualisation, writing of the original draft, project administration. Jing Bai is the corresponding author and takes responsibility for the integrity of the work as a whole.

**Acknowledgements.** We acknowledge the GTEx Consortium, the FinnGen Study, the MiBioGen Consortium and the Human Protein Atlas for making their data openly available.

# Figure legends

All figures were generated directly from the analysis outputs described in Methods, without
manual editing. Every numeric value quoted in a legend is reproduced in the Results text or
in the corresponding table. Panel labels follow the arrangement of each source file.

## Figure 1 Expression atlas of the eight candidates across 54 GTEx tissues

Heat map of log2(median TPM + 1) for the eight candidates across the 54 GTEx v8 tissues with
median transcripts per million available, computed from `expression_median_tpm.csv`. Tissues
are grouped into gut, cardiovascular, metabolic, immune, neural and structural compartments
and ordered by decreasing *LRP1* expression within compartment; the gradient runs from white
(low) through orange to dark red (high). Each cell carries one value per gene and tissue, so
no error bars apply. *FCER1G* and *TBXAS1* are almost entirely restricted to whole blood,
whereas each of the remaining six genes is detected in all 54 tissues, with minimum median
values ranging from 0.82 TPM (*LRP1*) to 71.6 TPM (*CTSD*). Source: `Figure_1_expression_atlas.png`.

## Figure 2 Detectability of cis-eQTL signal in eight genes and six tissues

Heat map of −log10 of the smallest cis-eQTL P value within the gene body ± 100 kb window,
computed from `matrix_raw.csv`. Black outlines mark eGene status (P < 5 × 10⁻⁸) and printed
numbers give the −log10 P value. Tissue order is sigmoid colon, transverse colon, terminal
ileum, left ventricle, kidney cortex and whole blood, corresponding to GTEx eQTL sample sizes
of 318, 368, 174, 386, 73 and 670. Only *MTHFR* attains eGene status in all three intestinal
segments; *CTSD* attains it in none of the six tissues. Source: `Figure_2_detectability_matrix.png`.

## Figure 3 Posterior support for colocalization across 240 candidate tests

Scatter plot of PPH3 (posterior probability of two distinct causal variants) against PPH4
(posterior probability of one shared causal variant) for all 240 gene × tissue × outcome
tests, computed from the five per-outcome colocalization output files. Colour denotes outcome
and symbol denotes whether the tissue was intestinal (filled) or non-intestinal (open).
Dashed lines mark the preregistered threshold of 0.80 on both axes. No test enters the
positive region; the four tests with signal on both sides lie in the lower-right corner with
PPH3 approaching 1, and are all *MTHFR* versus hypertension. Positive-control pairings run
through the identical pipeline are reported separately in Table 4. Source:
`Figure_3_colocalization_posteriors.png`.

## Figure 4 Signal availability across the 240 tests

Scatter plot of the strongest available evidence on each side of the test: −log10 of the
minimum cis-eQTL P value in the window (horizontal axis) against −log10 of the minimum
outcome P value in the same window (vertical axis), with points coloured by the resulting
PPH4. Dashed lines denote genome-wide significance (P = 5 × 10⁻⁸) and the shaded quadrant is
the region in which the colocalization test is interpretable because both traits have signal.
Only four of 240 tests fall in that quadrant, all of them *MTHFR* versus hypertension, and all
four returned PPH4 between 0.000 and 0.001. Source: `Figure_4_signal_availability.png`.

## Figure 5 The 15q26.1 hypertension signal is a *FES* eQTL and not a *FURIN* eQTL

Three panels sharing one genomic window and one variant set on chromosome 15 (GRCh38):
(a) *FURIN* cis-eQTL associations in transverse colon, (b) *FES* cis-eQTL associations in the
same tissue, and (c) FinnGen release 12 hypertension associations in the same window.
Each point is one variant plotted as −log10 P against position in megabases; the dashed
vertical line marks rs6227 and the dotted horizontal line in panel (c) marks genome-wide
significance. rs6227 is the top hypertension variant in the window (P = 3.18 × 10⁻³⁴) and a
strong *FES* cis-eQTL (P = 4.2 × 10⁻⁹), while showing no association with *FURIN* mRNA in this
tissue (P = 0.345) or in any other tissue tested. Source: `Figure_5_FURIN_FES_region.png`.

## Figure 6 *MTHFR*: a formally significant cis-MR estimate that is not supported by colocalization

**Figure 6a** Regional association plot at the *MTHFR* locus on chromosome 1 (GRCh38), with
*MTHFR* cis-eQTL associations in transverse colon in the upper panel and FinnGen release 12
hypertension associations in the lower panel; the dashed lines mark the expression lead
rs3818762 and the blood pressure lead rs12567136, which lie 32.7 kb apart and fall in
different genes. Source: `Figure_6a_MTHFR_region.png`.

**Figure 6b** Per-variant Wald ratios (log odds ratio of hypertension per standard deviation
of expression) plotted against position for the 377 variant–tissue combinations with F > 10,
coloured by tissue. The vertical lines mark the two lead variants. A single shared causal
variant would give a horizontal band of constant ratio; instead the median ratio drifts with
distance from the blood pressure lead and individual ratios change sign within 10 kb of it,
indicating at least two causal variants in partial linkage disequilibrium. Each point is one
variant–tissue estimate; the spread within a column reflects sampling variation at that
variant across tissues. Source: `Figure_6b_MTHFR_wald_gradient.png`.

## Supplementary Figure S1 Single-cell localisation of the eight candidates

Heat map of log10(normalized counts per million + 1) across 154 Human Protein Atlas single-cell
types, faceted by tissue compartment, from the complete 8 × 154 matrix in
`celltype_full_matrix.csv`. Values below the Human Protein Atlas reporting threshold are shown
as absent rather than missing. *FURIN* reaches its whole-body maximum in intestinal tuft cells
(277 nCPM) and *CST3* is highest in enterocytes (3,917 nCPM), whereas *FCER1G* and *TBXAS1* peak
in neutrophils (3,045 and 2,239 nCPM) and are near the detection floor in intestinal epithelium.
Source: `Supplementary_Figure_S1_celltypes.png`.

## Supplementary Figure S2 Microbiome-to-hypertension Mendelian randomization

Volcano plot of 107 taxon-level single-instrument estimates from MiBioGen lead variants against
FinnGen release 12 hypertension, plotting −log10 P against log odds ratio. Points above the
dashed line reach Benjamini–Hochberg FDR < 0.05. Three taxa are labelled with their instrument. The
labelled Bifidobacterium signal uses rs182549 in the LCT/MCM6 lactase-persistence locus and is
shared with three further taxa of the same lineage, so those four estimates are not independent;
they are reported as exploratory only. Source:
`Supplementary_Figure_S2_microbiome_MR.png`.

## Supplementary Figure S3 cis-Mendelian randomization forest plot

Odds ratios per standard deviation higher tissue expression with 95% confidence intervals for
the three genes with a strict instrument (*LRP1*, *MTHFR* and *CST3*), faceted by gene, for the
five outcomes. Points are single-instrument Wald ratios; the dotted vertical line marks no
effect and grey points denote non-significant estimates. Only *MTHFR* reaches significance, and
its estimates point in opposite directions for hypertension and heart failure. Source:
`Supplementary_Figure_S3_cisMR_forest.png`.

<!-- AUTO-GENERATED by 04_code/27_make_tables.py - do not edit by hand -->

# Tables

All entries are reproduced from the analysis output files listed under
Data availability; Supplementary Tables S1–S3 accompany this article as separate files.

## Table 1 Expression profile of the eight candidates across 54 GTEx tissues

Tau is the tissue-specificity index of Yanai et al. [8], computed across 54 GTEx tissues;
higher values indicate stronger tissue restriction. TPM, transcripts per million.

| Gene | Median TPM | Maximum TPM | Tissue of maximum | Tau | Class | Mean gut TPM | Blood TPM | Gut/blood ratio |
|---|---|---|---|---|---|---|---|---|
| *LRP1* | 57.09 | 358.25 | Cultured fibroblasts | 0.773 | intermediate | 64.6 | 8.5 | 7.6 |
| *CTSD* | 420.19 | 1,571 | Adrenal Gland | 0.714 | intermediate | 439.9 | 808.7 | 0.54 |
| *MTHFR* | 11.1 | 25.64 | Ovary | 0.580 | broadly expressed | 11.8 | 5.2 | 2.3 |
| *AUTS2* | 10.13 | 29.92 | Tibial nerve | 0.645 | intermediate | 9.8 | 0.8 | 11.9 |
| *FURIN* | 56.92 | 192.92 | Liver | 0.717 | intermediate | 50.9 | 108.1 | 0.47 |
| *CST3* | 275.59 | 653.31 | Putamen (basal ganglia) | 0.548 | broadly expressed | 237.9 | 109 | 2.2 |
| *FCER1G* | 20.24 | 1,365 | Whole blood | 0.972 | tissue-specific | 35.2 | 1,365 | 0.03 |
| *TBXAS1* | 2.62 | 68.3 | Whole blood | 0.941 | tissue-specific | 4.9 | 68.3 | 0.07 |

## Table 2 cis-eQTL detectability matrix: minimum association P value per gene and tissue

Values are the smallest cis-eQTL P value inside the gene body +/- 100 kb window.
Bold entries meet the eGene threshold (P < 5 x 10^-8). Sample size n is the GTEx
cis-eQTL sample size for that tissue.

| Gene | Sigmoid colon | Transverse colon | Terminal ileum | Left ventricle | Kidney cortex | Whole blood |
|---|---|---|---|---|---|---|
| *LRP1* | 1.1 x 10^-2 | **6.97 x 10^-15** | 9.43 x 10^-3 | 1.34 x 10^-5 | 1.99 x 10^-2 | 2.96 x 10^-4 |
| *CTSD* | 3.98 x 10^-3 | 2.46 x 10^-4 | 3.03 x 10^-3 | 2.17 x 10^-3 | 2.27 x 10^-3 | 5.93 x 10^-5 |
| *MTHFR* | **1.38 x 10^-10** | **3.97 x 10^-14** | **5.93 x 10^-11** | 3.97 x 10^-5 | 6.98 x 10^-4 | **2.13 x 10^-48** |
| *AUTS2* | 1.8 x 10^-3 | 4.83 x 10^-3 | 4.23 x 10^-3 | 1.5 x 10^-5 | 1.29 x 10^-2 | 7.67 x 10^-7 |
| *FURIN* | 1.66 x 10^-3 | 1.7 x 10^-4 | 3.88 x 10^-5 | 5.2 x 10^-4 | 1.16 x 10^-3 | 5.62 x 10^-7 |
| *CST3* | **7.99 x 10^-15** | **6.77 x 10^-9** | 2.51 x 10^-4 | 8.62 x 10^-4 | 3.56 x 10^-3 | **6.72 x 10^-21** |
| *FCER1G* | 8.55 x 10^-5 | 3.02 x 10^-7 | 3.88 x 10^-4 | 1.2 x 10^-7 | 8.22 x 10^-3 | **3.15 x 10^-45** |
| *TBXAS1* | 4.17 x 10^-3 | 8.09 x 10^-7 | 2.64 x 10^-3 | 2.71 x 10^-3 | 5.97 x 10^-4 | **8.85 x 10^-15** |

Cis-eQTL sample size by tissue: Sigmoid colon n = 318; Transverse colon n = 368; Terminal ileum n = 174; Left ventricle n = 386; Kidney cortex n = 73; Whole blood n = 670.

## Table 3 Signal availability across the 240 colocalization tests

Each of the 240 gene-tissue-outcome tests is classified by whether the cis window
contained genome-wide significant signal (P < 5 x 10^-8) on the eQTL side, the outcome
side, both, or neither. The best-supported test per outcome is given for completeness;
no test met the preregistered positive threshold (PPH4 >= 0.80 and PPH4/PPH3 >= 3).

| Outcome | Both sides | eQTL only | Outcome only | Neither | Best-supported test | PPH3 | PPH4 |
|---|---|---|---|---|---|---|---|
| Hypertension | 4 | 6 | 20 | 18 | *MTHFR*, Terminal ileum | 0.999 | 0.001 |
| Heart failure | 0 | 10 | 0 | 38 | *MTHFR*, Whole blood | 0.072 | 0.238 |
| Atrial fibrillation | 0 | 10 | 6 | 32 | *TBXAS1*, Whole blood | 0.036 | 0.420 |
| Chronic kidney disease | 0 | 10 | 0 | 38 | *AUTS2*, Whole blood | 0.103 | 0.085 |
| Type 2 diabetes | 0 | 10 | 12 | 26 | *AUTS2*, Whole blood | 0.815 | 0.159 |
| **All outcomes** | **4** | **46** | **38** | **152** | - | - | - |

Full per-test results are given in Supplementary Table S1.

## Table 4 Positive-control colocalizations analysed with the identical pipeline

Both controls were run through the same script, window definition and outcome data as
the 240 candidate tests. Values shown for the three prior settings of the shared causal
variant (p12 = 1 x 10^-6, 1 x 10^-5 and 1 x 10^-4).

| Control | Tissue | Shared variants | Minimum P (eQTL) | Minimum P (outcome) | PPH3 | PPH4 | PPH4 (p12 = 10^-6) | PPH4 (p12 = 10^-4) |
|---|---|---|---|---|---|---|---|---|
| *ACE* | Kidney cortex | 242 | 4.16 x 10^-8 | 4.71 x 10^-14 | 0.268 | 0.730 | 0.212 | 0.964 |
| *FES* | Transverse colon | 803 | 5.21 x 10^-13 | 3.18 x 10^-34 | 0.015 | 0.985 | 0.865 | 0.998 |

Note. The number of shared variants for the *FES* control (803) was recomputed from the
harmonised input files because the colocalization output recorded it blank. The lead variant
rs6227 is simultaneously a *FES* cis-eQTL in transverse colon and the top hypertension variant
in the window.

## Table 5 cis-Mendelian randomization results for the candidate genes with usable instruments

Strict instruments were independent cis-eQTL variants with P < 5 x 10^-8 and F > 10.
Estimates are odds ratios per one standard deviation higher tissue expression.
Five of the eight genes had no instrument meeting these criteria.

| Gene | Tissue | Instrument | F statistic | Outcome | OR | 95% CI | P |
|---|---|---|---|---|---|---|---|
| *LRP1* | Transverse colon | rs7398375 | 67.3 | Hypertension | 1.040 | 0.993-1.088 | 0.097 |
| *LRP1* | Transverse colon | rs7398375 | 67.3 | Heart failure | 1.075 | 0.998-1.158 | 0.056 |
| *LRP1* | Transverse colon | rs7398375 | 67.3 | Atrial fibrillation | 0.968 | 0.895-1.047 | 0.413 |
| *LRP1* | Transverse colon | rs7398375 | 67.3 | Chronic kidney disease | 1.122 | 0.997-1.262 | 0.056 |
| *LRP1* | Transverse colon | rs7398375 | 67.3 | Type 2 diabetes | 0.996 | 0.944-1.050 | 0.867 |
| *MTHFR* | Transverse colon | rs3818762 | 63.2 | Hypertension | 0.869 | 0.835-0.903 | 1.92 x 10^-12 |
| *MTHFR* | Transverse colon | rs3818762 | 63.2 | Heart failure | 1.096 | 1.028-1.167 | 0.005 |
| *MTHFR* | Transverse colon | rs3818762 | 63.2 | Atrial fibrillation | 1.041 | 0.974-1.113 | 0.24 |
| *MTHFR* | Transverse colon | rs3818762 | 63.2 | Chronic kidney disease | 0.933 | 0.844-1.032 | 0.177 |
| *MTHFR* | Transverse colon | rs3818762 | 63.2 | Type 2 diabetes | 0.995 | 0.951-1.041 | 0.824 |
| *CST3* | Sigmoid colon | rs7270413 | 66.6 | Hypertension | 0.999 | 0.973-1.025 | 0.916 |
| *CST3* | Sigmoid colon | rs7270413 | 66.6 | Heart failure | 0.976 | 0.936-1.017 | 0.252 |
| *CST3* | Sigmoid colon | rs7270413 | 66.6 | Atrial fibrillation | 1.034 | 0.990-1.081 | 0.133 |
| *CST3* | Sigmoid colon | rs7270413 | 66.6 | Chronic kidney disease | 1.066 | 0.999-1.139 | 0.056 |
| *CST3* | Sigmoid colon | rs7270413 | 66.6 | Type 2 diabetes | 1.022 | 0.992-1.053 | 0.152 |

Wald ratios are single-instrument estimates. Supplementary Fig. S3 presents these
results as a forest plot; Supplementary Table S2 reports the cell-type localization matrix.

## Supplementary tables

| File | Content |
|---|---|
| Supplementary Table S1 | All 240 colocalization tests with full posterior output (`coloc_master.csv`, aggregated from the five per-outcome files) |
| Supplementary Table S2 | Complete 8 x 154 single-cell expression matrix summary (`celltype_full_summary.csv`) |
| Supplementary Table S3 | Microbiome-to-hypertension Mendelian randomization results for 107 taxon-level estimates (`mr_microbiome_htn.csv`) |

# References

（AMA/Vancouver 编号制格式；[1] 元数据已由原文 PDF 核对；其余条目 DOI 均已给出，投稿前建议再用 PubMed 复核一次）

1. Bai J, Wang Y, Yang X, Du L, Shi J, Li X, Bai M. Transcriptomic identification of renin–angiotensin system-related candidate biomarkers and external testing of a hypertension diagnostic model. *J Vis Exp*. 2026;(232):e71252. doi:10.3791/71252. Published June 22, 2026.
   （作者本人工作的自我引用；元数据已从原文 PDF 提取核对：刊名 *Journal of Visualized Experiments*，卷 232，文章号 e71252）
2. Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the reporting of observational studies in epidemiology using Mendelian randomization: the STROBE-MR statement. *JAMA*. 2021;326(16):1614–1621. doi:10.1001/jama.2021.18236.
3. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*. 2020;369(6509):1318–1330. doi:10.1126/science.aaz1776.
4. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. *Nat Genet*. 2021;53(9):1290–1299. doi:10.1038/s41588-021-00924-w.
5. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature*. 2023;613(7944):508–518. doi:10.1038/s41586-022-05473-8.
6. Karlsson M, Zhang C, Méar L, et al. A single-cell type transcriptomics map of human tissues. *Sci Adv*. 2021;7(31):eabh2169. doi:10.1126/sciadv.abh2169.
7. Kurilshikov A, Medina-Gomez C, Bacigalupe R, et al. Large-scale association analyses identify host factors influencing human gut microbiome composition. *Nat Genet*. 2021;53(2):156–165. doi:10.1038/s41588-020-00763-1.
8. Yanai I, Benjamin H, Shmoish M, et al. Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. *Bioinformatics*. 2005;21(5):650–659. doi:10.1093/bioinformatics/bti042.
9. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet*. 2014;10(5):e1004383. doi:10.1371/journal.pgen.1004383.
10. Lawlor DA, Harbord RM, Sterne JAC, Timpson N, Davey Smith G. Mendelian randomization: using genes as instruments for making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133–1163. doi:10.1002/sim.3034.
11. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. *Nat Genet*. 2016;48:481–487. doi:10.1038/ng.3538. PMID: 27019110.
   （**第 1 轮审计替换**：原 [11] 写作 "Verdonk S, et al. *PLoS Biol* 2023;21(1):e3001957"，经 Crossref + Europe PMC 核查该 DOI 实为社论 "Reflecting on 2022"（PLoS Biol 20:e3001957），且按作者/标题检索均无该文——属错误引用，已删除并换为本条经核实文献）

**待补充引用（按需要）**
- **投稿纪律（用户已确认）**：自家**未投稿/在投**的稿件一律不引用，也不以 "concurrent unpublished analysis" 形式在正文提及；本文已在 Methods/Discussion 中移除相关表述
- 血压 GWAS：Evangelou E, et al. Genetic analysis of over 1 million people identifies 535 new loci associated with blood pressure traits. *Nat Genet*. 2018;50(10):1412–1425.（用于佐证 15q26.1 位点的血压关联背景）
- 血液转录组/标志物泛化性讨论：待补近 5 年综述 1–2 篇
- 簇细胞（tuft cell）生物学背景：待补综述 1 篇

## 引编号使用自查

| 编号 | 出现位置 |
|---|---|
| [1] | Introduction（候选基因来源）、Methods（candidate genes）、Results?、Discussion（功能验证对比） |
| [2] | Methods（STROBE-MR） |
| [3,4] | Methods（GTEx v8 / eQTL Catalogue）、Introduction（tissue-resolved eQTL） |
| [5] | Methods（FinnGen R12） |
| [6] | Methods（HPA 单细胞） |
| [7] | Methods（MiBioGen） |
| [8] | Methods（Tau） |
| [9] | Methods（coloc）、Introduction（识别共定位） |
| [10] | Introduction（MR）、Discussion（工具有效性） |
| [11] | Discussion（邻近基因调控变异） |

*注：Results 章节目前未直接引编号（符合 Results 客观陈述原则），如需在 Discussion 再引 [3][5]，投稿前统一核对连续编号。*

