# -*- coding: utf-8 -*-
"""Extract 1000G EUR genotypes for the two candidate loci and save dosages.

The 1000G high-coverage panel VCF carries a genotype column per sample, so a
region query is expensive.  Only the two loci that matter are pulled, and only the
633 EUR samples are kept.

Output (08_ld/):
  gt_<locus>.tsv    chr, pos, ref, alt, maf_eur, then one dosage (0/1/2) per EUR sample
  eur_samples.txt   sample order used in the columns

Usage:  python 01_extract_ld.py [locus ...]
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, r"F:\longxiawenjian\gaoxueyaqiguanzhou\P0_mr")
from _tabix import Tabix                                    # noqa: E402

BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
        "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

# gene body +/- 100 kb, GRCh38, identical to the eQTL windows in the manuscript
LOCI = {
    "MTHFR": ("1", 11685723, 11906455),
    "FURIN": ("15", 90768588, 90983564),
}


def eur_indices():
    """Column offset of each EUR sample in the VCF sample block."""
    samples, pop = [], {}
    with open(os.path.join(HERE, "panel.txt"), encoding="utf-8") as f:
        next(f)
        for ln in f:
            p = ln.split()
            if len(p) >= 7:
                samples.append(p[1])
                pop[p[1]] = p[6]
    keep = [(i, s) for i, s in enumerate(samples) if pop.get(s) == "EUR"]
    return samples, [i for i, _ in keep], [s for _, s in keep]


def main():
    loci = sys.argv[1:] or list(LOCI)
    samples, idx, names = eur_indices()
    print("panel samples: %d | EUR: %d" % (len(samples), len(idx)), flush=True)
    with open(os.path.join(HERE, "eur_samples.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(names) + "\n")

    for locus in loci:
        chrom, beg, end = LOCI[locus]
        t0 = time.time()
        tb = Tabix(BASE % chrom, tbi_path=os.path.join(HERE, "chr%s.vcf.gz.tbi" % chrom))
        rows = tb.query(chrom, beg, end)
        print("%s: %d variants fetched in %.0fs" % (locus, len(rows), time.time() - t0),
              flush=True)

        out = os.path.join(HERE, "gt_%s.tsv" % locus)
        n = 0
        with open(out, "w", encoding="utf-8", newline="\n") as f:
            for ln in rows:
                c = ln.split("\t")
                if len(c) < 10 or "SNV" not in c[2] and len(c[3]) != 1:
                    if len(c[3]) != 1 or len(c[4]) != 1:
                        continue                      # biallelic SNVs only
                pos, ref, alt, info = int(c[1]), c[3], c[4], c[7]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                gts = c[9:9 + len(samples)]
                if len(gts) != len(samples):
                    continue
                doses = []
                for i in idx:
                    g = gts[i][:3]
                    if len(g) < 3:
                        doses.append(-1)
                        continue
                    a, b = g[0], g[2]
                    if a == "." or b == ".":
                        doses.append(-1)
                    else:
                        doses.append((a != "0") + (b != "0"))
                ok = [d for d in doses if d >= 0]
                if len(ok) < 300:
                    continue
                af = sum(ok) / (2.0 * len(ok))
                maf = min(af, 1 - af)
                f.write("%s\t%d\t%s\t%s\t%.4f\t%s\n"
                        % (chrom, pos, ref, alt, maf,
                           "\t".join(str(d) for d in doses)))
                n += 1
        print("%s: %d biallelic SNVs written -> %s" % (locus, n, out), flush=True)


if __name__ == "__main__":
    main()
