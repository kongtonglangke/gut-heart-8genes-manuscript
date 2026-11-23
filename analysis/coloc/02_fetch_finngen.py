"""
02_fetch_finngen.py —— 8 基因窗口内全部变异 x FinnGen R12 多结局

为什么这样做：
  FinnGen 的 GCS 全量 sumstats 在本网络不可达（TCP 超时），但官方逐变异接口
  https://r12.finngen.fi/api/variant/{chr}-{pos}-{ref}-{alt}
  一次调用就返回该变异在 2,470 个表型上的 beta / sebeta / pval / maf / n_case / n_control / pip。
  因此"按变异抓取"等价于"区域全量抓取"，且只需要抓一遍——结局与组织无关。

关键去重：
  同一基因在不同组织里的变异高度重复，按 (chr,pos,ref,alt) 去重后再请求，
  8 个基因窗口合计约 9 千个唯一变异，而不是 8x6 组织重复的 5 万个。

缓存：
  每个基因一个 jsonl 检查点，中断后重跑自动跳过已完成变异。

用法：
  python 02_fetch_finngen.py                # 全部 8 个基因
  python 02_fetch_finngen.py CST3 FURIN     # 指定基因
"""
import os, sys, csv, json, time, threading, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
SRC = os.path.join(ROOT, "03_data_raw", "gtex_full")
OUT = os.path.join(ROOT, "03_data_raw", "outcome")
CKPT = os.path.join(OUT, "_ckpt")
os.makedirs(OUT, exist_ok=True)
os.makedirs(CKPT, exist_ok=True)

BASE = "https://r12.finngen.fi/api/variant/"

CHR = {"LRP1": "12", "CTSD": "11", "MTHFR": "1", "AUTS2": "7",
       "FURIN": "15", "CST3": "20", "FCER1G": "1", "TBXAS1": "7"}

# 预设结局（冻结于方案 v1 第 4 节 M3）
PHENOS = ["I9_HYPTENS", "I9_HYPTENSHD", "I9_HEARTFAIL", "I9_AF",
          "I9_CARDMYO", "N14_CHRONKIDNEYDIS",              # 主/次结局
          "T2D", "E4_HYPERLIPMIX"]                          # 特异性对照

NTHREAD = 8
lock = threading.Lock()


def vid(chromo, pos, ref, alt):
    """FinnGen 只认连字符格式；冒号格式实测 404。"""
    return "%s-%s-%s-%s" % (chromo, pos, ref, alt)


def collect_variants(gene):
    """把该基因在全部已抽组织里的变异按位点去重。"""
    seen, out = set(), []
    for fn in sorted(os.listdir(SRC)):
        if not (fn.startswith(gene + "_") and fn.endswith(".csv")):
            continue
        with open(os.path.join(SRC, fn), encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                k = (int(r["pos"]), r["ref"], r["alt"])
                if k in seen:
                    continue
                seen.add(k)
                out.append(k)
    return out


def load_ckpt(p):
    done = set()
    if os.path.exists(p):
        with open(p, encoding="utf-8") as fh:
            for ln in fh:
                try:
                    done.add(json.loads(ln)["v"])
                except Exception:
                    pass
    return done


def fetch_one(v, tries=3):
    for k in range(tries):
        try:
            req = urllib.request.Request(BASE + v, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read())
            hit = {x["phenocode"]: x for x in d.get("results", []) if x["phenocode"] in PHENOS}
            return v, hit
        except Exception as e:
            if k == tries - 1:
                return v, {"__err__": str(e)[:80]}
            time.sleep(2 * (k + 1))


def run_gene(gene):
    chromo = CHR[gene]
    vs = [vid(chromo, p, ref, alt) for (p, ref, alt) in collect_variants(gene)]
    ck = os.path.join(CKPT, "%s.jsonl" % gene)
    done = load_ckpt(ck)
    todo = [v for v in vs if v not in done]
    print("%-8s 唯一变异=%d 已完成=%d 待抓=%d" % (gene, len(vs), len(done), len(todo)))
    sys.stdout.flush()

    if todo:
        fh = open(ck, "a", encoding="utf-8")
        n = 0
        with ThreadPoolExecutor(max_workers=NTHREAD) as ex:
            futs = {ex.submit(fetch_one, v): v for v in todo}
            for f in as_completed(futs):
                v, hit = f.result()
                fh.write(json.dumps({"v": v, "r": hit}, ensure_ascii=False) + "\n")
                n += 1
                if n % 200 == 0:
                    fh.flush()
                    print("   %s %d/%d" % (gene, n, len(todo)))
                    sys.stdout.flush()
        fh.close()

    # 汇总写 csv
    rows = []
    with open(ck, encoding="utf-8") as fh:
        for ln in fh:
            try:
                o = json.loads(ln)
            except Exception:
                continue
            v = o["v"]
            for ph, x in o["r"].items():
                if ph == "__err__":
                    continue
                rows.append(dict(gene=gene, variant=v, pheno=ph,
                                 beta=x.get("beta"), se=x.get("sebeta"),
                                 p=x.get("pval"), maf=x.get("maf"),
                                 ncase=x.get("n_case"), ncontrol=x.get("n_control"),
                                 pip=x.get("pip")))
    dst = os.path.join(OUT, "%s_finngen.csv" % gene)
    if rows:
        with open(dst, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    print("%-8s 写入 %d 行 -> %s" % (gene, len(rows), dst))
    sys.stdout.flush()


if __name__ == "__main__":
    genes = sys.argv[1:] or list(CHR.keys())
    t0 = time.time()
    for g in genes:
        run_gene(g)
    print("总耗时 %.1f 分钟" % ((time.time() - t0) / 60))
