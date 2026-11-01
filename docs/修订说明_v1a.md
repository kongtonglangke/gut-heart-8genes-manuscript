# 审稿意见核验与修订记录 · Manuscript_GutHeart_8genes_v1（v1a）

生成日期：2026-11-01
核验依据：原文逐字检索 + 项目原始数据重算 + PubMed/代码核对 + 图件目视比对
产出版本：`Manuscript_GutHeart_8genes_v1.docx`（26 项自检连续三轮全通过）
原稿备份：`Manuscript_GutHeart_8genes_v1.src.md`（只读，SHA256 ab7e6243dea77789）

---

## 一、经核实**不成立**的意见（未改动，理由是原稿本就正确）

| # | 审稿意见 | 核验结果 | 证据 |
|---|---|---|---|
| 1 | "MTHFR" 多处误拼为 "MTHFER" | **不成立** | 全文 `MTHFER` = **0 处**；`MTHFR` = 36 处，拼写全部正确 |
| 2 | Abstract 中 "gt–heart axis" 缺字母 | **不成立** | `gt–heart` = **0 处**；`gut–heart` = 12 处，全部正确 |
| 3 | "GTEx v8" 与 "GTEx release 8" 混用 | **不成立** | 全文只有 `GTEx v8`（6 处）+ `the GTEx v8 release`（1 处，同义）；`release 8` = 0 处。`release 12` 专指 FinnGen，未混用 |
| 4 | 参考文献 [1] 标注 2026-06-22 属"未来日期" | **不成立** | PubMed 实核 **PMID 42406622**、DOI **10.3791/71252**，Bai J 等，*J Vis Exp* 2026;(232):e71252。稿件当前日期 2026-11-01，该日期**已过去 4 个多月**。此为审阅模型知识截止于 2025 年导致的误判 |
| 5 | Figure 1 配色 "white → orange → dark red" 需核对 | **不成立** | 目视配色条：浅白（≈2.5）→ 橙 → 深红（10.0），与图注一致 |
| 6 | Figure 2 黑框是否与显著结果一一对应 | **不成立** | 目视 + 数值：10 个黑框值分别为 8/10/10/13/14/14/20/45/48，全部 ≥ 8；未框单元格最大 6.9。eGene 阈值 −log10(5×10⁻⁸) = **7.30** → 一一对应，无漏无多 |
| 7 | Figure 3 是否与 Results 完全对应 | **不成立** | 原始 CSV 重算：240 测试中 **恰好 4 条双侧显著**，全部为 *MTHFR*×高血压，PPH3 = 0.9987 / 0.9994 / 0.99999 / 0.99999（"趋近 1"），PPH4 = 1.5×10⁻⁶–5.9×10⁻⁴。0.80 双虚线目视确认 |
| 8 | Figure 5 rs6227 竖线位置（chr15:90,882,002） | **不成立** | 图中竖线位于 ≈90.882 Mb，与 GRCh38 chr15:90,882,002 一致；c 面板点线为全基因组显著线 |
| 9 | Figure 6b 组织内点分布是否清晰 | **不成立** | 图可读；且已按数据重算确认：295 个变异中 **0 个**跨组织符号翻转，故同 x 位置的点仅按组织区分，表述成立 |
| 10 | Methods 未说明 eQTL Catalogue 统一处理 | **不成立（已存在）** | Methods 原文："obtained from the GTEx v8 release **as uniformly reprocessed by the eQTL Catalogue [3,4]**" |
| 11 | 参考文献 [2]–[11] 格式 | 认可 | 无需改动 |

---

## 二、经核实**成立并已修订**的意见

### 1. Table 1 表注缺 Tau 分档阈值（审稿人指出的真问题）
- 原表注仅写"higher values indicate stronger tissue restriction"，`Class` 列的 intermediate / broadly expressed / tissue-specific 无从自解释。
- 已补：**Tau ≥ 0.80 组织特异、Tau ≤ 0.60 广泛表达、其间为 intermediate**（与 Methods 第 72 行完全一致）。
- 复核：8 个基因的 `tau_class` 与阈值逐条自洽（LRP1 0.773 / CTSD 0.714 / AUTS2 0.645 / FURIN 0.717 = intermediate；MTHFR 0.580 / CST3 0.548 = broadly；FCER1G 0.972 / TBXAS1 0.941 = specific）。

### 2. Table 2 表注称"Bold entries"，但 Word 版无加粗（**审稿人指出+此前已知**）
- 根因：上一版按"永不加粗"规则把表格里的 `**` 标记一并剥离，导致表注指向落空。
- 已改为：**表内允许加粗（正文仍严格不加粗）**。依据是表格内的 `**` 是表注所依赖的语义标记（eGene 状态、合计行），不属于装饰性强调。
- 新增机器断言 `2c`：逐格解析 P 值，**加粗 ⇔ P < 5×10⁻⁸**。结果：10 处加粗、10 格显著，零错配；且与 Figure 2 的 10 个黑框逐格一致。

### 3. Figure 6b 图注「10 kb 内反号」与正文三个全负中位数并置易误读
- 数据重算：377 条 Wald 比值中 **42 条为正**，距高血压 lead 最近的正值变异为 rs35983202（**6.94 kb**）。中位数分箱 −0.139（0–20 kb, n=147）/ −0.141（20–50 kb, n=169）/ −0.063（>50 kb, n=61）**与原文完全一致，原文数字无误**。
- 问题在于 Figure 6b「10 kb 内反号」所指的具体变异对，正文未给出两点间距离，读者无法把两者接上。
- 已改：正文补 **"9.5 kb apart"**（rs12404124 chr1:11,813,812 与 rs198397 chr1:11,823,275，相距 9,463 bp）并补 **"42 of the 377…nearest 6.9 kb"**；图注改写为"中位随距离漂移（三档数值）+ 个别比值符号不守恒（该 9.5 kb 变异对）"，两个层面明确分离。

### 4. Table 5 表注应说明单工具变量法的代价
- 已补：单工具法避免 LD 参考面板不匹配带来的方差膨胀，但**牺牲统计功效与跨工具异质性/多效性检验能力**。同时删除表后重复的 "Wald ratios are single-instrument estimates." 一句。

### 5. Supplementary Figure S1 缺 HPA 报告阈值数值
- 已补：**0.1 nCPM**（矩阵 1232 个值中非零最小值恰为 0.100000，57 个 0 值即 HPA 返回的下限以下值）。
- 同时澄清"完整矩阵、未加阈值"与"零值显示为 absent"并非矛盾。

---

## 三、**两位审稿人都漏掉、但确实存在**的错误（本次主动查出并修正）

### S1 图注称 *CST3* "is highest in enterocytes" —— 事实错误
- 数据：*CST3* 全身上限为**视网膜色素上皮细胞 11,782.8 nCPM**；enterocytes 为 3,917.2 nCPM。原文把"肠道内最高"写成了"全身最高"。
- 已改为 "highest **among intestinal cell types** in enterocytes…although its whole-body maximum lies outside the intestine"。
- 复核：肠道 7 类细胞中 CST3 确以 enterocytes 最高（3,917.2 > colonocytes 2,002.8 > goblet 1,209.8 > tuft 945.3 > TA 910.6 > paneth 653.6 > stem 648.9）✓

### S1 图注称 *TBXAS1* "near the detection floor in intestinal epithelium" —— 夸大
- 数据：*TBXAS1* 肠道内最高为 tuft cells **80.7 nCPM**（其余 1.0–16.3），粒细胞峰值 2,238.6 nCPM（比值 3.6%）。80.7 nCPM 远高于检测下限 0.1 nCPM，谈不上"贴近下限"。
- 已改为分别陈述：两者均以中性粒细胞为峰值，*FCER1G* 在肠道上皮确实基本缺失（0.8 nCPM）；不再对 *TBXAS1* 作"贴近下限"的定性。

### 图注中的反引号会原样进入 Word
- 原 md 用 `` `file.png` `` 标注图件文件名，Word 版会带字面反引号。已在排版层剥离；同时把 `+/-` 统一为 `±`（与 Methods 用法一致，共 5 处）。

---

## 四、第 2 部分（生信提升建议）的处置

这些属于**新增分析**，不是文字问题，无法在改稿中"补上"。逐条判断：

| 建议 | 判断 | 处置 |
|---|---|---|
| ① pQTL 共定位（FURIN/CST3/MTHFR） | **最有价值，且确实可行**（deCODE / UKB-PPP / Fenland 均有 FURIN、CST3、MTHFR 的 cis-pQTL 汇总数据）。Discussion 已自陈该缺口 | **未做**（需新数据+新计算）。Limitations 中原句保留。**建议作为下一篇或本稿补做的第一优先** |
| ② 单细胞 eQTL / 细胞类型反卷积 | 可行但工作量大（OneK1K / sc-eQTLGen 无肠道 tuft 细胞 eQTL） | 未做；Limitations 第 4 条已明确该局限 |
| ③ fine-mapping + 条件分析 | **部分可行**。本机无 LD 面板，coloc.susie 在既往实测中不收敛（需本地 1000G LD） | **已在 Limitations 第 4 条明确写出该局限及其后果**，不再暗示已具备 |
| ④ Trans-eQTL / WGCNA 网络 | 超出现有数据范围 | 未采纳 |
| ⑤ 双向微生物 MR | 已有 MiBioGen 侧工具，但宿主基因→菌群方向样本量极小 | 未采纳 |
| ⑥ 敏感性分析：**MR-Egger / 加权中位数 / 留一法** | **方法学上不可行**——本稿每个基因-组织对只有 **1 个**合格工具变量，Egger/加权中位数/留一法均需 ≥2 个工具。审稿人此处建议有误 | **不采纳**，Limitations 第 2 条已如实说明"single-instrument Wald ratios without heterogeneity or pleiotropy diagnostics"。改为提出可行替代：**p12 先验敏感性 + eCAVIAR/SMR-HEIDI 正交验证** |
| ⑦ RegulomeDB / CADD / LINSIGHT 功能注释 | 可行但需引入新注释流程并新增结果段 | 未做；如需可作为下一次单独任务 |
| ⑧ 扩展组织与终点（肝/脂肪、脉压、LV 质量） | FinnGen R12 有脉压等表型，但工作量等同重做一轮 | 未采纳 |
| ⑨ tuft 细胞功能解读 | **文本层面已具备**（Discussion 已写 tuft 细胞为化学感应哨兵 + 必须以 FES 作对照） | 无需改 |
| ⑩ 多层证据整合评分（随机森林） | 8 个基因、样本极少，过拟合风险高于信息量 | 未采纳，理由充分 |

**结论**：第 2 部分没有可直接落地的文字修订；真正能提升的是 ①pQTL 共定位。若要做，建议单独立项，用与本稿完全相同的 coloc 流程，把 pQTL 作为暴露、FinnGen R12 高血压为结局，至少覆盖 *FURIN*、*CST3*、*MTHFR*。

---

## 五、修订后自检（`03_verify.py`，26 项）

连续三轮：**pass 26 / fail 0 / warn 0**

| 类别 | 结果 |
|---|---|
| 空段落 | 0 |
| 正文加粗 | 0 |
| 表内加粗 | 15 处，与源标记逐格对齐；Table 2 加粗 ⇔ P<5×10⁻⁸ 零错配 |
| 一级/二级标题 | 8 / 31，字号与间距合规 |
| 正文段落 | 75，12 pt / 宋体 / 1.5 倍行距 / 两端对齐 |
| 斜体 | 220 处，仅统计符号 + 基因名 + 期刊名；OR/CI/RR/HR/SD/SE/AUC **零斜体** |
| 三线表 | 6 张全部 `(True, [])` |
| 上标 | 100 处真 `vertAlign` |
| 中文字符 | 0 |
| 文字保真 | md 与 docx 去标记后逐字符一致 |
| 反引号 / `+/-` / 裸 `*` | 0 / 0 / 0 |

字数：Front matter 461 / Introduction 656 / Methods 1,182 / Results 1,399 / Discussion 1,355 / 图注 1,121 / 表格 1,061 / 参考文献 287，合计 **7,522**。
