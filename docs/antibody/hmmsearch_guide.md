# hmmsearch 完全指南：从序列谱到模板发现

> 面向无生物学背景的开发者，聚焦算法原理与 Protenix 集成实现

> **系列导航**：[mmCIF 格式](mmcif_guide.md) · [生物组装体](bioassembly_guide.md) · [MSA](msa_guide.md) · **hmmsearch** · [模板](template_guide.md) · [预处理](preprocessing_guide.md) · [模型输入](model_input_guide.md) · [损失与优化](loss_and_optimization_guide.md)

---

## 目录

1. [为什么需要 hmmsearch](#1-为什么需要-hmmsearch)
2. [从序列比对到谱比对：搜索方法的演进](#2-从序列比对到谱比对搜索方法的演进)
3. [隐马尔可夫模型（HMM）：核心直觉](#3-隐马尔可夫模型hmm核心直觉)
4. [Profile HMM：蛋白质家族的概率画像](#4-profile-hmm蛋白质家族的概率画像)
5. [从 MSA 构建 Profile：hmmbuild 的工作原理](#5-从-msa-构建-profilehmmbuild-的工作原理)
6. [用 Profile 搜索数据库：hmmsearch 的评分逻辑](#6-用-profile-搜索数据库hmmsearch-的评分逻辑)
7. [统计显著性：E-value 与 bit score](#7-统计显著性e-value-与-bit-score)
8. [HMMER 的加速管线：MSV/Viterbi/Forward 三级过滤](#8-hmmer-的加速管线msvviterbiforward-三级过滤)
9. [数据格式：Stockholm 与 A3M 的转换角色](#9-数据格式stockholm-与-a3m-的转换角色)
10. [Protenix 的 hmmsearch 集成：端到端管线](#10-protenix-的-hmmsearch-集成端到端管线)
11. [配置参数详解](#11-配置参数详解)
12. [输出解析：从 A3M 到 TemplateHit](#12-输出解析从-a3m-到-templatehit)
13. [工程边界与容错设计](#13-工程边界与容错设计)
14. [hmmsearch 在整体管线中的位置](#14-hmmsearch-在整体管线中的位置)
15. [设计哲学总结](#15-设计哲学总结)

---

## 1. 为什么需要 hmmsearch

### 模板搜索的核心问题

结构预测模型（如 Protenix 所基于的 AlphaFold 3 架构）需要回答一个问题：**蛋白质数据库（PDB）中，有没有已经被实验解析过的蛋白质结构，与我要预测的蛋白质"足够相似"，可以作为三维空间参考？**

这些"足够相似"的已知结构称为**模板（template）**。找到它们的过程称为**模板搜索**（template search）。模板搜索的质量直接影响预测精度——好的模板是结构预测中最强的信号来源之一（详见 [模板完全指南](template_guide.md)）。

### 什么叫"足够相似"

蛋白质由氨基酸链组成，每种蛋白质有一条特定的氨基酸序列。进化过程中，蛋白质的**三维结构**比**序列**保守得多——两个蛋白质的序列可能只有 25% 相同，却折叠成几乎一样的形状。

这意味着，模板搜索不能只找"序列高度相同"的蛋白质，还需要能发现那些"序列差异大但结构保守"的远同源蛋白。这正是 hmmsearch 存在的理由。

### 从 MSA 到 hmmsearch 的逻辑链

[MSA 完全指南](msa_guide.md) 中介绍了多序列比对（MSA）——将来自不同物种的同源蛋白质序列对齐排列，得到一个矩阵。这个矩阵携带了蛋白质家族的**进化统计信息**：哪些位置保守、哪些位置允许变化、哪些位置经常一起突变。

hmmsearch 的思路是：**把 MSA 中的统计信息压缩成一个数学模型（Profile HMM），用这个模型去搜索 PDB 数据库中的每一条序列，找到与整个蛋白质家族匹配的模板**。

用一句话概括逻辑链：

```
查询序列 → MSA 搜索 → 同源序列矩阵 → Profile HMM → hmmsearch → PDB 模板候选
```

---

## 2. 从序列比对到谱比对：搜索方法的演进

在理解 hmmsearch 之前，有必要先看看它要解决的问题，以及为什么更简单的方法不够用。

### 第一代：序列-序列比对（BLAST）

最直觉的搜索方式是拿查询序列直接与数据库中的每条序列逐一比对，计算相似度评分，取评分最高的作为候选。BLAST（Basic Local Alignment Search Tool）是这类方法的代表。

**优点**：速度快，实现简单。
**局限**：只能发现与查询序列**直接相似**的目标。如果查询序列与模板序列在进化上已经偏离很远（序列相似度低于 25%），BLAST 基本找不到它们。

### 第二代：谱-序列比对（PSI-BLAST / hmmsearch）

换一个思路：不用一条序列去搜，而是先构建一个**序列谱（Profile）**——一个描述"这个蛋白质家族在每个位置上倾向于出现哪种氨基酸"的统计模型——然后用这个谱去搜索数据库。

| 比较维度 | 序列-序列比对 | 谱-序列比对 |
|---------|------------|-----------|
| **搜索依据** | 一条具体序列 | 一个概率模型（综合数百条同源序列） |
| **灵敏度** | 只发现近同源 | 能发现远同源（序列相似度 < 25%） |
| **背景信息** | 没有家族统计 | 每个位置都有保守性/可变性的概率分布 |
| **类比** | 拿一张照片找人 | 拿一个人脸特征的统计画像找人 |

谱-序列比对的核心优势在于：它不问"这个位置和查询序列一不一样"，而问"这个位置出现这种氨基酸的**概率是多少**"。一个在所有同源序列中都是亮氨酸（L）的位置，如果候选序列也是 L，得分很高；如果候选是异亮氨酸（I，化学性质与 L 相近），得分仍然不低；如果是天冬氨酸（D，带负电荷，性质完全不同），得分就很低。这种**连续的概率评分**比简单的"匹配/不匹配"携带了远更丰富的信息。

hmmsearch 属于谱-序列比对的一种，但它使用的不是简单的频率矩阵，而是一种更强大的数学框架——**Profile Hidden Markov Model（Profile HMM）**。

### 为什么不停留在 PSI-BLAST

PSI-BLAST（Position-Specific Iterated BLAST）也是谱-序列比对，它通过迭代搜索构建 PSSM（Position-Specific Scoring Matrix）。与 Profile HMM 的区别在于：

**PSSM 只记录"每个位置出现每种氨基酸的概率"**，没有显式建模插入和缺失。
**Profile HMM 同时建模匹配、插入、缺失三种事件的概率**，对于长度不同的序列比对能力更强。

在蛋白质进化中，**插入（insertion）和缺失（deletion）**是常见现象——蛋白质环区经常发生长度变化。Profile HMM 能优雅地处理这些事件，而 PSSM 则用固定的 gap 罚分粗暴地惩罚所有插入/缺失，不区分不同位置的差异。

---

## 3. 隐马尔可夫模型（HMM）：核心直觉

### 适合开发者的类比

如果你有 NLP（自然语言处理）背景，HMM 就像一个有限状态机，每个状态会输出一个符号（observable），但你只能观察到输出的符号序列，看不到状态转移路径——"隐"指的就是状态序列不可观测。

更具体地说，HMM 定义了两层结构：

**状态层（隐藏的）**：一系列离散状态，按照**转移概率**在状态之间移动。你可以把它想象成一个有向图，节点是状态，边上标注了从一个状态移动到另一个状态的概率。

**输出层（可观测的）**：每个状态到达时，按照**发射概率**产生一个可观测符号。在蛋白质的场景下，符号就是 20 种氨基酸中的一种。

```
状态：  S₁  →  S₂  →  S₃  →  S₄  →  ...
        ↓       ↓       ↓       ↓
输出：  M       V       L       S        （氨基酸序列）
```

给定一条氨基酸序列（可观测输出），HMM 的任务是计算**这条序列是由该模型生成的概率有多大**。如果概率高，说明这条序列符合模型描述的蛋白质家族特征；如果概率低，说明这条序列不属于该家族。

### HMM 的三个核心问题

HMM 理论中有三个经典问题，在蛋白质搜索中各有对应：

| 问题 | 算法 | 在蛋白质搜索中的含义 |
|------|------|-------------------|
| **评估问题**：给定模型和观测序列，计算序列的生成概率 | Forward 算法 | 这条蛋白质序列与家族模型匹配得多好？ |
| **解码问题**：给定模型和观测序列，找到最可能的状态路径 | Viterbi 算法 | 序列的每个位置最可能对应模型的哪个状态？（即：对齐方式是什么？） |
| **学习问题**：给定观测数据，调整模型参数使概率最大化 | Baum-Welch 算法 | 如何从 MSA 数据中训练出最优的 Profile HMM？ |

hmmbuild（构建 Profile）对应的是学习问题；hmmsearch（搜索数据库）对应的是评估问题和解码问题的结合。

---

## 4. Profile HMM：蛋白质家族的概率画像

### 模型结构

Profile HMM 不是通用 HMM 的随意连接，而是一种**特定的拓扑结构**，专门为序列比对设计。它的状态分为三类，沿序列方向排列成一个链式结构：

**① Match 状态（Mᵢ）**：对应序列的第 i 个对齐位置。每个 Match 状态有自己的发射概率分布——描述该位置上 20 种氨基酸各自出现的概率。保守位置的分布集中（某种氨基酸概率 > 80%），可变位置的分布平坦（多种氨基酸概率相近）。

**② Insert 状态（Iᵢ）**：允许在第 i 个和第 i+1 个对齐位置之间插入额外的氨基酸。Insert 状态也有发射概率，通常接近均匀分布（插入的氨基酸没有特别偏好）。Insert 状态可以自环（连续插入多个氨基酸）。

**③ Delete 状态（Dᵢ）**：跳过第 i 个对齐位置，不发射任何氨基酸。这建模了序列中的缺失——某条序列在该位置没有对应的氨基酸。

```
    M₁ → M₂ → M₃ → M₄ → ...     （匹配链：对齐位置）
   ↗ ↘  ↗ ↘  ↗ ↘  ↗ ↘
  I₀   I₁   I₂   I₃   ...       （插入状态：额外氨基酸）
   ↗ ↘  ↗ ↘  ↗ ↘  ↗ ↘
    D₁ → D₂ → D₃ → D₄ → ...     （删除链：跳过位置）
```

（实际的转移比图示更丰富：M→M、M→I、M→D、I→M、I→I、D→M、D→D 等方向都有转移概率。）

### 一个具体的例子

假设我们有一个 3 个位置的 Profile HMM，描述一个极简蛋白质家族：

```
位置 1（M₁）发射概率：L=0.7, I=0.2, V=0.08, 其他=0.02
位置 2（M₂）发射概率：E=0.9, D=0.08, 其他=0.02
位置 3（M₃）发射概率：均匀分布（每种氨基酸约 0.05）
```

对于序列 "LEG"：
- 位置 1 发射 L：概率 0.7（高分——该位置偏好 L）
- 位置 2 发射 E：概率 0.9（高分——该位置偏好 E）
- 位置 3 发射 G：概率 ≈ 0.05（中等——该位置没有特别偏好）
- 总概率 ≈ 0.7 × 0.9 × 0.05 = 0.0315（还需乘以转移概率）

对于序列 "DWG"：
- 位置 1 发射 D：概率很低（该位置不偏好 D）
- 位置 2 发射 W：概率很低（该位置强烈偏好 E/D）
- 整体概率远低于 "LEG"

这就是 Profile HMM 的判别逻辑：它不是问"这条序列和查询序列像不像"，而是问"这条序列像不像这个家族"。

### 与简单频率矩阵（PSSM）的本质区别

PSSM 只记录了上面的"发射概率"部分。Profile HMM 在此基础上增加了两个维度：

**插入模型**：位置 1 和位置 2 之间可以插入任意多个氨基酸，概率由 I₁ 状态的自环转移概率控制。不同位置的插入概率可以不同——蛋白质环区（loop）经常发生插入，β 折叠核心几乎不发生。

**缺失模型**：可以通过 D₂ 状态跳过位置 2，概率由 M₁→D₂ 的转移概率决定。某些位置在整个家族中几乎从不缺失（M→D 概率极低），另一些位置经常在某些物种中缺失（M→D 概率较高）。

**关键区别**：PSSM 用固定的 gap 开放/延伸罚分处理所有插入缺失；Profile HMM 用**位置特异的**转移概率处理——每个位置的插入和缺失概率独立估计，反映了真实的进化约束。

---

## 5. 从 MSA 构建 Profile：hmmbuild 的工作原理

### 输入：多序列比对

hmmbuild 是 HMMER 套件中负责"学习"的工具。它的输入是一个 MSA——经过对齐的同源序列集合，通常来自之前的 jackhmmer 或 mmseqs2 搜索。

在 Protenix 中，hmmbuild 的输入来源是 `pairing.a3m` 和 `non_pairing.a3m` 中的 MSA 数据。构建流程会先将 A3M 格式转换为 Stockholm 格式（HMMER 的原生输入格式，见第 9 节），再调用 hmmbuild。

### 构建过程的三个阶段

**阶段一：确定模型长度（对齐列选择）**

MSA 矩阵中有些列几乎所有序列都有氨基酸（核心对齐列），有些列大部分序列是 gap（边缘或插入列）。hmmbuild 需要决定哪些列对应 Match 状态，哪些列归入 Insert 状态。

HMMER 有两种模式：
- **`--fast`（默认）**：自动根据 gap 比例判断——gap 比例 < 50% 的列作为 Match 位置，其余作为 Insert 位置
- **`--hand`（手动）**：由输入 MSA 的格式标记决定——Stockholm 格式中可以用 `#=GC RF` 行显式标注哪些列是对齐列（标记 `x`）、哪些是插入列（标记 `.`）

Protenix 在模板搜索中使用 `--hand` 模式（见 `protenix/data/tools/search.py:556`），这意味着 MSA 中 A3M 格式的大小写约定直接决定了 Profile 的结构——大写字母对应 Match 位置，小写字母对应 Insert 位置。

**阶段二：估计发射概率**

对每个 Match 状态 Mᵢ，统计 MSA 中该对齐列上各氨基酸的出现频率。但原始频率直接使用存在一个问题：如果某种氨基酸在该列中从未出现，它的概率为 0，意味着任何在该位置包含该氨基酸的序列都不可能匹配模型——这太严格了。

HMMER 使用**混合先验（mixture Dirichlet priors）**进行平滑：将观测频率与先验分布混合，确保每种氨基酸都有非零概率。先验分布基于蛋白质氨基酸替换的一般统计规律（例如疏水氨基酸之间更容易相互替换），使得即使某种氨基酸未在 MSA 中出现，只要它与观测到的氨基酸化学性质相似，也会得到一个合理的小概率。

**阶段三：估计转移概率**

对 MSA 中相邻列之间的 gap 模式进行统计：
- 多少序列在位置 i 有氨基酸且在位置 i+1 也有 → 贡献给 M→M 转移
- 多少序列在位置 i 有氨基酸但在位置 i+1 是 gap → 贡献给 M→D 转移
- 多少序列在位置 i 之前有插入 → 贡献给 M→I 和 I→I 转移

同样使用先验平滑处理。最终得到每个位置 9 种转移概率（M→M, M→I, M→D, I→M, I→I, D→M, D→D 以及起始/终止转移）。

### Protenix 中的调用

```python
# protenix/data/tools/search.py

class Hmmbuild(BinaryWrapper):
    def build_profile_from_sto(self, sto: str, model: str = "fast") -> str:
        """从 Stockholm 格式 MSA 构建 HMM Profile"""
        return self.build(sto, "stockholm", model)
```

实际命令行等效于：
```bash
hmmbuild --hand --amino output.hmm input.sto
```

输出是一个文本格式的 HMM Profile 文件，包含所有 Match/Insert/Delete 状态的发射和转移概率，以对数几率（log-odds）的形式编码。

---

## 6. 用 Profile 搜索数据库：hmmsearch 的评分逻辑

### 搜索的方向

hmmsearch 的工作方式是：**拿一个 Profile HMM，逐一扫描数据库中的每条序列，计算每条序列与 Profile 的匹配分数。**

在 Protenix 中，搜索的目标数据库是 `pdb_seqres.fasta`——PDB 中所有蛋白质链的序列集合。这个数据库中的每条序列都对应一个已知三维结构，找到高分匹配就意味着找到了潜在的结构模板。

### 评分原理：对数几率比（Log-Odds Ratio）

对于目标序列中的每个氨基酸，评分不是单纯看"概率有多大"，而是看"**这个氨基酸出现在这个位置的概率，相对于随机出现的概率，高了多少倍**"：

```
位置 i 的评分 = log( P(氨基酸 | Match 状态 Mᵢ) / P(氨基酸 | 随机模型) )
```

- 如果该位置强烈偏好该氨基酸（如保守位置出现了保守的残基），分数为正
- 如果该位置对该氨基酸没有特别偏好（概率接近随机），分数接近零
- 如果该位置排斥该氨基酸（如疏水位置出现了亲水残基），分数为负

整条序列的总分是所有位置分数之和，加上转移概率的贡献（走不同路径的对数概率之和）。

### 对齐的产生：Viterbi 解码

hmmsearch 不仅计算总分，还输出**最优对齐**——通过 Viterbi 算法找到概率最大的状态路径（序列的哪个位置对应 Profile 的哪个 Match/Insert/Delete 状态）。这个对齐信息就是后续模板坐标提取的基础。

Viterbi 算法的本质是动态规划。定义 `V[i][k]` 为序列前 i 个氨基酸匹配 Profile 前 k 个位置的最优路径概率，递推关系类似编辑距离（Needleman-Wunsch），但状态空间更丰富。时间复杂度为 O(L × M)，其中 L 是目标序列长度，M 是 Profile 长度。

### 对齐输出格式

hmmsearch 的对齐结果以 Stockholm 格式输出，记录了：
- 查询 Profile 的哪些位置被匹配
- 目标序列的哪些位置被匹配
- 哪些位置是插入（序列有但 Profile 无）
- 哪些位置是缺失（Profile 有但序列无）

Protenix 随后将 Stockholm 对齐结果转换为 A3M 格式（大写=匹配，小写=插入，`-`=缺失），保存为 `hmmsearch.a3m`。

---

## 7. 统计显著性：E-value 与 bit score

### 为什么需要统计显著性

hmmsearch 在数据库中找到了一条高分序列。但这个高分是真的"同源关系"，还是碰巧在随机序列中得到的？

如果数据库有 10 万条序列，即使每条序列与 Profile 完全无关，也会有一些序列因为统计涨落获得偏高的分数。统计显著性的任务就是区分**真实匹配**和**随机巧合**。

### Bit Score

原始的对数几率分数（单位取决于具体的对数底数和模型参数）不方便跨搜索比较。HMMER 将其标准化为 **bit score**：以 2 为底的对数几率单位。

bit score = S 意味着该序列由 Profile 模型生成的概率是由随机模型生成的概率的 2^S 倍。例如：
- bit score = 10：匹配 Profile 的概率是随机的 1024 倍
- bit score = 50：匹配 Profile 的概率是随机的约 10^15 倍
- bit score < 0：随机模型比 Profile 更能解释该序列

### E-value（期望值）

E-value 是最直观的显著性度量：**在当前大小的数据库中，纯粹因为随机原因而获得同等或更高分数的序列数量的数学期望**。

- E-value = 0.001：在该数据库中，随机得到如此高分的概率约为千分之一——这是一个非常可靠的匹配
- E-value = 1：期望有 1 条随机序列达到这个分数——可能是真的，也可能是巧合
- E-value = 100：期望有 100 条随机序列达到这个分数——几乎可以确定是噪声

E-value 与数据库大小成正比：同样的 bit score，在更大的数据库中 E-value 更大（因为有更多机会产生随机高分）。

### Protenix 为什么用宽松的 E-value 阈值

在 Protenix 的配置中（`protenix/data/tools/search.py:512-525`）：

```python
@dataclasses.dataclass(frozen=True, slots=True)
class HmmsearchConfig:
    e_value: float = 100
    inc_e: float = 100
    dom_e: float = 100
    incdom_e: float = 100
```

E-value 阈值设为 100——看起来非常宽松（允许大量低质量命中通过）。这是有意的设计：

**原因一**：hmmsearch 只是候选生成阶段，后续的[预过滤](template_guide.md#5-预过滤廉价筛选的逻辑)（对齐比例、长度、日期隔离等）和模型自身的注意力机制会进一步筛选。宽松的 E-value 确保不遗漏有价值的远同源模板。

**原因二**：蛋白质结构预测对"远同源"的容忍度比传统序列分析高。即使一个模板的序列匹配质量不高，它的整体折叠拓扑仍可能提供有用的空间先验。

**原因三**：E-value 阈值过严格可能导致某些蛋白质一个模板都找不到，而有模板（即使质量一般）几乎总是比完全没有模板好。

---

## 8. HMMER 的加速管线：MSV/Viterbi/Forward 三级过滤

### 计算瓶颈

Profile HMM 的评分需要对每条目标序列运行动态规划，计算复杂度为 O(L × M)。PDB 序列数据库有数十万条序列，如果对每条序列都运行完整的 Forward 和 Viterbi 算法，计算量巨大。

HMMER3（Protenix 使用的版本）通过一个**三级渐进过滤管线**解决这个问题：先用极快但粗糙的方法淘汰大部分不可能匹配的序列，只对少数通过初筛的序列运行昂贵的精确算法。

### F1：MSV 过滤（最快，最粗糙）

**MSV（Multiple Segment Viterbi）** 是一个极简化的模型：忽略所有插入和缺失状态，只保留 Match 状态链，计算一种"只匹配不 gap"的最优评分。

MSV 的优势在于它可以用 SIMD（Single Instruction, Multiple Data）指令高度向量化——同时处理 16 条或更多序列，充分利用 CPU 的并行能力。每条序列只需微秒级时间。

**默认过滤比例**：约保留数据库的 2%（即淘汰 98% 的序列）。在 Protenix 中对应参数 `filter_f1 = 0.1`，即保留 P-value < 0.1 的序列。

### F2：Viterbi 过滤（中等速度，中等精度）

对通过 F1 的序列运行完整的 **Viterbi 算法**（带插入和缺失），但仍然使用 SIMD 加速。Viterbi 考虑了 gap 信息，精度比 MSV 高很多。

**默认过滤比例**：约保留 F1 结果的 10%。在 Protenix 中对应参数 `filter_f2 = 0.1`。

### F3：Forward 过滤（最慢，最精确）

对通过 F2 的序列运行 **Forward 算法**。Forward 算法与 Viterbi 的区别在于：Viterbi 只找概率最大的路径，Forward 对**所有可能路径**的概率求和。Forward 分数比 Viterbi 分数更精确，但也更慢。

**默认过滤比例**：约保留 F2 结果的 10%。在 Protenix 中对应参数 `filter_f3 = 0.1`。

### 过滤管线的整体效果

```
数据库（~30 万条序列）
    │ F1: MSV filter (filter_f1 = 0.1)
    │ 淘汰 ~90%
    ↓
~3 万条序列
    │ F2: Viterbi filter (filter_f2 = 0.1)
    │ 淘汰 ~90%
    ↓
~3000 条序列
    │ F3: Forward filter (filter_f3 = 0.1)
    │ 淘汰 ~90%
    ↓
~300 条序列
    │ 完整的 Forward + 域定义 + 对齐
    ↓
最终命中列表（hmmsearch.a3m）
```

这个管线使得 HMMER3 能在几秒到几分钟内搜索完整个 PDB 序列数据库，同时几乎不损失灵敏度——F1/F2 阶段淘汰的序列，99.9% 以上确实不会通过最终的 Forward 评分。

### `--max` 模式

当 `filter_max = True` 时，HMMER 跳过所有三级过滤，对数据库中的每条序列都运行完整的 Forward 算法。这最大化了灵敏度但速度极慢，通常只在特殊调试场景中使用。Protenix 默认不启用此模式。

---

## 9. 数据格式：Stockholm 与 A3M 的转换角色

### 为什么有两种格式

Protenix 内部使用 **A3M 格式**（见 [MSA 完全指南第 3 节](msa_guide.md#3-a3m-格式msa-的存储方式)），而 HMMER 套件使用 **Stockholm 格式**作为输入输出。两种格式本质上表达相同的信息（对齐后的多序列），只是编码方式不同。

### A3M 格式回顾

A3M 的核心规则：**大写字母和 `-` 对应查询序列的对齐位置，小写字母表示插入**。

```
>query
MVLSEGEW
>hit1
MV-SEGEWak
>hit2
MVLSEaGEW
```

- `hit1` 在位置 3 有缺失（`-`），在末尾有两个插入（`ak`，小写）
- `hit2` 在位置 5-6 之间有一个插入（`a`，小写）

### Stockholm 格式

Stockholm 是一种**所有序列等宽对齐**的格式，插入用额外的列（dash 填充）表示：

```
# STOCKHOLM 1.0

#=GS query_0 DE query sequence
#=GS hit1_1  DE first hit

query_0  MVLSEG--EW
hit1_1   MV-SEGakEW
hit2_2   MVLSEaGEW-

#=GC RF  xxxx.xxxxx
//
```

注意几个关键区别：
- 所有序列的长度（含 gap）完全相同——每一列要么是对齐列，要么是插入列
- `#=GC RF` 行标注哪些列是对齐位置（`x`）、哪些是插入位置（`.`）
- `#=GS` 行携带序列的元数据描述

Stockholm 格式对 HMMER 很重要，因为 `--hand` 模式下 hmmbuild 直接读取 `#=GC RF` 标注来决定 Profile 的 Match/Insert 结构。

### Protenix 中的转换

```python
# protenix/data/tools/common.py

def convert_a3m_to_stockholm(a3m: str, max_seqs: Optional[int] = None) -> str:
    """A3M → Stockholm：展开插入列，生成 RF 标注"""
    ...

def convert_stockholm_to_a3m(sto_io, ...) -> str:
    """Stockholm → A3M：折叠插入列为小写字母"""
    ...
```

在 hmmsearch 管线中，转换发生两次：

1. **输入侧**：A3M（MSA 存储格式） → Stockholm（hmmbuild 输入格式）
2. **输出侧**：Stockholm（hmmsearch 对齐输出） → A3M（Protenix 内部处理格式）

这两次转换是无损的——信息完全保留，只是编码形式的切换。

---

## 10. Protenix 的 hmmsearch 集成：端到端管线

### 完整数据流

从 MSA 文件到模板候选列表，数据依次经过以下步骤：

```
pairing.a3m / non_pairing.a3m
        │
        │ ① 读取 A3M，截取前 N 条序列
        ↓
   A3M 字符串（≤ 300 条序列）
        │
        │ ② convert_a3m_to_stockholm()
        ↓
   Stockholm 字符串
        │
        │ ③ hmmbuild --hand --amino
        ↓
   HMM Profile 文件（文本格式）
        │
        │ ④ hmmsearch --noali -A output.sto profile.hmm pdb_seqres
        ↓
   Stockholm 对齐结果
        │
        │ ⑤ convert_stockholm_to_a3m()
        ↓
   hmmsearch.a3m（A3M 格式的命中列表）
        │
        │ ⑥ HmmsearchA3MParser.parse()
        ↓
   List[TemplateHit]（结构化命中对象）
```

### 各步骤的代码位置

**步骤 ①-⑤：搜索执行**

入口函数在 `protenix/data/tools/search.py:528`：

```python
def run_hmmsearch_with_a3m(
    database_path, hmmsearch_config, max_a3m_query_sequences, a3m
) -> str:
    """
    输入：A3M 格式的 MSA 字符串
    输出：A3M 格式的搜索结果字符串
    """
    cfg = dataclasses.asdict(hmmsearch_config)
    search_bin = cfg.pop("hmmsearch_binary_path")
    build_bin = cfg.pop("hmmbuild_binary_path")
    searcher = Hmmsearch(
        binary_path=search_bin,
        hmmbuild_binary_path=build_bin,
        database_path=database_path,
        **cfg,
    )
    sto = convert_a3m_to_stockholm(a3m, max_a3m_query_sequences)
    return searcher.query_with_sto(sto, model="hand")
```

这个函数将步骤 ①-⑤ 封装为一次调用：A3M 进，A3M 出。内部通过临时文件与 HMMER 命令行工具交互（hmmbuild 和 hmmsearch 都是独立的可执行文件，通过 `subprocess.run` 调用）。

**步骤 ⑥：结果解析**

在 `protenix/data/template/template_parser.py:591`：

```python
class HmmsearchA3MParser:
    @staticmethod
    def parse(query_seq, a3m_str, skip_first=True, max_hits=200):
        """
        输入：查询序列 + hmmsearch.a3m 内容
        输出：TemplateHit 列表
        """
```

### `max_a3m_query_sequences = 300` 的工程考量

构建 Profile 时不使用全部 MSA 序列，而是只取前 300 条。原因：

- Profile 的质量在 MSA 深度超过约 200 条后提升非常缓慢——额外的序列贡献的统计信息边际递减
- hmmbuild 的运行时间与输入序列数成正比，300 条是效率与质量的经验平衡点
- PDB 搜索数据库本身远小于 UniRef100，不需要超高灵敏度的 Profile 也能发现大多数模板

---

## 11. 配置参数详解

Protenix 通过 `HmmsearchConfig` 数据类管理 hmmsearch 的所有配置：

```python
@dataclasses.dataclass(frozen=True, slots=True)
class HmmsearchConfig:
    hmmsearch_binary_path: str      # hmmsearch 可执行文件路径
    hmmbuild_binary_path: str       # hmmbuild 可执行文件路径
    e_value: float = 100            # 全序列 E-value 阈值 (-E)
    inc_e: float = 100              # 全序列包含 E-value 阈值 (--incE)
    dom_e: float = 100              # 单域 E-value 阈值 (--domE)
    incdom_e: float = 100           # 单域包含 E-value 阈值 (--incdomE)
    alphabet: str = "amino"         # 序列类型：蛋白质用 "amino"
    filter_f1: Optional[float] = 0.1  # F1 (MSV) 过滤阈值
    filter_f2: Optional[float] = 0.1  # F2 (Viterbi) 过滤阈值
    filter_f3: Optional[float] = 0.1  # F3 (Forward) 过滤阈值
    filter_max: bool = False          # 是否跳过所有过滤（--max）
```

### E-value 参数族

HMMER 区分两组 E-value 概念：

| 参数 | HMMER 选项 | 含义 |
|------|-----------|------|
| `e_value` | `-E` | **报告阈值**：E-value 低于此值的命中会出现在输出中 |
| `inc_e` | `--incE` | **包含阈值**：E-value 低于此值的命中被标记为"显著"（significant） |
| `dom_e` | `--domE` | **域级报告阈值**：单个域匹配的 E-value 阈值 |
| `incdom_e` | `--incdomE` | **域级包含阈值**：单个域匹配的显著性阈值 |

"域（domain）"的概念：一条蛋白质可能包含多个独立折叠的区域（域），hmmsearch 会分别对每个匹配区域打分。全序列 E-value 是所有域分数的综合，域级 E-value 是单个匹配区域的分数。

### 过滤参数族

`filter_f1`、`filter_f2`、`filter_f3` 分别控制三级过滤管线的严格程度（见第 8 节）。值越小，过滤越严格（淘汰更多序列），搜索越快但可能遗漏远同源命中。

Protenix 的默认值 0.1 是一个相对宽松的设置——保留 P-value < 0.1 的序列通过每级过滤。HMMER 的默认值为 `f1=0.02, f2=0.001, f3=1e-5`，Protenix 选择更宽松的阈值以确保不遗漏潜在模板。

### `--noali` 选项

Protenix 在调用 hmmsearch 时添加了 `--noali` 选项（`protenix/data/tools/search.py:262`），这告诉 HMMER 不在标准输出中打印对齐的文本格式（只通过 `-A` 选项将对齐写入 Stockholm 文件）。这是一个纯粹的性能优化——省略冗长的文本输出可以减少 I/O 开销。

---

## 12. 输出解析：从 A3M 到 TemplateHit

### hmmsearch.a3m 的结构

hmmsearch 的输出（经过 Stockholm→A3M 转换后）是一个多序列文件，每条记录对应一个匹配 Profile 的 PDB 链：

```
>7wnn_A/33-316 [subseq from] mol:protein length:336  3-hydroxyisobutyrate dehydrogenase
--VTIIGLGAMGTALANAFLDAGHSTTVWNRTAARATALAARGAHHAE...
>5ocm_A/5-289 [subseq from] mol:protein length:291  NAD_Gly3P_dh
-DVTVLGLGLMGQALAGAFLKDGHATTVWNRSEGKAGQLAEQGAVLA...
```

每条记录包含：
- **描述行**（`>` 开头）：`{pdb_id}_{chain}/{start}-{end} ... mol:protein length:{total_length} {description}`
- **序列行**：对齐后的氨基酸序列（大写=匹配，小写=插入，`-`=缺失）

### 描述行的正则解析

`HmmsearchA3MParser._parse_description()` 使用正则表达式解析描述行：

```python
pattern = r"^>?([a-z0-9]+)_(\w+)/([0-9]+)-([0-9]+).*protein length:([0-9]+) *(.*)$"
```

各捕获组的含义：

| 组 | 示例值 | 含义 |
|----|-------|------|
| `[1]` | `7wnn` | PDB ID（4 字符） |
| `[2]` | `A` | 链标识符 |
| `[3]` | `33` | 匹配区域在该链序列中的起始位置 |
| `[4]` | `316` | 匹配区域的结束位置 |
| `[5]` | `336` | 该链的完整序列长度 |
| `[6]` | `3-hydroxy...` | 蛋白质描述文本 |

解析结果存入 `HitMetadata` 数据类：

```python
@dataclasses.dataclass(frozen=True)
class HitMetadata:
    pdb_id: str     # "7wnn"
    chain: str      # "A"
    start: int      # 33
    end: int        # 316
    length: int     # 336
    text: str       # 描述文本
```

### 位置索引的构建

对每条命中序列，解析器需要建立**查询序列与命中序列之间的位置映射**。这通过 `_get_indices()` 方法完成：

```python
@staticmethod
def _get_indices(seq: str, start: int) -> List[int]:
    indices = []
    curr = start
    for char in seq:
        if char == "-":
            indices.append(-1)       # gap：无对应位置
        elif char.islower():
            curr += 1               # 插入：跳过，不生成索引
        else:
            indices.append(curr)     # 匹配：记录当前位置
            curr += 1
    return indices
```

对查询序列和命中序列分别运行此方法，得到两个索引数组（`indices_query` 和 `indices_hit`），长度相同，逐位对应。当两个数组的同一位置都不是 -1 时，表示查询的该位置与命中的该位置对齐。

### TemplateHit 的最终结构

```python
@dataclasses.dataclass(frozen=True)
class TemplateHit:
    index: int              # 在搜索结果中的排名（从 1 开始）
    name: str               # "7wnn_A"（PDB ID + 链名）
    aligned_cols: int       # 对齐的列数
    sum_probs: Optional[float]  # HMMER 置信度分数（hmmsearch 不输出，为 None）
    query: str              # 查询序列（原始）
    hit_sequence: str       # 命中序列（全大写化后）
    indices_query: List[int]  # 查询序列位置索引
    indices_hit: List[int]    # 命中序列位置索引
```

`query_to_hit_mapping` 属性提供便捷的位置映射字典：

```python
@functools.cached_property
def query_to_hit_mapping(self) -> Mapping[int, int]:
    """查询位置 → 命中位置的映射（只包含双方都有对齐的位置）"""
    mapping = {}
    for q_idx, h_idx in zip(self.indices_query, self.indices_hit):
        if (q_idx != -1) and (h_idx != -1):
            mapping[q_idx] = h_idx
    return mapping
```

这个映射是后续[坐标提取](template_guide.md#6-坐标提取对齐是关键)的起点——告诉下游流程"查询序列的第 k 个残基应该从模板的第 mapping[k] 个残基提取坐标"。

### `max_hits = 200` 的防护逻辑

解析器有一个硬性上限：最多收集 200 个有效命中。超过此数就停止解析。

这个限制专门应对一类工程问题：**高度保守的蛋白质家族**（如抗体的 Fc 区域）可能在 PDB 中有数万条相似序列，全部解析会使 DataLoader 线程在单个样本上阻塞数分钟。200 是一个安全上界，因为下游特征化最多使用约 20 个模板（通常只用 4 个），200 个候选提供了足够的冗余。

### `mol:protein` 过滤

解析器在处理每条命中前检查描述行是否包含 `mol:protein`——这是 PDB 序列数据库中蛋白质链的标记。不含此标记的条目（如 RNA 链、DNA 链、配体）直接跳过。

---

## 13. 工程边界与容错设计

### 描述行格式的多样性

PDB 是一个有数十年历史的数据库，不同年代、不同实验室提交的数据格式存在差异。hmmsearch 输出的描述行格式取决于输入数据库（`pdb_seqres.fasta`）中各条目的原始 header。实际遇到的非标准格式包括：

| 异常类型 | 示例 | 问题 |
|---------|------|------|
| 长度字段非数字 | `length:VETVIIGH...` | 正则 `[0-9]+` 无法匹配 |
| 缺少空格分隔 | `mol:proteinMVGLD...` | 无法正确识别 `protein length:` 子串 |
| Tab 分隔、大写 ID | `SRR6185503_6370620\t121\t` | `[a-z0-9]+` 不匹配大写字母 |

### 容错策略

Protenix 对解析失败的命中采取**跳过而非中止**的策略：

```python
# 解析描述行时，捕获异常并跳过
try:
    meta = HmmsearchA3MParser._parse_description(h_desc)
except ValueError:
    logger.warning(f"Skipping unparseable hmmsearch hit: {h_desc}")
    continue
```

这个设计选择有清晰的工程理由：
- 一次 hmmsearch 通常返回数十到数百个命中，个别格式异常不应影响其他有效命中
- 格式异常的命中往往来自非标准或低质量的 PDB 条目，跳过它们几乎不影响模板质量
- 日志警告确保问题可追踪，不是"静默丢失"

### 临时文件管理

HMMER 是命令行工具，通过文件 I/O 与 Protenix 交互。`Hmmsearch` 类使用 `tempfile.TemporaryDirectory` 管理临时文件：

```python
def query_with_hmm(self, hmm: str) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        hmm_p, sto_p = f"{tmp}/q.hmm", f"{tmp}/o.sto"
        pathlib.Path(hmm_p).write_text(hmm)
        cmd = [self.path, "--noali", "--cpu", "8"] + self.flags + ["-A", sto_p, hmm_p, self.db]
        run_shell(cmd, "hmmsearch")
        with open(sto_p) as f:
            return convert_stockholm_to_a3m(f, remove_gaps=False, linewidth=60)
```

临时目录在 `with` 块退出后自动清理，即使 hmmsearch 进程失败也不会留下孤立文件。`--cpu 8` 指定 HMMER 使用 8 个线程进行并行计算。

### 错误传播

`run_shell()` 函数在子进程返回非零退出码时抛出 `RuntimeError`，附带 stdout 和 stderr 内容：

```python
def run_shell(cmd, name, **kwargs):
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{name} failed\nStdout: {e.stdout}\nStderr: {e.stderr}") from e
```

这确保了 HMMER 的错误信息不会被 Python 调用栈吞没，而是以清晰的形式暴露给上层。

---

## 14. hmmsearch 在整体管线中的位置

### 上下游关系

hmmsearch 在 Protenix 的数据流水线中扮演**桥梁角色**——连接 MSA（进化信息）和模板（结构信息）两个子系统：

```
┌─────────────────────────────────────────────────────────┐
│                    数据准备阶段                           │
│                                                         │
│  查询序列 ─→ jackhmmer/mmseqs2 ─→ MSA (pairing.a3m,    │
│                                     non_pairing.a3m)    │
│                    │                                    │
│                    │ 取前 300 条                          │
│                    ↓                                    │
│              hmmbuild ─→ Profile HMM                    │
│                    │                                    │
│                    ↓                                    │
│         hmmsearch (pdb_seqres) ─→ hmmsearch.a3m         │
│                                                         │
└─────────────────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────┐
│                   特征提取阶段                            │
│                                                         │
│  hmmsearch.a3m                                          │
│      │                                                  │
│      │ HmmsearchA3MParser                               │
│      ↓                                                  │
│  List[TemplateHit]                                      │
│      │                                                  │
│      │ 预过滤（对齐比例、长度、日期）                       │
│      ↓                                                  │
│  过滤后命中                                              │
│      │                                                  │
│      │ Kalign 重比对 + mmCIF 坐标提取                     │
│      ↓                                                  │
│  ATOM37 坐标 + 掩码                                      │
│      │                                                  │
│      │ 几何特征化（distogram, unit vectors）              │
│      ↓                                                  │
│  模板特征张量 ──→ 配对表示（pair representation）          │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

### 与 HMMER 家族其他工具的关系

HMMER 套件包含多个工具，Protenix 使用了其中三个：

| 工具 | 用途 | 在 Protenix 中的角色 |
|------|------|---------------------|
| **jackhmmer** | 迭代序列搜索（序列→序列数据库） | MSA 构建（搜索 UniRef 等大型序列库） |
| **hmmbuild** | 从 MSA 构建 Profile HMM | 模板搜索的 Profile 构建步骤 |
| **hmmsearch** | 用 Profile 搜索序列数据库 | 模板搜索（搜索 PDB 序列库） |
| **nhmmer** | 核酸序列搜索（类似 jackhmmer 的 DNA/RNA 版本） | RNA MSA 构建 |
| **hmmalign** | 将序列对齐到 Profile | RNA 序列的重比对 |

jackhmmer 和 hmmsearch 的核心区别：
- **jackhmmer**：输入是一条序列，自己迭代构建 Profile 并搜索，目标是**找到所有同源序列**（构建 MSA）
- **hmmsearch**：输入是一个已有的 Profile，单次搜索，目标是**在特定数据库中找到匹配**（发现结构模板）

### 蛋白质 vs RNA：不同的搜索路径

蛋白质和 RNA 的模板搜索使用不同的 HMMER 工具和数据库：

| 维度 | 蛋白质 | RNA |
|------|--------|-----|
| **MSA 搜索工具** | jackhmmer | nhmmer |
| **模板搜索工具** | hmmsearch | nhmmer（直接序列搜索） |
| **搜索数据库** | pdb_seqres.fasta | PDB RNA 序列库 |
| **alphabet 参数** | `amino` | `rna` |
| **是否构建 Profile** | 是（hmmbuild） | 否（nhmmer 自行处理） |

---

## 15. 设计哲学总结

Protenix 对 hmmsearch 的集成体现了几条一贯的工程理念：

**① 谱搜索优于序列搜索**：用 Profile HMM 而非原始查询序列搜索模板数据库，以灵敏度换取更广的模板覆盖范围。这在远同源检测上的优势是决定性的——序列相似度低于 25% 的有价值模板，只有谱搜索才能可靠发现。

**② 宽进严出**：hmmsearch 使用宽松的 E-value 阈值（100），尽可能保留候选；后续通过对齐比例、长度、日期等多条规则逐步过滤。这种"漏斗式"策略确保了灵敏度和精确度的平衡——宁可在初筛阶段多保留一些噪声，也不在早期丢失有价值的信号。

**③ 封装外部工具的代价最小化**：HMMER 是 C 语言编写的高性能工具，Protenix 不重新实现其算法，而是通过 `subprocess` + 临时文件的方式调用。封装层只负责格式转换（A3M ↔ Stockholm）和参数传递，核心计算完全由 HMMER 完成。这是典型的"不重复造轮子"原则。

**④ 格式转换是胶水而非开销**：A3M 和 Stockholm 之间的转换看似冗余，实际上它使得 Protenix 可以在内部统一使用 A3M 格式（紧凑、对查询序列长度归一化），同时仍然利用 HMMER 的原生 Stockholm 接口。转换代价是微不足道的字符串操作，远小于 hmmsearch 本身的计算量。

**⑤ 防御式解析**：PDB 数据的格式多样性是工程现实，解析器不假设数据完美。对无法解析的条目跳过并记录日志，而非让整个流程崩溃。这种容错设计在处理大规模生物信息学数据时是必不可少的。

**⑥ 资源约束的显式管理**：`max_hits=200` 限制解析数量，`max_a3m_query_sequences=300` 限制 Profile 构建的输入规模，`--cpu 8` 限制 HMMER 线程数。每个资源约束都有明确的工程理由，而非任意选择的"魔法数字"。

---

## 附录：关键概念速查

| 概念 | 定义 |
|------|------|
| **Profile HMM** | 隐马尔可夫模型的一种特定结构，用 Match/Insert/Delete 状态链建模蛋白质家族 |
| **hmmbuild** | HMMER 套件中的工具，从 MSA 训练 Profile HMM |
| **hmmsearch** | HMMER 套件中的工具，用 Profile HMM 搜索序列数据库 |
| **发射概率** | Profile 每个位置产生各种氨基酸的概率分布 |
| **转移概率** | 从一个状态到另一个状态的概率（M→M, M→I, M→D 等） |
| **Viterbi 算法** | 找到给定序列最可能经过的状态路径（最优对齐） |
| **Forward 算法** | 计算序列在所有可能路径下的总生成概率 |
| **Bit score** | 标准化的对数几率分数，以 2 为底 |
| **E-value** | 数据库中随机获得同等或更高分数的序列的期望数量 |
| **MSV** | Multiple Segment Viterbi，HMMER 的快速预过滤算法 |
| **Stockholm 格式** | HMMER 的原生序列对齐格式，所有序列等宽对齐 |
| **A3M 格式** | 紧凑的序列对齐格式，小写字母表示插入 |
| **pdb_seqres** | PDB 中所有蛋白质/核酸链的序列数据库 |
| **TemplateHit** | Protenix 中描述一个模板命中的数据结构 |
| **HitMetadata** | 从 hmmsearch 描述行解析出的 PDB ID、链名、位置等元数据 |
| **对齐列（aligned columns）** | 查询序列与命中序列之间实际对齐（非 gap）的位置数 |
| **sum_probs** | HMMER 输出的匹配概率总和，hmmsearch 管线中不可用（为 None） |

---

> 上一篇：[MSA 完全指南](msa_guide.md) · 下一篇：[结构模板（Template）完全指南](template_guide.md)
