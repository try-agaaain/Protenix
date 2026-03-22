# mmCIF / CIF.gz 文件完全指南

> 面向无生物学背景的开发者，结合 Protenix 项目代码讲解

> **系列导航**：**mmCIF 格式** · [生物组装体](bioassembly_guide.md) · [MSA](msa_guide.md) · [hmmsearch](hmmsearch_guide.md) · [模板](template_guide.md) · [预处理](preprocessing_guide.md) · [模型输入](model_input_guide.md) · [损失与优化](loss_and_optimization_guide.md)

---

## 目录

1. [背景：为什么要存储分子结构？](#1-背景为什么要存储分子结构)
2. [PDB 数据库简介](#2-pdb-数据库简介)
3. [CIF 与 mmCIF 的关系](#3-cif-与-mmcif-的关系)
4. [.cif.gz 是什么](#4cifgz-是什么)
5. [mmCIF 文件格式详解](#5-mmcif-文件格式详解)
6. [关键数据块逐段解读](#6-关键数据块逐段解读)
7. [原子坐标：_atom_site 详解](#7-原子坐标_atom_site-详解)
8. [Protenix 如何使用 mmCIF 文件](#8-protenix-如何使用-mmcif-文件)
9. [数据处理全流程](#9-数据处理全流程)
10. [常见生物学概念速查](#10-常见生物学概念速查)
11. [工具与资源](#11-工具与资源)

---

## 1. 背景：为什么要存储分子结构？

### 蛋白质是什么？

蛋白质是生命体内执行各种功能的分子机器。它由一串**氨基酸**（amino acid）按顺序连接而成，这条链会在三维空间中折叠成特定形状。

```
氨基酸序列（1D）：M-V-L-S-P-A-D-K...
         ↓ 折叠
三维空间结构（3D）：每个原子都有 x, y, z 坐标
```

蛋白质的**三维结构**决定了它的功能。例如：
- 抗体（antibody）的形状决定它能识别哪种病原体
- 酶的形状决定它能催化哪种化学反应

### 为什么需要存储结构？

科学家通过 X 射线晶体学、核磁共振（NMR）、冷冻电镜（cryo-EM）等实验手段解析蛋白质结构，解析出每个原子在三维空间中的坐标后，需要一种标准格式保存下来，供全球研究者共享、复用。

---

## 2. PDB 数据库简介

**蛋白质数据银行（Protein Data Bank, PDB）** 是全球最大的生物大分子结构公共数据库，成立于 1971 年。

| 属性 | 说明 |
|------|------|
| 网址 | https://www.rcsb.org |
| 条目数 | 约 220,000+ 结构（截至 2025 年） |
| 内容 | 蛋白质、DNA、RNA 及其复合物的三维坐标 |
| 标识符 | 每条记录有唯一的 **4 位 PDB ID**，如 `102M`、`2LWU`、`7ZZX` |

Protenix 项目广泛使用 PDB 数据，例如 `scripts/msa/data/mmcif/` 目录下存放了多个结构文件：
```
102m.cif  1k1a.cif  5zyh.cif  7zzx.cif  8f9h.cif
```
这些文件名中的 `102m`、`1k1a` 等就是 PDB ID（大小写不敏感）。

---

## 3. CIF 与 mmCIF 的关系

### CIF（Crystallographic Information File）

CIF 是**晶体学信息文件**，最初由国际晶体学联合会（IUCr）于 1991 年制定，用于存储晶体学实验数据。它是一种纯文本、自描述的数据格式。

### mmCIF（macromolecular CIF）

mmCIF 是 CIF 的**宏分子扩展版本**，专为大分子（蛋白质、核酸等）设计。PDB 在 1997 年将其作为新的标准格式，并于 2014 年起将 mmCIF 作为所有新结构提交的**官方格式**（取代旧的 PDB 格式）。

```
CIF（通用晶体学格式）
 └── mmCIF（宏分子扩展）
      └── PDBx/mmCIF（PDB 官方使用的 mmCIF 变体）
```

Protenix 代码中明确声明遵循 `mmcif_pdbx.dic`（PDBx 字典）：
```
_audit_conform.dict_name       mmcif_pdbx.dic
_audit_conform.dict_version    5.397
```

### 与旧 PDB 格式的对比

| 特性 | 旧 PDB 格式（.pdb） | mmCIF 格式（.cif） |
|------|---------------------|-------------------|
| 原子数上限 | ~99,999 个 | 无限制 |
| 链数上限 | 62 条 | 无限制 |
| 数据完整性 | 列宽固定，易截断 | 键值对，完整表达 |
| 可扩展性 | 差 | 好（基于字典） |
| 可读性 | 较直观 | 更结构化 |

---

## 4. .cif.gz 是什么

`.cif.gz` = `.cif` 文件经过 **gzip 压缩**后的结果。

- **为什么压缩？** 一个蛋白质复合物的 CIF 文件可能有几 MB 到几十 MB，压缩后通常缩小 5～10 倍
- **如何读取？** 用标准 gzip 库解压后，内容与普通 `.cif` 文件完全相同

Protenix 代码处理 `.cif.gz` 的方式（`protenix/data/core/parser.py`）：

```python
import gzip
import biotite.structure.io.pdbx as pdbx

def _parse(self, file_path: str):
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as f:
            cif_file = pdbx.CIFFile.read(f)
    else:
        cif_file = pdbx.CIFFile.read(file_path)
```

---

## 5. mmCIF 文件格式详解

### 5.1 基本语法规则

mmCIF 是纯文本格式，基于以下规则：

#### ① 数据块（Data Block）

每个文件以 `data_` 开头，后跟 PDB ID，声明一个数据块：

```
data_102M
```

#### ② 单值键值对

```
_entry.id   102M
_exptl.method   'X-RAY DIFFRACTION'
```

格式：`_category.item   value`
- `_category` 是数据类别（如 `_entry`、`_exptl`）
- `.item` 是该类别下的字段名
- 字符串值如含空格需用单引号包裹

#### ③ 循环表（Loop）

用于表示多行表格数据：

```
loop_
_atom_site.group_PDB        # 列名1
_atom_site.id               # 列名2
_atom_site.type_symbol      # 列名3
...
ATOM  1  N  ...             # 数据行1
ATOM  2  C  ...             # 数据行2
```

`loop_` 关键字后跟列名定义，再跟数据行，`#` 表示注释或段落分隔。

#### ④ 特殊值

| 符号 | 含义 |
|------|------|
| `?` | 未知值（unknown） |
| `.` | 不适用（not applicable） |
| `'...'` | 含空格的字符串 |
| `"..."` | 含特殊字符的字符串 |

### 5.2 文件整体结构

一个完整的 mmCIF 文件按以下顺序组织：

```
data_{PDB_ID}
│
├── 元数据（Metadata）
│   ├── _entry.id                    # PDB ID
│   ├── _audit_conform.*             # 字典版本声明
│   ├── _database_2.*                # 交叉引用数据库
│   └── _pdbx_audit_revision_history.*  # 版本历史
│
├── 实验信息（Experiment）
│   ├── _exptl.*                     # 实验方法
│   ├── _refine.*                    # 精化参数（X-ray）
│   └── _em_3d_reconstruction.*      # 冷冻电镜参数
│
├── 生物实体定义（Entities）
│   ├── _entity.*                    # 实体基本信息
│   ├── _entity_poly.*               # 聚合物实体（蛋白/DNA/RNA）
│   ├── _entity_poly_seq.*           # 聚合物序列
│   └── _chem_comp.*                 # 化学组分描述
│
├── 不对称单元（Asymmetric Unit）
│   ├── _struct_asym.*               # 链定义
│   └── _pdbx_struct_assembly.*      # 生物组装体信息
│
├── 原子坐标（Atomic Coordinates）
│   └── _atom_site.*                 # ← 最核心的数据，见第7节
│
└── 连接信息（Connectivity）
    ├── _struct_conn.*               # 共价键、二硫键等连接
    └── _chem_comp_bond.*            # 化学组分内部键
```

---

## 6. 关键数据块逐段解读

### 6.1 _entry：结构标识

```
_entry.id   102M
```

最简单的字段——该结构在 PDB 中的唯一标识符。Protenix 在解析时直接提取：

```python
# protenix/data/core/parser.py
file_id = mmcif_dict["_entry.id"]
```

### 6.2 _exptl：实验方法

```
_exptl.entry_id    102M
_exptl.method      'X-RAY DIFFRACTION'
```

告诉我们结构是用什么方法解析的。Protenix 支持的方法（`protenix/data/constants.py`）：

```python
CRYSTALLIZATION_METHODS = {
    "X-RAY DIFFRACTION",         # X 射线晶体学（最常见）
    "NEUTRON DIFFRACTION",       # 中子衍射
    "ELECTRON CRYSTALLOGRAPHY",  # 电子晶体学
    "ELECTRON MICROSCOPY",       # 冷冻电镜（cryo-EM）
    "SOLUTION NMR",              # 核磁共振
    "SOLID-STATE NMR",
    "FIBER DIFFRACTION",
}
```

### 6.3 _entity 与 _entity_poly：分子实体

**实体（Entity）** 是 mmCIF 中的核心概念：同一种分子（无论在结构中出现多少次）被视为同一个实体。

```
loop_
_entity.id
_entity.type
_entity.pdbx_description
1  polymer  'MYOGLOBIN'       # 实体1：蛋白质
2  non-polymer  'HEME'        # 实体2：血红素（配体）
3  water  'water'             # 实体3：水分子
```

对于聚合物实体，`_entity_poly` 给出更详细的序列信息：

```
_entity_poly.entity_id                    1
_entity_poly.type                         'polypeptide(L)'
_entity_poly.pdbx_seq_one_letter_code     MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG
```

`pdbx_seq_one_letter_code` 就是蛋白质的**氨基酸序列**，用单字母码表示（M=甲硫氨酸，V=缬氨酸，L=亮氨酸...）。

聚合物类型（Protenix 定义）：

| mmCIF 类型 | Protenix 内部名称 | 说明 |
|-----------|-----------------|------|
| `polypeptide(L)` | `PROTEIN_CHAIN` | 蛋白质（L 型氨基酸，最常见） |
| `polypeptide(D)` | `PROTEIN_CHAIN` | 蛋白质（D 型氨基酸） |
| `polyribonucleotide` | `RNA_CHAIN` | RNA |
| `polydeoxyribonucleotide` | `DNA_CHAIN` | DNA |

### 6.4 _struct_asym：链定义

**链（Chain/Asymmetric Unit）** 表示结构中的一条物理链：

```
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A  1  'MYOGLOBIN CHAIN A'
B  2  'HEME GROUP'
```

同一个实体（`entity_id`）可能对应多条链（比如同一个蛋白的两个拷贝），每条链有唯一的 `id`（如 A、B、C...）。

### 6.5 版本历史

```
loop_
_pdbx_audit_revision_history.ordinal
_pdbx_audit_revision_history.revision_date
1  1998-04-08    # 初次发布
2  2008-03-24    # 修订
...
```

Protenix 会读取最新的 `revision_date` 作为结构的发布日期，用于数据过滤。

---

## 7. 原子坐标：`_atom_site` 详解

这是 mmCIF 文件中**最重要、数据量最大**的部分，记录了每个原子在三维空间中的精确坐标。

### 7.1 字段定义

```
loop_
_atom_site.group_PDB         # ATOM 或 HETATM
_atom_site.id                # 原子序号
_atom_site.type_symbol       # 元素符号（N, C, O, S...）
_atom_site.label_atom_id     # 原子名称（N, CA, C, O, CB...）
_atom_site.label_alt_id      # 替代构象ID（通常为 .）
_atom_site.label_comp_id     # 残基名称（MET, VAL, LEU...）
_atom_site.label_asym_id     # 链ID（mmCIF 编号）
_atom_site.label_entity_id   # 实体ID
_atom_site.label_seq_id      # 残基序号（mmCIF 编号）
_atom_site.pdbx_PDB_ins_code # 插入码（通常为 .）
_atom_site.Cartn_x           # X 坐标（埃，Å）
_atom_site.Cartn_y           # Y 坐标（埃，Å）
_atom_site.Cartn_z           # Z 坐标（埃，Å）
_atom_site.occupancy         # 占有率（0.0~1.0，通常为 1.00）
_atom_site.B_iso_or_equiv    # B 因子（温度因子，反映原子振动程度）
_atom_site.pdbx_formal_charge # 形式电荷
_atom_site.auth_seq_id       # 残基序号（作者编号）
_atom_site.auth_comp_id      # 残基名称（作者编号）
_atom_site.auth_asym_id      # 链ID（作者编号）
_atom_site.auth_atom_id      # 原子名称（作者编号）
_atom_site.pdbx_PDB_model_num # 模型编号（NMR 有多个模型）
```

### 7.2 实际数据示例（来自 102m.cif）

```
ATOM   1    N  N   . MET A 1 1   ? 24.512 8.259   -9.688  1.00 33.83 ? 0   MET A N   1
ATOM   2    C  CA  . MET A 1 1   ? 24.523 9.740   -9.865  1.00 32.90 ? 0   MET A CA  1
ATOM   3    C  C   . MET A 1 1   ? 25.889 10.228  -10.330 1.00 31.90 ? 0   MET A C   1
ATOM   4    O  O   . MET A 1 1   ? 26.886 9.516   -10.198 1.00 32.07 ? 0   MET A O   1
```

逐列解读第1行：
| 字段 | 值 | 含义 |
|------|-----|------|
| group_PDB | `ATOM` | 标准氨基酸原子（HETATM 为非标准/配体） |
| id | `1` | 第1个原子 |
| type_symbol | `N` | 氮元素 |
| label_atom_id | `N` | 主链氮原子 |
| label_comp_id | `MET` | 甲硫氨酸（Methionine） |
| label_asym_id | `A` | A 链 |
| label_seq_id | `1` | 第1个残基 |
| Cartn_x/y/z | `24.512 / 8.259 / -9.688` | 三维坐标，单位：埃（Å，1Å = 0.1纳米） |
| occupancy | `1.00` | 占有率100%（原子确实在此位置） |
| B_iso_or_equiv | `33.83` | B 因子，越大表示原子位置越不确定 |

### 7.3 ATOM vs HETATM

| 标记 | 含义 | 示例 |
|------|------|------|
| `ATOM` | 标准聚合物原子（蛋白/DNA/RNA） | 氨基酸、核苷酸 |
| `HETATM` | 非标准原子（配体、修饰残基、水） | 血红素、ATP、Zn²⁺、H₂O |

### 7.4 双重编号系统（label vs auth）

mmCIF 中每个残基/链有**两套编号**：

| 前缀 | 来源 | 特点 |
|------|------|------|
| `label_` | mmCIF 标准编号 | 连续、无间隙、程序友好 |
| `auth_` | 作者原始编号 | 可能跳号、有插入码，与文献一致 |

Protenix 内部主要使用 `auth_` 编号与 PDB 文献对应，同时保留 `label_` 用于内部索引。

---

## 8. Protenix 如何使用 mmCIF 文件

### 8.1 核心解析模块

```
protenix/data/core/parser.py       # 主要解析器（两个类）
protenix/data/core/ccd.py          # 化学组分字典（CCD）解析
protenix/data/core/featurizer.py   # 结构特征提取
protenix/data/pipeline/data_pipeline.py  # 流水线入口
```

### 8.2 两个解析器

**① MMCIFParser（基于 Biotite 库）**

用于训练数据准备，将整个结构解析为 `AtomArray`（原子数组）：

```python
# 调用方式（data_pipeline.py）
parser = MMCIFParser(mmcif_file_path)
atom_array = parser.get_structure()
```

**② TemplateParser（基于 Biopython）**

用于模板搜索结果的解析，返回 `MmcifObject`：

```python
@dataclass(frozen=True)
class MmcifObject:
    file_id: str              # PDB ID
    header: PdbHeader         # 元数据（方法、分辨率、发布日期）
    structure: Structure      # Biopython 结构对象
    chain_to_seqres: dict     # 链ID → 氨基酸序列
    seqres_to_structure: dict # 序列位置 → 结构位置映射
    raw_string: str           # 原始文件内容
```

关于生物组装体的详细逻辑，请参阅 [生物组装体完全指南](bioassembly_guide.md)。

### 8.3 解析流程（以训练数据为例）

```
.cif 或 .cif.gz 文件
        │
        ▼
  MMCIFParser._parse()
  ┌─────────────────────────────────────────┐
  │  1. 处理压缩文件（gzip.open if .gz）      │
  │  2. 用 Biotite 读取 CIF 结构             │
  │  3. 提取 bioassembly（生物学活性组装体）   │
  │  4. 过滤水分子、非标准残基               │
  └─────────────────────────────────────────┘
        │
        ▼
  AtomArray（Biotite 原子数组）
  ┌──────────────────────────────────────────────┐
  │ 每行=一个原子，列=属性：                       │
  │ element, atom_name, res_name, chain_id,       │
  │ res_id, coord (x,y,z), b_factor, occupancy   │
  └──────────────────────────────────────────────┘
        │
        ▼
  AtomArrayTokenizer
  ┌──────────────────────────────────────────────┐
  │  将原子粒度转换为 Token 粒度                   │
  │  蛋白质：每个残基=1个token（以Cα原子为代表）   │
  │  核酸：每个碱基=1个token                      │
  │  配体：每个原子=1个token                      │
  └──────────────────────────────────────────────┘
        │
        ▼
  TokenArray + 特征字典（pickle 文件）
  用于模型训练
```

### 8.4 提取的关键特征

| 特征 | 来源字段 | 用途 |
|------|---------|------|
| 原子坐标 | `_atom_site.Cartn_x/y/z` | 结构真值（监督信号） |
| 残基类型 | `_atom_site.label_comp_id` | 序列信息 |
| 链信息 | `_atom_site.auth_asym_id` | 多链结构区分 |
| 实体类型 | `_entity_poly.type` | 区分蛋白/DNA/RNA/配体 |
| 氨基酸序列 | `_entity_poly.pdbx_seq_one_letter_code` | MSA 输入 |
| 分辨率 | `_refine.ls_d_res_high` | 数据质量过滤 |
| 发布日期 | `_pdbx_audit_revision_history` | 时序过滤（防数据泄露） |

### 8.5 CCD：化学组分字典

除了结构文件，Protenix 还使用**化学组分字典（Chemical Component Dictionary, CCD）**，这是 PDB 官方维护的小分子/配体数据库：

```python
# protenix/data/core/ccd.py
def get_component_atom_array(ccd_code: str, keep_leaving_atoms=False, keep_hydrogens=False):
    """
    根据 CCD 代码获取标准化的配体/修饰残基原子信息
    例如：ccd_code="HEM" → 血红素的标准原子坐标
    """
```

CCD 以 `.cif` 格式存储每种化学组分的理想几何结构，供 Protenix 在处理非标准残基时参考。

---

## 9. 数据处理全流程

### 9.1 训练数据准备

```
PDB 数据库（约 22 万个结构）
        │
        ▼ scripts/prepare_training_data.py
  读取 .cif/.cif.gz
        │
        ├── MMCIFParser → bioassembly_dict
        │     ├── pdb_id, assembly_id
        │     ├── sequences（实体序列字典）
        │     ├── atom_array（原子坐标）
        │     ├── token_array（token 化后的表示）
        │     ├── resolution（分辨率）
        │     └── release_date（发布日期）
        │
        ├── MSA Pipeline（多序列比对）
        │     └── 找到同源蛋白序列，提供进化信息
        │
        └── Template Pipeline（模板搜索）
              ├── hmmsearch 搜索结构模板
              ├── 解析匹配到的 .cif 文件
              └── 提取模板结构特征
        │
        ▼
  .pkl.gz 文件（压缩的 pickle）
  用于模型训练
```

### 9.2 推理时的 CIF 文件

推理时，Protenix 不直接读取 CIF，而是读取 JSON 格式的输入（见 `docs/infer_json_format.md`）。但在模板搜索步骤中，系统会：

1. 通过 hmmsearch 找到最相似的已知结构
2. 从本地缓存或 PDBe API 下载对应的 `.cif` 文件
3. 提取模板坐标作为结构预测的参考

相关代码（`protenix/data/template/template_utils.py`）：

```python
def _fetch_or_read_cif(self, pdb_id: str) -> str:
    """
    先查本地缓存（mmcif/目录），未命中则从 PDBe API 下载
    """
    local_path = os.path.join(self.mmcif_dir, f"{pdb_id}.cif")
    if os.path.exists(local_path):
        return open(local_path).read()
    else:
        # 从 PDBe 下载
        url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
        return requests.get(url).text
```

---

## 10. 常见生物学概念速查

| 概念 | 简单解释 | mmCIF 中的字段 |
|------|---------|--------------|
| **氨基酸（Amino Acid）** | 蛋白质的基本构件，共20种标准类型 | `label_comp_id`（三字母码：MET、VAL...） |
| **残基（Residue）** | 蛋白质链中的一个氨基酸单元 | `label_seq_id` 编号 |
| **链（Chain）** | 一条完整的聚合物分子链 | `label_asym_id` / `auth_asym_id` |
| **实体（Entity）** | 同种分子的抽象（不管有几个拷贝） | `label_entity_id` |
| **Cα（Alpha Carbon）** | 氨基酸主链上的碳原子，常用代表该残基的空间位置 | `label_atom_id = CA` |
| **配体（Ligand）** | 与蛋白质结合的小分子（如药物、金属离子） | `group_PDB = HETATM` |
| **分辨率（Resolution）** | X 射线衍射中的精度指标，单位埃（Å），越小越精确 | `_refine.ls_d_res_high` |
| **B 因子（B-factor）** | 原子位置的不确定性/振动程度，越大越"模糊" | `B_iso_or_equiv` |
| **占有率（Occupancy）** | 原子实际处于该位置的概率（0.0~1.0） | `occupancy` |
| **生物组装体（Bioassembly）** | 蛋白质在生物体中实际发挥功能的四级结构 | `_pdbx_struct_assembly.*` |
| **MSA** | 多序列比对，通过找同源序列提供进化保守性信息 | 外部计算，非 CIF 字段 |
| **插入码（Insertion Code）** | 历史编号问题留下的字母后缀（如残基"47A"） | `pdbx_PDB_ins_code` |

### 20 种标准氨基酸（单字母码 → 三字母码）

```
A=ALA  C=CYS  D=ASP  E=GLU  F=PHE
G=GLY  H=HIS  I=ILE  K=LYS  L=LEU
M=MET  N=ASN  P=PRO  Q=GLN  R=ARG
S=SER  T=THR  V=VAL  W=TRP  Y=TYR
```

在 CIF 文件中，`label_comp_id` 使用三字母码；在 `_entity_poly.pdbx_seq_one_letter_code` 中使用单字母码。

---

## 11. 工具与资源

### 可视化工具

| 工具 | 平台 | 说明 |
|------|------|------|
| [RCSB PDB 3D Viewer](https://www.rcsb.org) | 浏览器 | 在线查看任意 PDB 结构 |
| [PyMOL](https://pymol.org) | 桌面 | 专业分子可视化，可直接打开 CIF |
| [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) | 桌面 | UCSF 出品，支持 CIF |
| [Mol\*](https://molstar.org) | 浏览器 | 现代化 3D 分子查看器 |

### Python 库

| 库 | 用途 | Protenix 使用情况 |
|----|------|-----------------|
| [Biotite](https://www.biotite-python.org) | 结构解析与操作 | **主要使用**（`MMCIFParser`） |
| [Biopython](https://biopython.org) | 结构解析 | 用于模板解析（`TemplateParser`） |
| [gemmi](https://gemmi.readthedocs.io) | 快速 CIF 解析 | 部分使用 |

### 数据来源

| 资源 | 说明 |
|------|------|
| [RCSB PDB](https://www.rcsb.org) | 美国 PDB 主库，下载 .cif.gz |
| [PDBe](https://www.ebi.ac.uk/pdbe) | 欧洲 PDB 镜像，Protenix 模板下载来源 |
| [CCD](https://www.wwpdb.org/data/ccd) | 化学组分字典，Protenix 配体处理依据 |

### 快速查询 PDB 结构

```bash
# 下载指定 PDB ID 的 CIF 文件（以 102m 为例）
curl -O https://files.rcsb.org/download/102m.cif.gz

# 解压查看
gunzip -k 102m.cif.gz
head -50 102m.cif
```

---

## 附录：Protenix 项目 mmCIF 相关文件索引

| 文件路径 | 作用 |
|---------|------|
| `protenix/data/core/parser.py` | 核心解析器（`MMCIFParser` + `TemplateParser`） |
| `protenix/data/core/ccd.py` | 化学组分字典解析 |
| `protenix/data/core/featurizer.py` | 结构→特征转换 |
| `protenix/data/pipeline/data_pipeline.py` | 训练数据流水线入口 |
| `protenix/data/template/template_utils.py` | 模板 CIF 获取与处理 |
| `protenix/data/template/template_parser.py` | 模板搜索结果解析 |
| `protenix/data/constants.py` | 氨基酸编码、聚合物类型等常量 |
| `protenix/utils/file_io.py` | `save_structure_cif()` 保存结构为 CIF |
| `scripts/prepare_training_data.py` | 批量处理训练数据脚本 |
| `scripts/msa/data/mmcif/` | 示例 CIF 文件（102m、1k1a 等） |
| `examples/2lwu.cif` | 推理示例 CIF 文件 |

---

> 下一篇：[生物组装体（Biological Assembly）完全指南](bioassembly_guide.md)
