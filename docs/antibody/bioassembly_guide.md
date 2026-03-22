# 生物组装体（Biological Assembly）完全指南

> 面向无生物学背景的开发者，结合 Protenix 项目代码深度解析

> **系列导航**：[mmCIF 格式](mmcif_guide.md) · **生物组装体** · [MSA](msa_guide.md) · [模板](template_guide.md) · [预处理](preprocessing_guide.md) · [模型输入](model_input_guide.md) · [损失与优化](loss_and_optimization_guide.md)

---

## 目录

1. [核心问题：蛋白质在细胞中"长什么样"？](#1-核心问题蛋白质在细胞中长什么样)
2. [不对称单元与生物组装体的区别](#2-不对称单元与生物组装体的区别)
3. [mmCIF 中的 Assembly 定义](#3-mmcif-中的-assembly-定义)
4. [对称操作：从单元扩展到组装体](#4-对称操作从单元扩展到组装体)
5. [bioassembly_dict：Protenix 的核心数据结构](#5-bioassembly_dictprotenix-的核心数据结构)
6. [构建流程：get_bioassembly() 全解析](#6-构建流程get_bioassembly-全解析)
7. [链重命名与副本追踪](#7-链重命名与副本追踪)
8. [等价分子识别](#8-等价分子识别)
9. [过滤策略](#9-过滤策略)
10. [采样：链级别与界面级别](#10-采样链级别与界面级别)
11. [完整数据流](#11-完整数据流)
12. [各场景示例对比](#12-各场景示例对比)

---

## 1. 核心问题：蛋白质在细胞中"长什么样"？

本文建立在 [mmCIF 文件完全指南](mmcif_guide.md) 的基础之上，深入解析生物组装体的概念与 Protenix 的实现。

### 晶体学带来的"误解"

科学家解析蛋白质结构时，常用 **X 射线晶体学**方法：让蛋白质形成晶体，再用 X 射线照射晶体，通过衍射图样反推结构。

晶体的本质是**周期性重复排列的分子**，为了节省空间，晶体学家只记录最小重复单元——**不对称单元（Asymmetric Unit）**，而不是整个晶体。

```
晶体（周期性重复）
 ┌──┬──┬──┐
 │AU│AU│AU│
 ├──┼──┼──┤    AU = Asymmetric Unit（不对称单元）
 │AU│AU│AU│    这是 CIF 文件原始存储的内容
 └──┴──┴──┘
```

**关键区分**：不对称单元 ≠ 蛋白质在细胞中的真实形态。

### 生物组装体（Biological Assembly）

生物组装体是蛋白质**在生物体中实际发挥功能的形式**。

最典型的例子：**血红蛋白**
- 不对称单元：可能只有 1 条链（α 链）
- 生物组装体：4 条链（2 条 α 链 + 2 条 β 链）组成的四聚体，才能真正携带氧气

```
不对称单元中存储的            生物组装体（功能形式）
      α链              →    α₁ β₁
                            α₂ β₂
```

更多典型例子：

| 蛋白质 | 不对称单元 | 生物组装体 | 功能 |
|--------|----------|-----------|------|
| 肌红蛋白（102M） | 1 条链 | 1 条链（单体） | 肌肉储氧 |
| 血红蛋白 | 1 条链 | 4 条链（四聚体） | 血液携氧 |
| 抗体（IgG） | 半个分子 | 2 条重链 + 2 条轻链 | 识别抗原 |
| 核糖体大亚基 | 若干链 | 30+ 条链的复合物 | 蛋白质合成 |

---

## 2. 不对称单元与生物组装体的区别

### 概念对比

| 维度 | 不对称单元（Asymmetric Unit） | 生物组装体（Biological Assembly） |
|------|------------------------------|----------------------------------|
| **定义** | 晶体中的最小重复单元 | 蛋白质在生物体中的功能形态 |
| **来源** | 直接来自 `_atom_site`（实验测量） | 由不对称单元经对称操作**生成** |
| **链数** | 通常较少，晶体压缩 | 等于真实功能所需链数 |
| **多样性** | 同一蛋白可能结晶为不同形式 | 通常唯一（标注为 Assembly 1） |
| **Protenix 方法** | `get_asym_unit()` | `get_bioassembly(assembly_id="1")` |

### 为什么 Protenix 使用生物组装体？

Protenix 的目标是**预测蛋白质在生物体中的真实结构**，不是晶体中的压缩表示。使用生物组装体有以下好处：

1. **训练信号更准确**：模型学习的是功能状态，而非人工晶体包装
2. **接口信息完整**：多条链的接触面（interface）包含重要的相互作用信息
3. **对称性约束**：同源链（同一蛋白的多个拷贝）应有相似结构，提供额外监督

---

## 3. mmCIF 中的 Assembly 定义

mmCIF 文件中有三个关键数据块共同定义生物组装体：

### 3.1 `_pdbx_struct_assembly`：组装体概要

```
_pdbx_struct_assembly.id                   1
_pdbx_struct_assembly.details              author_defined_assembly
_pdbx_struct_assembly.method_details       ?
_pdbx_struct_assembly.oligomeric_details   monomeric
_pdbx_struct_assembly.oligomeric_count     1
```

| 字段 | 含义 | 示例值 |
|------|------|-------|
| `id` | 组装体编号（从1开始） | `1`, `2` |
| `details` | 来源说明 | `author_defined_assembly`（作者定义）/ `software_defined_assembly`（软件计算） |
| `oligomeric_details` | 寡聚体描述 | `monomeric`（单体）/ `dimeric`（二聚体）/ `tetrameric`（四聚体） |
| `oligomeric_count` | **聚合物链的数量** | `1`, `2`, `4`, `?`（未知） |

**Protenix 如何使用 `oligomeric_count`：**

```python
# protenix/data/core/parser.py
def num_assembly_polymer_chains(self, assembly_id: str = "1") -> int:
    for _assembly_id, _chain_count in zip(
        self.cif.block["pdbx_struct_assembly"]["id"].as_array(),
        self.cif.block["pdbx_struct_assembly"]["oligomeric_count"].as_array(),
    ):
        if _assembly_id == assembly_id:
            try:
                chain_count += int(_chain_count)  # 直接用这个数
            except ValueError:
                return None  # '?' → 无法确定，跳过此结构
```

若 `oligomeric_count > 1000` 或为 `?`，Protenix 跳过该结构（太大或信息不完整）。

### 3.2 `_pdbx_struct_assembly_gen`：生成指令

告诉程序"哪些链"用"哪个操作"生成"哪个组装体"：

```
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   1
_pdbx_struct_assembly_gen.asym_id_list      A,B,C,D
```

| 字段 | 含义 |
|------|------|
| `assembly_id` | 属于哪个组装体 |
| `oper_expression` | 使用哪个（些）对称操作，如 `1` 或 `(1,2)(3,4)` |
| `asym_id_list` | 要变换的链列表（逗号分隔） |

复杂表达式示例（用于大型对称复合物）：
```
oper_expression   (1-5)      # 操作1到5依次应用
oper_expression   (1,2)(3)   # 先做操作1和2，再做操作3（复合变换）
```

### 3.3 `_pdbx_struct_oper_list`：对称操作定义

定义每个操作的具体变换矩阵（来自102m.cif）：

```
_pdbx_struct_oper_list.id                   1
_pdbx_struct_oper_list.type                 'identity operation'
_pdbx_struct_oper_list.name                 1_555
_pdbx_struct_oper_list.symmetry_operation   x,y,z
_pdbx_struct_oper_list.matrix[1][1]         1.0000000000
_pdbx_struct_oper_list.matrix[1][2]         0.0000000000
_pdbx_struct_oper_list.matrix[1][3]         0.0000000000
_pdbx_struct_oper_list.vector[1]            0.0000000000
_pdbx_struct_oper_list.matrix[2][1]         0.0000000000
_pdbx_struct_oper_list.matrix[2][2]         1.0000000000
_pdbx_struct_oper_list.matrix[2][3]         0.0000000000
_pdbx_struct_oper_list.vector[2]            0.0000000000
_pdbx_struct_oper_list.matrix[3][1]         0.0000000000
_pdbx_struct_oper_list.matrix[3][2]         0.0000000000
_pdbx_struct_oper_list.matrix[3][3]         1.0000000000
_pdbx_struct_oper_list.vector[3]            0.0000000000
```

这是一个**恒等变换**（identity operation）：3×3 旋转矩阵 = 单位矩阵，3D 平移向量 = 零向量，即不做任何旋转/平移，链保持原位。

对称操作的完整形式：

```
新坐标 = 旋转矩阵 × 原坐标 + 平移向量

    ⎡ x' ⎤   ⎡ m11 m12 m13 ⎤   ⎡ x ⎤   ⎡ v1 ⎤
    ⎢ y' ⎥ = ⎢ m21 m22 m23 ⎥ × ⎢ y ⎥ + ⎢ v2 ⎥
    ⎣ z' ⎦   ⎣ m31 m32 m33 ⎦   ⎣ z ⎦   ⎣ v3 ⎦
```

常见操作类型：

| type 值 | 含义 |
|---------|------|
| `identity operation` | 恒等，不动 |
| `crystal symmetry operation` | 晶体空间群对称操作 |
| `point symmetry operation` | 点对称操作（旋转） |
| `helical symmetry operation` | 螺旋对称 |

---

## 4. 对称操作：从单元扩展到组装体

### expand_assembly() 的工作原理

Protenix 使用 Biotite 库的 `expand_assembly()` 方法（封装在 `protenix/data/core/parser.py`，约第 606 行）：

```python
def expand_assembly(self, structure: AtomArray, assembly_id: str = "1") -> AtomArray:
    assembly = struc.AtomArray(0)      # 空的目标组装体
    assembly_1_mask = []               # 记录哪些原子来自 Assembly 1

    for id, op_expr, asym_id_expr in zip(
        assembly_gen_category["assembly_id"].as_array(str),
        assembly_gen_category["oper_expression"].as_array(str),
        assembly_gen_category["asym_id_list"].as_array(str),
    ):
        if id == assembly_id:
            # 1. 解析操作表达式 → 得到操作ID列表
            operations = _parse_operation_expression(op_expr)

            # 2. 选取对应的链
            asym_ids = asym_id_expr.split(",")
            sub_structure = structure[np.isin(structure.label_asym_id, asym_ids)]

            # 3. 应用旋转+平移变换
            sub_assembly = _apply_transformations(sub_structure, transformations, operations)

            # 4. 合并到组装体
            assembly = assembly + sub_assembly

            # 5. 标记是否属于 Assembly 1
            if id == "1":
                assembly_1_mask.extend([True] * len(sub_assembly))

    assembly.set_annotation("assembly_1", np.array(assembly_1_mask))
    return assembly
```

### 扩展过程示意

以一个假想的二聚体为例，不对称单元中只有链 A，通过旋转180°生成链 A 的镜像拷贝：

```
不对称单元（CIF _atom_site）
┌──────────┐
│  链 A    │  原子坐标 (x, y, z)
└──────────┘

_pdbx_struct_assembly_gen:
  assembly_id = 1
  oper_expression = 1,2       # 两个操作
  asym_id_list = A

操作1（identity）：新坐标 = 原坐标        → 链 A  （原位）
操作2（旋转180°）：新坐标 = R × 原坐标    → 链 A  （旋转后）

扩展后的生物组装体
┌──────────┬──────────┐
│  链 A    │  链 A    │  ← 两个拷贝，暂时同名
└──────────┴──────────┘
           ↓ unique_chain_and_add_ids()
┌──────────┬──────────┐
│  链 A    │  链 A.1  │  ← 重命名，后者加后缀
└──────────┴──────────┘
```

### assembly_1 标注的意义

当一个 CIF 文件定义了多个组装体（Assembly 1、Assembly 2 等），Protenix 默认取 **Assembly 1**，并在最终的 `AtomArray` 中为每个原子打上 `assembly_1` 布尔标注。

后续的 `too_many_chains_filter`（链数量过滤器）会优先保留 `assembly_1 == True` 的原子，确保当链数超上限时，Assembly 1 的链被优先保留。

---

## 5. `bioassembly_dict`：Protenix 的核心数据结构

`get_bioassembly()` 返回一个字典，是后续所有处理的基础：

```python
bioassembly_dict = {
    # ── 基础元数据 ──────────────────────────────────────────
    "pdb_id": str,                       # PDB ID（小写），如 "102m"
    "assembly_id": str,                  # 组装体编号，默认 "1"
    "release_date": str,                 # 首次发布日期，如 "1998-04-08"
    "resolution": float,                 # 分辨率（埃），NMR 时为 -1.0

    # ── 序列与类型信息 ─────────────────────────────────────
    "sequences": dict[str, str],         # entity_id → 规范化氨基酸序列
                                         # 例：{"1": "MVLSE...", "2": "ATCG..."}
    "entity_poly_type": dict[str, str],  # entity_id → 聚合物类型
                                         # 例：{"1": "polypeptide(L)", "2": "polyribonucleotide"}

    # ── 结构数据 ───────────────────────────────────────────
    "atom_array": AtomArray | None,      # Biotite 原子数组（完整展开后的组装体）
                                         # None 表示解析失败或被过滤
    "token_array": TokenArray,           # Token 化后的结构表示

    # ── 统计信息 ───────────────────────────────────────────
    "num_assembly_polymer_chains": int,  # 组装体中的聚合物链总数
    "num_prot_chains": int,              # 其中的蛋白质链数量
    "num_tokens": int,                   # Token 总数（中心原子掩码之和）
}
```

### AtomArray 中的关键注释字段

Protenix 在标准 Biotite AtomArray 之上，添加了大量**自定义注释列**：

| 注释字段 | 类型 | 来源步骤 | 含义 |
|---------|------|---------|------|
| `chain_id` | `str` | 原始 / 重命名后 | 链标识符，如 `"A"`, `"A.1"` |
| `label_entity_id` | `str` | mmCIF | 实体 ID |
| `asym_id_int` | `int` | `unique_chain_and_add_ids` | 链的整数编号（0, 1, 2...） |
| `entity_id_int` | `int` | `unique_chain_and_add_ids` | 实体的整数编号 |
| `sym_id_int` | `int` | `unique_chain_and_add_ids` | 同一实体的第几个拷贝（0, 1, 2...） |
| `mol_id` | `int` | `find_equiv_mol_and_assign_ids` | 共价连接分子的唯一 ID |
| `entity_mol_id` | `int` | `find_equiv_mol_and_assign_ids` | 等价分子共享相同值 |
| `mol_atom_index` | `int` | `find_equiv_mol_and_assign_ids` | 原子在分子内的索引 |
| `is_resolved` | `bool` | 解析阶段 | 原子坐标是否被实验确定 |
| `centre_atom_mask` | `bool` | tokenization | 是否为 token 的代表原子（Cα 等） |
| `assembly_1` | `bool` | `expand_assembly` | 是否属于 Assembly 1 |
| `token_mol_type` | `str` | tokenization | token 分子类型 |

---

## 6. 构建流程：`get_bioassembly()` 全解析

`get_bioassembly()` 是 Protenix 中最复杂的方法之一，约跨越 `parser.py` 的 745-913 行。以下是完整的 10 步流程：

```
输入：.cif 或 .cif.gz 文件
                │
    ┌───────────▼───────────┐
    │  Step 1: 链数校验      │  oligomeric_count > 1000 或 ? → 返回空
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 2: 加载不对称单元  │  get_structure() → AtomArray
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 3: 预过滤（8项）  │  见下文详细列表
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 4: 添加注释       │  token_mol_type, centre_atom_mask 等
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 5: 展开组装体     │  expand_assembly(assembly_id="1")
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 6: 坐标归零      │  未解析原子坐标 → (0, 0, 0)
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 7: 链重命名      │  unique_chain_and_add_ids()
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 8: 后过滤（2项）  │  too_many_chains_filter + remove_clashing_chains
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 9: 等价标注      │  find_equiv_mol_and_assign_ids()
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 10: 参考空间编码  │  add_ref_space_uid() + add_ref_info_and_res_perm()
    └───────────┬───────────┘
                │
            bioassembly_dict
```

### Step 3：预过滤（8 项清洗操作）

在展开组装体**之前**，先对不对称单元做清洗：

| 序号 | 操作 | 说明 |
|------|------|------|
| 1 | 移除水分子 | `HOH`/`WAT` 不参与预测 |
| 2 | 移除氢原子 | 氢原子多由软件推算，实验分辨率通常不足 |
| 3 | 移除全未知链 | 所有残基均未知（`UNK`）的链无信息价值 |
| 4 | 移除 Cα 过远链 | 相邻 Cα 距离过大，说明序列/坐标有问题 |
| 5 | 修正精氨酸 | 部分 PDB 文件中精氨酸（ARG）原子顺序有误 |
| 6 | 补全缺失原子 | 添加 `_struct_conn` 连接但缺失的原子 |
| 7 | MSE → MET 转换 | 硒代甲硫氨酸（MSE）是晶体学标记，替换为标准甲硫氨酸 |
| 8 | 移除元素 X | `type_symbol = X` 的虚拟原子（未知元素）无物理意义 |

---

## 7. 链重命名与副本追踪

### 问题背景

展开组装体后，原来 A 链的多个拷贝**都叫 A**，这在代码中会造成歧义。

### unique_chain_and_add_ids() 的解决方案

（`parser.py` 约第 3049 行）

**重命名规则：**
- 第一次出现的链：保持原名（如 `A`）
- 第 N+1 次出现：加上计数后缀（如 `A.1`、`A.2`）

```python
# 示例：原始链序列 [A, B, A, B, C]
# 重命名后：      [A, B, A.1, B.1, C]

chain_counter = Counter()
for chain in original_chains:
    cnt = chain_counter[chain]
    new_name = chain if cnt == 0 else f"{chain}.{cnt}"
    chain_counter[chain] += 1
```

**同时添加三个整数编号：**

```
链序列：  [A,  B,  A.1, B.1, C]
asym_id_int:  [0,  1,   2,   3,  4]  ← 每条链的唯一序号（0 起）
entity_id_int:[0,  1,   0,   1,  2]  ← 同一实体共享相同值
sym_id_int:   [0,  0,   1,   1,  0]  ← 同一实体的第几个拷贝
```

`sym_id_int` 非常重要：它让模型知道"这是第几份拷贝"，从而利用对称性约束进行训练。

---

## 8. 等价分子识别

### 为什么需要识别等价分子？

同一个生物组装体中，可能有多条完全相同的蛋白质链（如血红蛋白的两条 α 链）。识别它们有助于：
- 数据增强（等价链提供相同的训练信号）
- 损失函数设计（对等价链的预测误差做对称平均）
- 约束建模（等价链应有相同结构）

### find_equiv_mol_and_assign_ids() 的工作逻辑

（`parser.py` 约第 2869 行）

```
第一步：打断非典型聚合物链间的键
        （保留正常共价键，断开晶体堆积人工键）

第二步：用共价连接性（covalent connectivity）
        将所有原子划分为独立分子
        → 每个分子获得唯一 mol_id（0, 1, 2, ...）

第三步：按 entity_id 组合分组
        相同 entity_id 组合的分子 = 可能等价

第四步：比较原子名称序列
        若两分子的原子名称序列完全一致
        → 它们是等价分子，共享相同 entity_mol_id

第五步：写入注释
        mol_id：        每个共价连接单元的唯一 ID
        entity_mol_id： 等价分子共享相同值
        mol_atom_index：原子在分子内的顺序位置（0 起）
```

### 示例：同二聚体蛋白

```
组装体：链 A（蛋白X） + 链 A.1（蛋白X，完全相同）

mol_id:         [0, 0, 0, ..., 1, 1, 1, ...]
                  ← 链A原子 →    ← 链A.1原子 →

entity_mol_id:  [5, 5, 5, ..., 5, 5, 5, ...]
                  ← 相同！两个分子等价，共享 entity_mol_id=5 →

mol_atom_index: [0, 1, 2, ..., 0, 1, 2, ...]
                  ← 从头计数 →  ← 从头计数 →
```

---

## 9. 过滤策略

Protenix 有两道主要过滤关卡：

### 9.1 链数量过滤（too_many_chains_filter）

**触发时机：** 展开组装体后（Step 8）

**规则：**
- 最多保留 **20 条**聚合物链
- 最多 **5120 个** token
- 优先保留 `assembly_1 == True` 的链（`core_indices`）

```python
atom_array, removed_chains = Filter.too_many_chains_filter(
    atom_array,
    core_indices=core_indices,   # Assembly 1 的原子优先
    max_chains_num=20,
    max_tokens_num=5120,
)
```

**为什么限制 20 链 / 5120 token？**
Transformer 模型的计算量随序列长度二次方增长。超过这个阈值，GPU 显存会溢出，训练无法进行。

### 9.2 碰撞链过滤（remove_clashing_chains）

**触发时机：** 链数量过滤之后

**规则：** 若一条链有超过 **30%** 的原子与其他链的原子距离 < **1.7 Å**，则删除该链。

1.7 Å 大约是两个非氢原子之间的最小物理距离——比这更近说明坐标有误或晶体包装人工痕迹。

```python
atom_array, removed_chain_ids = Filter.remove_clashing_chains(
    atom_array,
    core_indices=core_indices,
    clashing_threshold=1.7,        # 距离阈值（埃）
    clashing_percentage=0.3,       # 30% 原子碰撞即判定为问题链
)
```

### 过滤前后统计

Protenix 在 `dataset.py` 中记录了过滤统计，可在训练日志中观察：

```
num_assembly_polymer_chains:  解析出的原始链数
num_prot_chains:              保留的蛋白质链数
num_tokens:                   最终 token 数量
atom_array is None:           True 表示整个结构被过滤掉
```

---

## 10. 采样：链级别与界面级别

Protenix 将每个生物组装体拆分为多个**采样单元（Sample）**，用于训练：

### 10.1 链采样（Chain Sample）

每条聚合物链作为一个独立的预测目标：

```python
# make_chain_indices() 生成的记录结构
{
    "type": "chain",
    "chain_1_id": "A",
    "entity_1_id": "1",
    "mol_1_type": "prot",           # 分子类型：prot / nuc / ligand
    "cluster_1_id": "MVLS...",      # 聚类 ID（蛋白质用序列，配体用残基名）

    # 单链时，chain_2 字段为空
    "chain_2_id": "",
    "entity_2_id": "",
    "mol_2_type": "",
    "cluster_2_id": "",

    # 元数据
    "pdb_id": "2lwu",
    "assembly_id": "1",
    "release_date": "2013-08-21",
    "num_tokens": 72,
    "num_prot_chains": 1,
    "resolution": -1.0,             # NMR，无分辨率
    "mol_type_group": "intra_prot", # 用于批次采样分组
    "eval_type": "intra_prot",      # 用于评估分类
}
```

### 10.2 界面采样（Interface Sample）

两条链之间的接触面作为预测目标，要求两链间至少有一对原子距离 < **5 Å**：

```python
# make_interface_indices() 生成的记录结构
{
    "type": "interface",
    "chain_1_id": "A",
    "chain_2_id": "B",
    "mol_1_type": "prot",
    "mol_2_type": "ligand",
    "cluster_id": "MVLS...:HEM",    # 两链 cluster_id 排序后用冒号连接
    "mol_type_group": "prot_ligand",
    "eval_type": "prot_ligand",
    # ... 其余元数据同上
}
```

### 10.3 界面检测算法

使用 Biotite 的 **CellList**（空间哈希网格）加速邻域搜索：

```python
# 构建 5Å 格子
cell_list = struc.CellList(atom_array, cell_size=5)

# 对链 A 的每个已解析原子，查询 5Å 内的所有原子
neighbors = cell_list.get_atoms(chain_A_coords, radius=5)

# 找到邻近原子所在的链
neighbor_chains = np.unique(atom_array.chain_id[neighbors])

# 每对相邻链 → 一个 interface 采样记录
```

### 10.4 mol_type_group 分类

采样记录按分子类型组合分类，用于**加权采样**（确保训练数据多样性）：

| mol_type_group | 含义 |
|----------------|------|
| `intra_prot` | 纯蛋白质（单链或多链接口） |
| `prot_nuc` | 蛋白质-核酸接口 |
| `prot_ligand` | 蛋白质-配体接口 |
| `nuc_nuc` | 核酸-核酸接口 |
| `nuc_ligand` | 核酸-配体接口 |
| `ligand_ligand` | 配体-配体接口 |

---

## 11. 完整数据流

从一个 `.cif.gz` 文件到训练可用的数据，完整路径如下：

```
.cif.gz 文件
    │
    ▼  data_pipeline.py: get_data_from_mmcif()
    │
    ├──► MMCIFParser(mmcif_file)
    │         │
    │         ▼
    │    num_assembly_polymer_chains()  ← 读 _pdbx_struct_assembly
    │         │  > 1000 或 ? → 跳过
    │         ▼
    │    get_structure()                ← 加载 _atom_site
    │         │
    │         ▼
    │    8 项预过滤                     ← 清洗不对称单元
    │         │
    │         ▼
    │    expand_assembly("1")          ← 读 _pdbx_struct_assembly_gen
    │         │                          读 _pdbx_struct_oper_list
    │         │                          应用旋转+平移矩阵
    │         ▼
    │    unique_chain_and_add_ids()    ← A,A → A,A.1；添加 sym_id_int 等
    │         │
    │         ▼
    │    too_many_chains_filter()      ← 最多 20 链 / 5120 token
    │    remove_clashing_chains()     ← 移除碰撞链
    │         │
    │         ▼
    │    find_equiv_mol_and_assign_ids() ← 标记等价分子
    │    add_ref_space_uid()             ← 参考空间编码
    │    add_ref_info_and_res_perm()     ← 参考构象信息
    │         │
    │         ▼
    │    bioassembly_dict              ← 核心数据结构
    │
    ├──► AtomArrayTokenizer(atom_array)
    │         │
    │         ▼
    │    token_array                  ← 蛋白质残基级别，配体原子级别
    │
    └──► make_indices(bioassembly_dict)
              │
              ├──► make_chain_indices()     ← 每条链一条记录
              └──► make_interface_indices() ← 每对相邻链一条记录
                        │
                        ▼
              sample_indices_list          ← 训练采样索引

最终输出（写入 .pkl.gz）：
    (sample_indices_list, bioassembly_dict)
```

---

## 12. 各场景示例对比

### 场景一：单体蛋白（如 102M，肌红蛋白）

```
_pdbx_struct_assembly.oligomeric_count  1
_pdbx_struct_assembly_gen.oper_expression  1  （恒等操作）
_pdbx_struct_assembly_gen.asym_id_list    A,B,C,D

不对称单元：链 A（蛋白）+ 链 B（血红素）+ 链 C（水）+ 链 D（...）
扩展后：同上（恒等操作，不复制）
过滤后：移除水，保留链 A + 链 B

bioassembly_dict:
  num_assembly_polymer_chains = 1
  num_prot_chains = 1
  entity_poly_type = {"1": "polypeptide(L)"}
  sequences = {"1": "MVLSEGEWQLVL..."}

sample_indices_list:
  [{"type": "chain", "chain_1_id": "A", "mol_type": "prot", ...},
   {"type": "interface", "chain_1_id": "A", "chain_2_id": "B", ...}]
           ↑蛋白质单链                    ↑蛋白质-配体接口
```

### 场景二：同二聚体（两条相同链）

```
不对称单元：链 A（蛋白X）
_pdbx_struct_assembly.oligomeric_count  2
_pdbx_struct_assembly_gen.oper_expression  1,2  （恒等 + 旋转180°）

展开后：
  链 A  （原位，assembly_1=True）
  链 A  （旋转后，assembly_1=True）

unique_chain_and_add_ids():
  链 A  → chain_id="A",   asym_id_int=0, sym_id_int=0
  链 A  → chain_id="A.1", asym_id_int=1, sym_id_int=1

find_equiv_mol_and_assign_ids():
  两个分子原子名称序列相同 → entity_mol_id 相同（如都为 3）

sample_indices_list:
  [{"type": "chain", "chain_1_id": "A"},
   {"type": "chain", "chain_1_id": "A.1"},
   {"type": "interface", "chain_1_id": "A", "chain_2_id": "A.1"}]
```

### 场景三：大型复合物（如病毒衣壳，>20 链）

```
_pdbx_struct_assembly.oligomeric_count  60   （60聚体）

展开后：60 条链
→ too_many_chains_filter 触发：
   - core_indices = assembly_1 的原子
   - 保留最多 20 条链，优先 assembly_1
   - 若 token 数超 5120，进一步截断

num_assembly_polymer_chains = 60  （原始值，记录在字典中）
实际 atom_array 中：≤ 20 条链
```

### 场景四：NMR 结构（多模型）

```
_exptl.method  'SOLUTION NMR'
NMR 结构无晶体，因此：
  - 只取第一个模型（_atom_site.pdbx_PDB_model_num == 1）
  - 没有晶体对称操作（或仅有 identity）
  - resolution = -1.0（无分辨率概念）
  - bioassembly_dict["resolution"] = -1.0
```

---

## 附录：关键代码位置速查

| 功能 | 文件 | 大致行号 |
|------|------|---------|
| `get_bioassembly()` 主方法 | `protenix/data/core/parser.py` | 745–913 |
| `num_assembly_polymer_chains()` | `protenix/data/core/parser.py` | 135–158 |
| `expand_assembly()` | `protenix/data/core/parser.py` | 606–684 |
| `unique_chain_and_add_ids()` | `protenix/data/core/parser.py` | 3049–3105 |
| `find_equiv_mol_and_assign_ids()` | `protenix/data/core/parser.py` | 2869–3002 |
| `add_ref_space_uid()` | `protenix/data/core/parser.py` | 2761–2787 |
| `too_many_chains_filter()` | `protenix/data/core/filter.py` | 342–400 |
| `remove_clashing_chains()` | `protenix/data/core/filter.py` | 400–529 |
| `make_chain_indices()` | `protenix/data/core/parser.py` | 1460–1532 |
| `make_interface_indices()` | `protenix/data/core/parser.py` | 1534–1589 |
| 数据流水线入口 | `protenix/data/pipeline/data_pipeline.py` | 45–100 |

---

> 上一篇：[mmCIF / CIF.gz 文件完全指南](mmcif_guide.md) · 下一篇：[多序列比对（MSA）完全指南](msa_guide.md)
