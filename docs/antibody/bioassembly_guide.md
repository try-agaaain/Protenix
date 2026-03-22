# 生物组装体（Biological Assembly）完全指南

> 面向无生物学背景的开发者，结合 Protenix 项目代码深度解析

> **系列导航**：[mmCIF 格式](mmcif_guide.md) · **生物组装体** · [MSA](msa_guide.md) · [模板](template_guide.md) · [预处理](preprocessing_guide.md) · [模型输入](model_input_guide.md) · [损失与优化](loss_and_optimization_guide.md)

---

## 目录

1. [核心问题：蛋白质在细胞中"长什么样"？](#1-核心问题蛋白质在细胞中长什么样)
2. [不对称单元与生物组装体的区别](#2-不对称单元与生物组装体的区别)
3. [mmCIF 中的 Assembly 定义](#3-mmcif-中的-assembly-定义)
4. [对称操作：从单元扩展到组装体](#4-对称操作从单元扩展到组装体)
5. [扩展之后的三个核心挑战](#5-扩展之后的三个核心挑战)
6. [构建流程全览](#6-构建流程全览)
7. [解决身份问题：链重命名与等价标注](#7-解决身份问题链重命名与等价标注)
8. [解决质量问题：展开后的过滤](#8-解决质量问题展开后的过滤)
9. [解决视角问题：训练采样策略](#9-解决视角问题训练采样策略)
10. [最终输出：bioassembly_dict](#10-最终输出bioassembly_dict)
11. [各场景示例对比](#11-各场景示例对比)

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

## 5. 扩展之后的三个核心挑战

到第 4 节为止，我们已经从 CIF 文件中还原了蛋白质在生物体中的功能形态。但这只是一个包含三维坐标的原子集合，距离能被模型用于训练还有三个关键问题需要解决。

### 挑战一：身份混乱——扩展后的链"谁是谁"？

对称操作会复制链。以同二聚体为例，扩展后两条链都叫"A"。但模型不仅需要区分名字，还需要精确理解三层语义：

| 需要区分的信息 | 为什么模型需要知道 | 对应的下游模块 |
|-------------|---------------|------------|
| 这是两条**不同的物理链** | 注意力计算中必须区分不同链，否则两条链的特征会混淆 | 相对位置编码、注意力机制 |
| 它们是同一种蛋白质的**拷贝** | 相同序列的不同拷贝应有相似结构，可提供额外训练约束 | 对称性损失函数 |
| 它们在结构上**等价** | 交换两条链后是同一个结构，损失函数不应因交换而产生虚假惩罚 | 对称置换机制 |

如果分不清这三层身份，位置编码会错位、损失函数会产生虚假惩罚、对称性训练信号会丢失。

### 挑战二：数据质量——哪些数据不能用？

PDB 中的结构数据并非完美，扩展后可能出现：
- **链数爆炸**：病毒衣壳等大型对称复合物有 60 条以上链，远超 GPU 显存容量（Transformer 计算量随序列长度二次方增长）
- **坐标碰撞**：对称操作产生的某些链与其他链在空间上严重重叠，属于晶体包装产物而非真实结构

这些问题若不处理，模型会从物理上不可能的数据中学习，导致预测能力下降。

### 挑战三：训练视角——一个结构有多少种学习目标？

一个包含蛋白质 A + 蛋白质 B + 配体 C 的复合物，蕴含多种学习视角：

```
单链视角：  A 的内部折叠、B 的内部折叠      → 学习"蛋白质如何折叠"
接口视角：  A-B 结合界面、A-C 结合界面      → 学习"分子如何结合"
```

模型训练时需要明确"这个样本重点学什么"。如果把整个复合物混在一起训练，梯度信号会在所有视角间稀释；如果将每种视角拆开组织，训练更高效、评估也能精确到每种任务。

---

接下来的几节按"**流程总览 → 身份 → 质量 → 视角 → 打包**"的顺序，逐一介绍 Protenix 如何解决这些挑战。

---

## 6. 构建流程全览

`get_bioassembly()` 是 Protenix 中最复杂的方法之一（`parser.py` 约 745-913 行）。以下是完整的处理流程，每一步标注了它在解决哪个挑战：

```
输入：.cif 或 .cif.gz 文件
                │
    ┌───────────▼───────────┐
    │  Step 1: 链数校验      │  oligomeric_count > 1000 或 ? → 返回空
    │  【质量】早期拦截       │
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 2: 加载不对称单元  │  get_structure() → AtomArray
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 3: 预过滤（8项）  │  移除水/氢/未知链等噪声数据
    │  【质量】基础清洗       │
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
    │  【质量】统一缺失值     │
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 7: 链重命名      │  unique_chain_and_add_ids()
    │  【身份】消除命名冲突   │
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 8: 后过滤（2项）  │  too_many_chains_filter + remove_clashing_chains
    │  【质量】控制规模与碰撞 │
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 9: 等价标注      │  find_equiv_mol_and_assign_ids()
    │  【身份】识别等价副本   │
    └───────────┬───────────┘
                │
    ┌───────────▼───────────┐
    │  Step 10: 参考空间编码  │  add_ref_space_uid() + add_ref_info_and_res_perm()
    └───────────┬───────────┘
                │
            bioassembly_dict
                │
    ┌───────────▼───────────┐
    │  采样索引生成           │  make_chain_indices() + make_interface_indices()
    │  【视角】拆分训练目标   │
    └───────────┬───────────┘
                │
        sample_indices_list
```

### Step 3：预过滤（8 项清洗操作）

在展开组装体**之前**对不对称单元做清洗——越早清除噪声，展开时复制的无用数据越少：

| 序号 | 操作 | 为什么需要 |
|------|------|---------|
| 1 | 移除水分子 | 溶剂水是结晶条件的产物，不参与生物功能预测 |
| 2 | 移除氢原子 | X 射线分辨率通常不足以确定氢原子位置，保留只增加噪声 |
| 3 | 移除全未知链 | 所有残基均标注为 UNK 的链没有序列信息，无法提供训练信号 |
| 4 | 移除 Cα 过远链 | 相邻 Cα 距离超过 10 Å 说明链"断裂"，序列与坐标不对应 |
| 5 | 修正精氨酸 | 部分 PDB 文件中精氨酸原子标注顺序有误，影响坐标正确性 |
| 6 | 补全缺失原子 | 结构连接表中有键但原子缺失，补全以维持拓扑完整性 |
| 7 | MSE → MET | 硒代甲硫氨酸是晶体学实验标记，替换为标准甲硫氨酸 |
| 8 | 移除元素 X | 未知元素的虚拟原子无物理意义 |

---

## 7. 解决身份问题：链重命名与等价标注

本节对应流程中的 Step 7 和 Step 9，解决[挑战一](#挑战一身份混乱扩展后的链谁是谁)中提出的"扩展后的链谁是谁"问题。

### 7.1 链重命名：给每条链一个唯一名字

**问题**：展开组装体后，同一链的多个副本名称相同（都叫"A"），代码中无法区分。

**方案**（`unique_chain_and_add_ids()`，`parser.py` 约第 3049 行）：
- 第一次出现的链保持原名（如 `A`）
- 后续副本加计数后缀（如 `A.1`、`A.2`）

```
原始：  [A, B, A, B, C]
重命名：[A, B, A.1, B.1, C]
```

### 7.2 三个整数 ID：编码三层身份信息

仅仅重命名只解决了"名字不同"的表面问题，但模型需要的是**语义层面**的身份区分。因此，重命名的同时为每个原子写入三个整数 ID：

```
链：          [A,  B,  A.1, B.1, C]
asym_id_int:  [0,  1,   2,   3,  4]    ← 每条链的唯一序号
entity_id_int:[0,  1,   0,   1,  2]    ← 同种蛋白质共享相同值
sym_id_int:   [0,  0,   1,   1,  0]    ← 同种蛋白质的第几个副本
```

这三个 ID 分别服务于不同的下游需求：

| ID | 编码的信息 | 下游如何使用 |
|----|---------|-----------|
| `asym_id_int` | "这是第几条物理链" | **相对位置编码**：模型通过此 ID 判断两个 token 是否属于同一条链。同一链内的 token 有序列邻近关系，跨链的没有，编码方式必须不同 |
| `entity_id_int` | "这是哪种蛋白质" | **对称置换**：损失计算时，entity_id 相同的链可以互换——如果交换后损失更小，就采用交换后的对应关系，避免虚假惩罚（详见[损失与优化指南](loss_and_optimization_guide.md)第 13 节） |
| `sym_id_int` | "这是该蛋白质的第几份副本" | **位置编码 + 数据增强**：让模型区分同种蛋白质的不同副本，训练时可随机打乱副本编号以消除位置偏见 |

### 7.3 等价分子标注：识别"结构上完全相同的副本"

**问题**：`entity_id_int` 只表明两条链的序列相同，但不保证它们在当前组装体中的**原子组成**完全一致。例如，同一种蛋白质的两个副本可能因实验条件差异，一个有完整侧链坐标，另一个部分缺失。

**为什么还需要更精确的等价判断？** 在损失函数的对称置换中，Protenix 允许对等价链做交换——如果预测的链 A 结构更像真值的 A.1，就把 A 和 A.1 对调后再算损失。但交换的前提是两条链的原子组成**完全一致**，否则交换后原子对不上，损失计算会出错。

**方案**（`find_equiv_mol_and_assign_ids()`，`parser.py` 约第 2869 行）：

1. 通过共价连接性将所有原子划分为独立分子 → 每个分子获得唯一 `mol_id`
2. 按 `entity_id` 分组找到"候选等价分子"
3. 逐一比较原子名称序列——只有名称序列**完全一致**的分子才标记为等价，共享相同的 `entity_mol_id`
4. `mol_atom_index` 记录原子在分子内的位置，使等价分子间的原子可以一一对应

**示例**（同二聚体）：

```
链 A（蛋白X）    + 链 A.1（蛋白X 的副本）
原子名称序列都是：N, CA, C, O, CB, ...（完全一致）

mol_id:          [0, 0, 0, ..., 1, 1, 1, ...]     ← 两个独立分子
entity_mol_id:   [5, 5, 5, ..., 5, 5, 5, ...]     ← 等价！可安全交换
mol_atom_index:  [0, 1, 2, ..., 0, 1, 2, ...]     ← 一一对应，交换后原子不会错位
```

---

## 8. 解决质量问题：展开后的过滤

本节对应流程中的 Step 8，解决[挑战二](#挑战二数据质量哪些数据不能用)中的数据质量问题。

### 8.1 链数量过滤——控制计算规模

**问题**：大型对称复合物（如病毒衣壳，60 聚体）展开后链数远超 GPU 处理能力。

**方案**（`too_many_chains_filter()`）：
- 最多保留 **20 条**聚合物链 / **5120 个** token
- 优先保留 `assembly_1 == True` 的链——这些是 CIF 文件标注为首选组装体的链，包含最核心的结构信息

### 8.2 碰撞链过滤——移除不合理坐标

**问题**：对称操作有时会产生坐标与其他链严重重叠的"幽灵链"——两条链的原子几乎占据同一空间位置，这在物理上不可能发生。

**方案**（`remove_clashing_chains()`）：
- 若一条链有超过 **30%** 的原子与其他链原子的距离 < **1.7 Å**（非氢原子的最小物理距离），判定为碰撞链并移除
- 碰撞双方中，优先保留 Assembly 1 的链和原子数更多的链

### 8.3 为什么预过滤和后过滤分成两个阶段？

**预过滤**（Step 3）在展开**之前**执行，清除水、氢、未知链等显而易见的噪声——越早清除，展开时复制的无用数据越少。

**后过滤**（Step 8）在展开**之后**执行，因为链碰撞和链数量问题只有在展开完成后才会出现。

---

## 9. 解决视角问题：训练采样策略

本节解决[挑战三](#挑战三训练视角一个结构有多少种学习目标)——如何从一个复合物中提取多种有意义的训练样本。

### 9.1 为什么不直接用整个结构训练？

一个包含蛋白质-蛋白质-配体的复合物同时包含链内折叠、链间结合、配体对接等多种信息。不区分视角会导致三个问题：
- **评估无法精细化**：无法单独衡量"蛋白质折叠预测是否准确"和"配体对接预测是否准确"
- **采样无法均衡**：PDB 中蛋白质-蛋白质复合物远多于蛋白质-配体复合物，不做分类就无法通过加权采样保证数据多样性
- **裁剪无法聚焦**：训练时会将大结构裁剪为 384 个 token 的子集（详见[预处理完全指南](preprocessing_guide.md)），采样索引告诉裁剪模块"应该围绕哪条链/哪个界面裁剪"

### 9.2 两类采样视角

**链采样（Chain Sample）**：每条聚合物链产生一个样本，学习目标是**链内部折叠**。

**界面采样（Interface Sample）**：每对空间接近（代表原子距离 < 5 Å）的链产生一个样本，学习目标是**链间结合方式**。界面检测使用空间哈希网格（CellList）加速邻域搜索。

一个复合物会产生多个采样记录。例如：

```
蛋白质 A + 蛋白质 B + 配体 C（A-B 接触，A-C 接触，B-C 不接触）

产生的采样记录：
  ① {"type": "chain",     "chain_1": "A"}                     → 学习 A 的折叠
  ② {"type": "chain",     "chain_1": "B"}                     → 学习 B 的折叠
  ③ {"type": "interface", "chain_1": "A", "chain_2": "B"}     → 学习 A-B 蛋白结合
  ④ {"type": "interface", "chain_1": "A", "chain_2": "C"}     → 学习 A-配体结合
```

### 9.3 按分子类型分组——均衡训练数据

每个采样记录标注 `mol_type_group`，用于控制训练批次的组成多样性：

| mol_type_group | 含义 | 训练价值 |
|----------------|------|---------|
| `intra_prot` | 纯蛋白质 | 学习蛋白质折叠 |
| `prot_nuc` | 蛋白质-核酸 | 学习转录/翻译相关结合 |
| `prot_ligand` | 蛋白质-配体 | 学习药物结合（药物设计核心） |
| `nuc_nuc` | 核酸-核酸 | 学习 RNA/DNA 结构 |
| `nuc_ligand` | 核酸-配体 | 学习核酸-小分子相互作用 |
| `ligand_ligand` | 配体-配体 | 较罕见 |

通过对 `mol_type_group` 加权采样，训练时可以确保模型在每种任务上都获得足够的学习机会，即使 PDB 中各类型的结构数量极不均衡。

---

## 10. 最终输出：bioassembly_dict

经过上述所有步骤，`get_bioassembly()` 返回一个字典。下表为每个字段标注了它的来源和用途：

```python
bioassembly_dict = {
    # ── 元数据（用于日期过滤和数据质量控制）─────────────────
    "pdb_id": str,                       # PDB ID，如 "102m"
    "assembly_id": str,                  # 组装体编号，默认 "1"
    "release_date": str,                 # 发布日期 → 模板搜索的时间隔离
    "resolution": float,                 # 分辨率 → 过滤低质量数据、控制置信度损失

    # ── 序列信息（用于 MSA 搜索和序列编码）─────────────────
    "sequences": dict[str, str],         # entity_id → 氨基酸序列 → MSA 搜索输入
    "entity_poly_type": dict[str, str],  # entity_id → 聚合物类型 → 区分蛋白/DNA/RNA

    # ── 结构数据（核心，承载所有坐标和身份标注）─────────────
    "atom_array": AtomArray | None,      # 展开、过滤、标注后的完整原子数组
    "token_array": TokenArray,           # Token 化后的表示 → 模型输入的基础

    # ── 统计信息（用于数据集级别的过滤和采样决策）──────────
    "num_assembly_polymer_chains": int,  # 原始聚合物链总数 → 数据集统计
    "num_prot_chains": int,              # 蛋白质链数量 → 采样分组依据
    "num_tokens": int,                   # Token 总数 → 裁剪决策依据
}
```

其中 `atom_array` 承载了前述所有身份标注（`asym_id_int`、`entity_id_int`、`sym_id_int`、`mol_id`、`entity_mol_id` 等），这些标注在后续的预处理、损失计算和评估中被反复读取。

连同采样索引一起，最终输出写入 `.pkl.gz` 文件：`(sample_indices_list, bioassembly_dict)`。

---

## 11. 各场景示例对比

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
