# Cellular-basis-of-sex-role-reversal
Code of article：“Cellular basis of sex-role reversal: male pregnancy in seahorses” by Yali Liu, Han Jiang, Yuanxiang Miao, Luonan Chen* , Axel Meyer* , Qiang Lin* et al.

![pic](https://github.com/miaoyx-SHTU/Cellular-basis-of-sex-role-reversal/blob/main/Image/fig1.png)

# 🐴 Seahorse Placental Development Single-Cell Analysis

<div align="center">

![Seahorse](https://img.shields.io/badge/Species-Seahorse-blue)
![scRNA-seq](https://img.shields.io/badge/Technology-scRNA--seq-green)
![ATAC-seq](https://img.shields.io/badge/Technology-ATAC--seq-orange)
![R](https://img.shields.io/badge/Language-R-blue)
![Python](https://img.shields.io/badge/Language-Python-yellow)

*A comprehensive single-cell multi-omics study revealing the molecular mechanisms underlying seahorse placental development*

</div>

## 📖 研究概述 / Overview

本研究通过单细胞RNA测序(scRNA-seq)和单细胞ATAC测序(scATAC-seq)技术，系统性地解析了海马(Seahorse)胎盘发育过程中的细胞类型多样性和基因调控网络。研究涵盖了多个发育时期的样本，并与人类和小鼠的胎盘发育进行了跨物种比较分析。

This study systematically characterizes cell type diversity and gene regulatory networks during seahorse placental development using single-cell RNA sequencing (scRNA-seq) and single-cell ATAC sequencing (scATAC-seq). The research includes multiple developmental timepoints and performs cross-species comparative analysis with human and mouse placental development.

## 🗂️ 项目结构 / Project Structure

```
sh_code/
├── 📂 Notebooks/          # Jupyter notebooks for analysis workflows
├── 📂 R code/            # R scripts for specialized analyses  
├── 📂 Image/             # Figures and visualization outputs
└── 📄 manuscript-NATECOLEVOL-25020487A.pdf  # Associated research paper
```

## 📓 Notebooks 文件夹详细说明

### 🔬 主要分析流程

#### `SCpipline.ipynb` - 单细胞RNA测序数据处理主流程
- **功能**: 完整的scRNA-seq数据预处理、质量控制和细胞类型注释流程
- **主要步骤**:
  - 10X Genomics数据读取和Seurat对象创建
  - 多批次数据整合 (BP1/BP2, ZQ1/ZQ2, HY1/HY2等发育时期)
  - 基因注释转换（海马基因ID到基因名）
  - 线粒体和核糖体基因比例计算
  - Harmony批次效应校正
  - UMAP/t-SNE降维可视化
  - 聚类分析和标记基因识别
- **输出**: 处理后的Seurat对象、质量控制图表、细胞类型注释结果

#### `Run_SAMAP.ipynb` - 跨物种细胞类型比较分析
- **功能**: 使用SAMap算法进行海马与人类、小鼠的跨物种细胞类型映射
- **核心分析**:
  - 海马-小鼠胎盘细胞类型对应关系分析
  - 海马-人类胎盘细胞类型对应关系分析
  - 生成Sankey图展示细胞类型映射结果
  - 识别同源基因对和细胞类型特异性标记基因
- **特色**: 创新性地将海马胎盘细胞与哺乳动物胎盘进行比较，揭示胎盘功能的进化保守性

#### `Cell_commu.ipynb` - 细胞间通讯网络分析
- **功能**: 基于CellPhoneDB算法分析细胞间配体-受体相互作用
- **分析内容**:
  - 构建细胞类型间通讯网络
  - 识别关键的配体-受体对
  - 生成弦图(Chord plot)和Sankey图可视化
  - 分析不同发育时期的通讯模式变化
- **可视化**: 包含网络图、热图、小提琴图等多种展示方式

#### `Cross_species_Evolutionary_tree.ipynb` - 跨物种进化系统发育分析
- **功能**: 构建海马、薛氏海马、洪福海马的细胞类型进化树
- **技术特点**:
  - 基于同源基因的表达模式构建系统发育树
  - 集成聚类(Ensemble clustering)提高树构建的稳定性
  - Bootstrap支持度评估
  - 识别系统发育树各节点的特征基因
- **创新性**: 首次在海马类群中进行细胞类型的系统发育分析

### 🎯 专业化分析

#### `DNB&chordplot.ipynb` - 动态网络生物标志物分析
- **功能**: 识别发育过程中的关键调控网络和生物标志物
- **方法**: 动态网络分析结合弦图可视化

#### `Enrich_for_SAMap.ipynb` - SAMap结果功能富集分析
- **功能**: 对跨物种比较结果进行GO/KEGG功能富集分析
- **输出**: 功能注释和通路富集结果

#### `Raw_BPHYZQ.ipynb` - 原始数据探索性分析
- **功能**: 初步的数据质量评估和探索性数据分析
- **内容**: 基础统计、数据分布检查、初步可视化

## 🧬 R Code 文件夹详细说明

### 🔬 ATAC-seq分析流程

#### `ATAC_code0.R` - ATAC-seq数据预处理和注释
- **主要功能**:
  - 基于ArchR框架的scATAC-seq数据处理
  - 海马基因组注释文件构建
  - Peak calling和质量控制
  - 与scRNA-seq数据整合进行细胞类型注释
- **技术细节**:
  - 自定义海马参考基因组(Seahorse.bs.genome)
  - Doublet检测和过滤
  - 迭代LSI降维和Harmony批次校正
  - 基因得分矩阵计算

#### `ATAC_code1.R` - 进阶ATAC-seq分析
- **核心分析**:
  - Peak calling和motif富集分析
  - 转录因子结合位点预测
  - 拟时间轨迹分析
  - ARE (Androgen Response Element) motif特异性分析
- **特色功能**:
  - 鉴定发育轨迹中的关键调控因子
  - 分析雄激素响应元件的作用模式

#### `seahorse.R` - 整合分析脚本
- **功能**: 整合scRNA-seq和scATAC-seq数据的综合分析
- **包含**:
  - 多组学数据整合
  - 联合聚类分析
  - 调控网络推断

## 📊 主要研究发现

### 🔍 细胞类型鉴定
- 鉴定出**14个主要细胞类型**，包括：
  - 滋养层细胞亚型 (TGC subtypes)
  - 基底细胞 (Basal cells) 
  - 蜕膜细胞 (Decidual cells)
  - 内皮细胞 (Endothelial cells)
  - 成纤维细胞 (Fibroblasts)
  - 免疫细胞 (Immune cells)

### 🧬 跨物种比较
- 发现海马胎盘滋养层细胞与哺乳动物胎盘的**功能同源性**
- 识别出保守的胎盘发育调控基因和通路
- 揭示了不同物种胎盘结构的进化关系

### 🔗 细胞通讯网络
- 构建了发育过程中的**细胞间通讯图谱**
- 识别关键的信号通路和调控因子
- 分析了不同发育阶段的通讯模式变化

## 🛠️ 技术方法

### 📈 生物信息学工具
- **R packages**: Seurat, ArchR, Harmony, CellPhoneDB, ggtree
- **Python packages**: SAMap, scanpy, pandas, matplotlib
- **统计方法**: Wilcoxon检验, 相关性分析, 聚类分析

### 🔬 实验技术
- **scRNA-seq**: 10X Genomics Chromium平台
- **scATAC-seq**: 10X Genomics单细胞ATAC测序
- **样本**: 多个海马发育时期 (BP, HY, ZQ等)

## 📋 使用说明

### 环境配置
```r
# R dependencies
install.packages(c("Seurat", "harmony", "dplyr", "ggplot2"))
BiocManager::install(c("ArchR", "GenomicRanges"))

# Python dependencies  
pip install scanpy samap pandas matplotlib seaborn
```

### 运行流程
1. **数据预处理**: 运行 `SCpipline.ipynb`
2. **跨物种分析**: 运行 `Run_SAMAP.ipynb` 
3. **细胞通讯**: 运行 `Cell_commu.ipynb`
4. **进化分析**: 运行 `Cross_species_Evolutionary_tree.ipynb`
5. **ATAC分析**: 运行 `ATAC_code0.R` 和 `ATAC_code1.R`

## 📚 引用信息

如果您使用本代码进行研究，请引用以下论文：

```bibtex
@article{seahorse_placenta_2024,
  title={Single-cell multi-omics reveals molecular mechanisms of seahorse placental development},
  author={[Authors]},
  journal={Nature Ecology & Evolution},
  year={2024},
  doi={[DOI]}
}
```

## 👥 贡献者

- 主要分析: [研究团队]
- 代码开发: [开发团队]
- 数据收集: [实验团队]

## 📞 联系方式

如有问题或建议，请通过以下方式联系：
- Email: [contact@email.com]
- Issues: 请在GitHub上提交issue

## 📄 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

---

<div align="center">

**🌟 如果这个项目对您有帮助，请给我们一个Star！**

![GitHub stars](https://img.shields.io/github/stars/username/repo?style=social)

</div>
