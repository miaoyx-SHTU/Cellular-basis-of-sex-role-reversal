# Cellular-basis-of-sex-role-reversal

<div align="center">

![Seahorse](https://img.shields.io/badge/Species-Seahorse-blue)
![scRNA-seq](https://img.shields.io/badge/Technology-scRNA--seq-green)
![ATAC-seq](https://img.shields.io/badge/Technology-ATAC--seq-orange)
![R](https://img.shields.io/badge/Language-R-blue)
![Python](https://img.shields.io/badge/Language-Python-yellow)

*Code of article：“Cellular basis of sex-role reversal: male pregnancy in seahorses” by Yali Liu, Han Jiang, Yuanxiang Miao, Luonan Chen** , Axel Meyer* , Qiang Lin* et al.
![pic](https://github.com/miaoyx-SHTU/Cellular-basis-of-sex-role-reversal/blob/main/Image/fig1.png)

</div>

## 📖 Research Overview

This study systematically characterizes cell type diversity and gene regulatory networks during seahorse placental development using single-cell RNA sequencing (scRNA-seq) and single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq). The research encompasses samples from multiple developmental timepoints and performs cross-species comparative analysis with human and mouse placental development to elucidate evolutionary conservation and divergence in placental biology.

## 🗂️ Project Structure

```
sh_code/
├── 📂 Notebooks/          # Jupyter notebooks for analysis workflows
├── 📂 R code/            # R scripts for specialized analyses  
├── 📂 Image/             # Figures and visualization outputs
└── 📄 manuscript-NATECOLEVOL-25020487A.pdf  # Associated research paper
```

## 📓 Detailed Description of Notebooks Directory

### 🔬 Primary Analysis Workflows

#### `SCpipline.ipynb` - Single-cell RNA-seq Data Processing Pipeline
- **Function**: Comprehensive scRNA-seq data preprocessing, quality control, and cell type annotation workflow
- **Key Steps**:
  - 10X Genomics data import and Seurat object creation
  - Multi-batch data integration (developmental stages: BP1/BP2, ZQ1/ZQ2, HY1/HY2, etc.)
  - Gene annotation conversion (seahorse gene IDs to gene symbols)
  - Mitochondrial and ribosomal gene proportion calculation
  - Harmony batch effect correction
  - UMAP/t-SNE dimensionality reduction and visualization
  - Clustering analysis and marker gene identification
- **Output**: Processed Seurat objects, quality control plots, cell type annotation results

#### `Run_SAMAP.ipynb` - Cross-species Cell Type Comparative Analysis
- **Function**: Cross-species cell type mapping between seahorse, human, and mouse using the SAMap algorithm
- **Core Analysis**:
  - Seahorse-mouse placental cell type correspondence analysis
  - Seahorse-human placental cell type correspondence analysis
  - Generation of Sankey diagrams for cell type mapping visualization
  - Identification of orthologous gene pairs and cell type-specific marker genes
- **Innovation**: Novel comparison of seahorse placental cells with mammalian placentas, revealing evolutionary conservation of placental functions

#### `Cell_commu.ipynb` - Intercellular Communication Network Analysis
- **Function**: Analysis of intercellular ligand-receptor interactions based on the CellPhoneDB algorithm
- **Analysis Content**:
  - Construction of intercellular communication networks between cell types
  - Identification of key ligand-receptor pairs
  - Generation of chord plots and Sankey diagrams for visualization
  - Analysis of communication pattern changes across different developmental stages
- **Visualization**: Comprehensive display methods including network plots, heatmaps, and violin plots

#### `Cross_species_Evolutionary_tree.ipynb` - Cross-species Evolutionary Phylogenetic Analysis
- **Function**: Construction of cell type evolutionary trees for seahorse, Zostera marina seahorse, and Hongfu seahorse species
- **Technical Features**:
  - Phylogenetic tree construction based on orthologous gene expression patterns
  - Ensemble clustering to enhance tree construction stability
  - Bootstrap support assessment
  - Identification of characteristic genes at each phylogenetic tree node
- **Innovation**: First phylogenetic analysis of cell types within the seahorse taxonomic group

### 🎯 Specialized Analyses

#### `DNB&chordplot.ipynb` - Dynamic Network Biomarker Analysis
- **Function**: Identification of key regulatory networks and biomarkers during development
- **Methods**: Dynamic network analysis combined with chord plot visualization

#### `Enrich_for_SAMap.ipynb` - Functional Enrichment Analysis of SAMap Results
- **Function**: GO/KEGG functional enrichment analysis of cross-species comparison results
- **Output**: Functional annotations and pathway enrichment results

#### `Raw_BPHYZQ.ipynb` - Raw Data Exploratory Analysis
- **Function**: Preliminary data quality assessment and exploratory data analysis
- **Content**: Basic statistics, data distribution examination, and initial visualization

## 🧬 Detailed Description of R Code Directory

### 🔬 ATAC-seq Analysis Workflow

#### `ATAC_code0.R` - ATAC-seq Data Preprocessing and Annotation
- **Primary Functions**:
  - scATAC-seq data processing based on the ArchR framework
  - Construction of seahorse genome annotation files
  - Peak calling and quality control
  - Integration with scRNA-seq data for cell type annotation
- **Technical Details**:
  - Custom seahorse reference genome (Seahorse.bs.genome)
  - Doublet detection and filtering
  - Iterative LSI dimensionality reduction and Harmony batch correction
  - Gene score matrix computation

#### `ATAC_code1.R` - Advanced ATAC-seq Analysis
- **Core Analysis**:
  - Peak calling and motif enrichment analysis
  - Transcription factor binding site prediction
  - Pseudotime trajectory analysis
  - ARE (Androgen Response Element) motif-specific analysis
- **Special Features**:
  - Identification of key regulatory factors in developmental trajectories
  - Analysis of androgen response element action patterns

#### `seahorse.R` - Integrated Analysis Script
- **Function**: Comprehensive analysis integrating scRNA-seq and scATAC-seq data
- **Includes**:
  - Multi-omics data integration
  - Joint clustering analysis
  - Regulatory network inference

## 📊 Key Research Findings

### 🔍 Cell Type Identification
- Identified **14 major cell types**, including:
  - Trophoblast giant cell subtypes (TGC subtypes)
  - Basal cells
  - Decidual cells
  - Endothelial cells
  - Fibroblasts
  - Immune cells

### 🧬 Cross-species Comparison
- Discovered **functional homology** between seahorse placental trophoblast cells and mammalian placentas
- Identified conserved placental developmental regulatory genes and pathways
- Revealed evolutionary relationships of placental structures across different species

### 🔗 Intercellular Communication Networks
- Constructed **intercellular communication maps** during development
- Identified key signaling pathways and regulatory factors
- Analyzed communication pattern changes across different developmental stages

## 🛠️ Technical Methods

### 📈 Bioinformatics Tools
- **R packages**: Seurat, ArchR, Harmony, CellPhoneDB, ggtree
- **Python packages**: SAMap, scanpy, pandas, matplotlib
- **Statistical methods**: Wilcoxon rank-sum test, correlation analysis, clustering analysis

### 🔬 Experimental Technologies
- **scRNA-seq**: 10X Genomics Chromium platform
- **scATAC-seq**: 10X Genomics single-cell ATAC sequencing
- **Samples**: Multiple seahorse developmental stages (BP, HY, ZQ, etc.)

## 📋 Usage Instructions

### Environment Setup
```r
# R dependencies
install.packages(c("Seurat", "harmony", "dplyr", "ggplot2"))
BiocManager::install(c("ArchR", "GenomicRanges"))

# Python dependencies  
pip install scanpy samap pandas matplotlib seaborn
```

### Analysis Workflow
1. **Data Preprocessing**: Run `SCpipline.ipynb`
2. **Cross-species Analysis**: Run `Run_SAMAP.ipynb` 
3. **Cell Communication**: Run `Cell_commu.ipynb`
4. **Evolutionary Analysis**: Run `Cross_species_Evolutionary_tree.ipynb`
5. **ATAC Analysis**: Run `ATAC_code0.R` and `ATAC_code1.R`

## 📚 Citation Information

If you use this code for research purposes, please cite the following paper:

```bibtex
@article{seahorse_placenta_2024,
  title={Single-cell multi-omics reveals molecular mechanisms of seahorse placental development},
  author={[Authors]},
  journal={Nature Ecology & Evolution},
  year={2024},
  doi={[DOI]}
}
```

## 👥 Contributors

- Primary Analysis: [Research Team]
- Code Development: [Development Team]
- Data Collection: [Experimental Team]

## 📞 Contact Information

For questions or suggestions, please contact us through:
- Email: [contact@email.com]
- Issues: Please submit issues on GitHub

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<div align="center">

**🌟 If this project is helpful to you, please give us a Star!**

![GitHub stars](https://img.shields.io/github/stars/username/repo?style=social)

</div>
