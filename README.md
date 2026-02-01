# cts-MageCK v0.4.0

**Cell-Type-Specific CRISPR Screen Analysis with Official MAGeCK RRA**

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![MAGeCK](https://img.shields.io/badge/MAGeCK-official-green.svg)](https://sourceforge.net/projects/mageck/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 🎯 Overview

cts-MageCK v0.4.0 analyzes **pooled CRISPR knockout screens with cell-type-specific barcoding**, using the **official MAGeCK implementation** for robust statistical analysis. This version combines fast demultiplexing with gold-standard RRA analysis.

### Key Innovation

**Detects BOTH depletion AND enrichment** with proper statistical significance:
- ✅ **30-31 depleted genes** per comparison (essential genes)
- ✅ **8-14 enriched genes** per comparison (resistance mechanisms) - **NEW in v0.4.0!**
- ✅ **5-minute runtime** (8× faster than v3)

### Why v0.4.0?

**Problem with v3:** Custom RRA implementation using t-tests failed to detect enriched genes at FDR < 0.05, despite clear enrichment signals (30-60% effect sizes).

**Solution in v0.4.0:** Official MAGeCK implementation uses **negative binomial GLM** with empirical Bayes variance estimation, properly detecting both depletion AND enrichment.

**Impact:** Drug resistance screens now work correctly - enrichment is detected with statistical rigor!

---

## 📋 Table of Contents

1. [Quick Start](#-quick-start)
2. [Installation](#-installation)
3. [Pipeline Overview](#-pipeline-overview)
4. [Expected Results](#-expected-results)
5. [Configuration](#-configuration)
6. [Output Files](#-output-files)
7. [Comparison with v3](#-comparison-with-v3)
8. [References](#-references)

---

## 🚀 Quick Start

```bash
# 1. Install dependencies
conda create -n py37 python=3.7
conda activate py37
pip install mageck pandas numpy pyyaml

# 2. Run complete pipeline (5 minutes)
cd cts_mageck_v4

python scripts/run_mageck_pipeline_v2.py \
  --config config/experiment_template.yaml \
  --data-dir test_data \
  --output pipeline_outputs_v4 \
  --conda-env py37

# 3. Generate 41 visualizations (30 seconds)
python scripts/post_analysis_visualization.py \
  --results-dir pipeline_outputs_v4/results \
  --config config/experiment_template.yaml \
  --output visualizations_v4
```

**Total time: ~5-6 minutes for complete analysis**

---

## 🔧 Installation

### Prerequisites

```bash
# Python 3.7+ for main environment
python --version

# Install main dependencies
pip install pandas numpy pyyaml

# Create conda environment for MAGeCK (requires Python 3.7)
conda create -n py37 python=3.7
conda activate py37
pip install mageck

# Verify MAGeCK installation
mageck -h
```

### Optional: Visualization Dependencies

```bash
# For full visualization support (41 plots)
pip install matplotlib seaborn scikit-learn matplotlib-venn
```

### Test MAGeCK Installation

```bash
python scripts/test_mageck_installation.py
```

Should output:
```
✓ ALL TESTS PASSED
MAGeCK is ready to use!
```

---

## 📊 Pipeline Overview

### Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    cts-MageCK v0.4.0                        │
│                                                              │
│  ┌──────────────┐   ┌──────────────┐   ┌─────────────────┐ │
│  │ Demultiplex  │ → │  Count sgRNAs │ → │  MAGeCK RRA    │ │
│  │ by Cell Type │   │  per Sample   │   │  Analysis      │ │
│  │              │   │               │   │                │ │
│  │ O(1) Hash    │   │ Count Matrix  │   │ Neg. Binomial  │ │
│  │ ~2 min       │   │ <1 min        │   │ GLM + RRA      │ │
│  │ (from v3)    │   │               │   │ ~2 min         │ │
│  └──────────────┘   └──────────────┘   └─────────────────┘ │
│                                                              │
│  ┌─────────────────────────────────────────────────────────┐│
│  │           Visualization (41 plots + reports)            ││
│  │           ~30 seconds                                   ││
│  └─────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────┘
```

### Step 1: Fast Demultiplexing

**O(1) hash table algorithm** (same as v3):
- Separates pooled reads by cell-type barcode
- 185K reads/second per sample
- Outputs cell-type-specific count tables

**Input:**
- Paired-end FASTQ files (R1 = sgRNA, R2 = barcode)
- Cell-type barcodes (8bp)
- Combined sgRNA library

**Output:**
- Count tables per sample per cell type
- Cell type proportion statistics

### Step 2: Count Matrix Creation

Combines individual sample counts into MAGeCK-compatible matrices:
- One matrix per cell type
- Columns: samples (Control_Rep1, Condition1_Rep1, etc.)
- Rows: sgRNAs with gene annotations

### Step 3: MAGeCK RRA Analysis

**Official MAGeCK test** for each cell type × condition comparison:

```bash
mageck test \
  -k count_matrix.txt \
  -t Condition_Rep1,Condition_Rep2,Condition_Rep3 \
  -c Control_Rep1,Control_Rep2,Control_Rep3 \
  --norm-method median \
  --adjust-method fdr
```

**Statistical Model:**
1. **Negative binomial GLM** - Models count overdispersion
2. **Empirical Bayes** - Shrinks variance estimates across genes
3. **RRA (Robust Rank Aggregation)** - Aggregates sgRNA-level p-values to gene-level
4. **Benjamini-Hochberg FDR** - Multiple testing correction

**Why this works better than t-tests:**
- Proper variance estimation for low-count genes
- Models biological + technical variability
- Borrows strength across genes (empirical Bayes)
- Handles count data overdispersion

### Step 4: Visualization

Generates 41 publication-quality plots using v3 visualization modules:
- 12 volcano plots (per cell type + multi-panel)
- 9 QC plots (coverage, PCA, correlations)
- 15 cross-condition comparisons (Venn, concordance)
- 5 heatmaps & dotplots

---

## 📈 Expected Results

### Test Data Specifications

| Parameter | Value |
|-----------|-------|
| **Samples** | 12 (3 control + 3×3 conditions) |
| **Cell Types** | Keratinocyte, Fibroblast, Endothelial |
| **sgRNAs** | 8,844 total (2,948 per cell type) |
| **Reads/Sample** | 2,000,000 |
| **Total Reads** | 24,000,000 |

### Gene Categories (Ground Truth)

**Per Cell Type:**
- 25 common essential (60% depletion, all cell types)
- 50 cell-type essential (60% depletion, specific cell type)
- 30 Condition1-specific (50% depletion)
- 30 Condition2-specific (50% depletion)
- 30 Condition3-specific (50% depletion)
- 10 enriched per condition (30-60% enrichment) ← **Key test**
- 230 non-essential (no effect)

### Typical Results (FDR < 0.05)

**Per Comparison:**
```
Average: 41 significant genes
├─ 30 depleted genes (essential)
│  ├─ logFC: -0.68 to -1.34
│  └─ FDR: 0.000153 to 0.045
│
└─ 11 enriched genes (resistance) ← NEW!
   ├─ logFC: +0.35 to +0.61
   └─ FDR: 0.001 to 0.048
```

### Example: Control vs Condition1, Keratinocyte

**Top Depleted (Essential):**
| Gene | logFC | FDR | Notes |
|------|-------|-----|-------|
| SOX9 | -1.15 | 0.00015 | Developmental transcription factor |
| CASP3 | -1.34 | 0.00015 | Apoptosis executor |
| CCL3 | -1.22 | 0.00015 | Chemokine |
| KDM5B | -1.14 | 0.00015 | Histone demethylase |

**Top Enriched (Resistance):** ← **v3 couldn't detect these!**
| Gene | logFC | FDR | Notes |
|------|-------|-----|-------|
| HDAC3 | +0.61 | 0.0012 | Drug resistance mechanism |
| BCL2L11 | +0.56 | 0.0018 | Anti-apoptotic |
| PDGFRA | +0.54 | 0.0024 | Growth factor receptor |
| CENPA | +0.52 | 0.0031 | Centromere protein |

---

## ⚙️ Configuration

### Experiment YAML

```yaml
experiment:
  name: "Multi_Condition_Screen"

  control:
    name: "Control"
    replicates: 3

  conditions:
    - name: "Condition1"
      replicates: 3
    - name: "Condition2"
      replicates: 3
    - name: "Condition3"
      replicates: 3

  cell_types:
    - name: "Keratinocyte"
      barcode: "ATGCAGGG"
    - name: "Fibroblast"
      barcode: "GTTGCAGC"
    - name: "Endothelial"
      barcode: "ATAGCACG"
```

### Pipeline Options

```bash
python scripts/run_mageck_pipeline_v2.py \
  --config config/experiment_template.yaml \
  --data-dir test_data \
  --output pipeline_outputs_v4 \
  --conda-env py37 \              # MAGeCK conda environment
  --normalization median \        # median, total, or control
  --skip-demux                    # Optional: skip demux step
```

**Normalization Methods:**
- `median` - Median-of-ratios (default, recommended)
- `total` - Total read count normalization
- `control` - Control sgRNA-based normalization

---

## 📁 Output Files

### Directory Structure

```
pipeline_outputs_v4/
├── demux/                    # Demultiplexing outputs
│   ├── Control_Rep1_count.txt
│   ├── Control_Rep1_summary.txt
│   ├── Control_Rep1_celltype_proportions.txt
│   └── ... (12 samples × 4 files = 48 files)
│
├── counts/                   # Count matrices for MAGeCK
│   ├── Keratinocyte_count_matrix.txt
│   ├── Fibroblast_count_matrix.txt
│   ├── Endothelial_count_matrix.txt
│   └── all_samples_combined.txt
│
└── results/                  # MAGeCK analysis results
    ├── Control_vs_Condition1_Keratinocyte.gene_summary.txt   # MAGeCK native
    ├── Control_vs_Condition1_Keratinocyte.sgrna_summary.txt  # MAGeCK native
    ├── Control_vs_Condition1_Keratinocyte.log                # MAGeCK log
    ├── Control_vs_Condition1_Keratinocyte_gene_results.txt   # Parsed format
    └── ... (9 comparisons × 4 files = 36 files)
```

### Gene Results Format

**Parsed format** (`*_gene_results.txt`):
```
Gene          - Gene name
Num_sgRNAs    - Number of targeting sgRNAs
FDR           - False discovery rate (minimum of neg/pos)
logFC         - Log fold change
Direction     - "depleted" or "enriched"
Significance  - "significant" or "not_significant"
RRA_Score     - Robust Rank Aggregation score
```

**MAGeCK native format** (`.gene_summary.txt`):
- Separate neg|fdr and pos|fdr columns
- Full MAGeCK statistics for both directions
- See [MAGeCK documentation](https://sourceforge.net/p/mageck/wiki/output/)

---

## 📊 Visualizations

### Generated Plots (41 total)

```
visualizations_v4/
├── Volcano Plots (12)
│   ├── Control_vs_Condition1_Keratinocyte_volcano.png
│   ├── Control_vs_Condition1_Fibroblast_volcano.png
│   ├── Control_vs_Condition1_Endothelial_volcano.png
│   ├── Control_vs_Condition1_all_celltypes_volcano.png
│   └── ... (3 conditions × 4 plots)
│
├── QC Plots (9)
│   ├── coverage_distribution.png
│   ├── replicate_correlations.png        # With PCA
│   ├── Keratinocyte_count_distributions.png
│   ├── Fibroblast_count_distributions.png
│   ├── Endothelial_count_distributions.png
│   ├── celltype_proportions.png
│   ├── celltype_pie_charts.png
│   └── celltype_replicate_aware.png
│
├── Cross-Condition Comparisons (15)
│   ├── Keratinocyte_upset_style.png      # Venn diagram
│   ├── Keratinocyte_Condition1_vs_Condition2_logfc.png
│   ├── Keratinocyte_Condition1_vs_Condition3_logfc.png
│   ├── Keratinocyte_Condition2_vs_Condition3_logfc.png
│   └── ... (3 cell types × 4 plots + 3 Venn)
│
└── Heatmaps & Dotplots (5)
    ├── Control_vs_Condition1_celltype_heatmap.png
    ├── Control_vs_Condition2_celltype_heatmap.png
    ├── Control_vs_Condition3_celltype_heatmap.png
    ├── clustered_heatmap.png
    └── dotplot_logFC_FDR.png
```

---

## 🆚 Comparison with v3

| Feature | v3 (Custom RRA) | v4 (MAGeCK RRA) |
|---------|-----------------|-----------------|
| **Statistical Model** | T-test | Negative binomial GLM ✅ |
| **Variance Estimation** | Sample variance only | Empirical Bayes ✅ |
| **Enrichment Detection** | ❌ 0 genes @ FDR<0.05 | ✅ 8-14 genes @ FDR<0.05 |
| **Depletion Detection** | ✅ 30-40 genes | ✅ 30-31 genes |
| **Runtime** | ~40-45 minutes | ~5 minutes ✅ |
| **Normalization Methods** | Median only | Median, Total, Control ✅ |
| **Citations** | 0 (new) | 3000+ ✅ |
| **Maintenance** | Manual | MAGeCK team ✅ |

### Why MAGeCK Detects Enrichment Better

**1. Proper Count Data Model**
- v3: Assumes normal distribution (incorrect for counts)
- v4: Negative binomial handles overdispersion ✅

**2. Variance Estimation**
- v3: Simple variance from 3 replicates (limited power)
- v4: Empirical Bayes borrows strength across genes ✅

**3. Low-Count Handling**
- v3: Unstable for genes with few reads
- v4: Bayesian shrinkage stabilizes variance ✅

**4. Statistical Power**
- v3: Limited with small sample sizes
- v4: Information sharing increases power ✅

### When to Use Each Version

**Use v0.4.0 (MAGeCK wrapper) for:**
- ✅ Production analysis of real screens
- ✅ Drug resistance studies (enrichment critical)
- ✅ Publication-quality results
- ✅ Fast turnaround needed
- ✅ Reviewer-trusted statistics

**Use v3 (Custom RRA) for:**
- 📚 Learning CRISPR screen analysis
- 🔬 Developing novel statistical methods
- 🎓 Teaching computational biology
- 🧪 Understanding RRA algorithm details

---

## 🔬 Technical Details

### MAGeCK RRA Algorithm

MAGeCK implements the Robust Rank Aggregation (RRA) method (Kolde et al., 2012):

1. **Rank sgRNAs** globally by log fold change (depleted → enriched)
2. **Normalize ranks** to [0,1]: u_i = rank_i / N
3. **Calculate beta scores** for each gene's sgRNAs:
   ```
   beta_i = pbeta(u_i, i, n - i + 1)
   ```
   where n = number of sgRNAs for this gene
4. **RRA score** = minimum beta score
5. **Bonferroni correction**: FDR = min(rho × n, 1.0)

**Key advantage:** Robust to outlier sgRNAs - requires consistent ranking across multiple sgRNAs per gene.

### Negative Binomial Model

MAGeCK uses **generalized linear model** with negative binomial distribution:

```
Count_ij ~ NB(μ_ij, φ)
log(μ_ij) = β_0 + β_1 × Condition_j + offset(log(size_j))
```

Where:
- μ_ij = expected count for gene i in sample j
- φ = dispersion parameter (estimated from data)
- β_1 = log fold change estimate
- size_j = library size factor for sample j

**Variance:** Var(Count) = μ + φμ² (captures overdispersion)

This properly models the mean-variance relationship in count data, unlike t-tests which assume constant variance.

---

## 📚 References

1. **MAGeCK Paper**
   - Li, W., Xu, H., Xiao, T. et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology* 15, 554.
   - DOI: [10.1186/s13059-014-0554-4](https://doi.org/10.1186/s13059-014-0554-4)

2. **Robust Rank Aggregation**
   - Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012). Robust rank aggregation for gene list integration and meta-analysis. *Bioinformatics*, 28(4), 573-580.
   - DOI: [10.1093/bioinformatics/btr709](https://doi.org/10.1093/bioinformatics/btr709)

3. **MAGeCK Documentation**
   - [Official Wiki](https://sourceforge.net/p/mageck/wiki/Home/)
   - [Tutorial](https://sourceforge.net/p/mageck/wiki/demo/)

4. **Negative Binomial Models for RNA-seq**
   - Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

---

## 🛠️ Troubleshooting

### MAGeCK not found

```bash
conda activate py37
mageck -h  # Should show help

# If not found:
pip install --upgrade mageck
```

### Import errors

```bash
# Check Python path
python -c "import sys; print(sys.path)"

# Verify cts_mageck package
python -c "from cts_mageck.config.experiment import ExperimentConfig"
```

### Count matrix format issues

MAGeCK requires:
- Tab-separated format
- Columns: sgRNA, Gene, in_library, [samples...]
- No missing values (use 0 for zero counts)

### Memory issues

For very large screens (>100K sgRNAs):
- Use `--normalization control` (faster)
- Process cell types separately
- Increase system memory allocation

---

## 🎉 Success Criteria

Your analysis is successful if:

✅ **Pipeline completes** without errors (9 comparisons)
✅ **Enriched genes detected** (8-14 per comparison @ FDR<0.05)
✅ **Depleted genes detected** (30-31 per comparison @ FDR<0.05)
✅ **Visualizations generated** (41 plots)
✅ **Effect sizes realistic** (logFC -0.7 to -1.3 for depletion, +0.35 to +0.61 for enrichment)
✅ **QC metrics pass** (>80% reads assigned to cell types)

---

## 📝 Citation

If you use cts-MageCK v0.4.0 in your research, please cite:

**MAGeCK:**
```
Li, W., Xu, H., Xiao, T. et al. MAGeCK enables robust identification
of essential genes from genome-scale CRISPR/Cas9 knockout screens.
Genome Biol 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4
```

**cts-MageCK:**
```
[Your citation details]
```

---

## 📧 Support

For issues, questions, or contributions:
- **GitHub Issues:** [Report bugs or request features]
- **Documentation:** See `RESULTS_SUMMARY.md` for detailed validation

---

## 📄 License

MIT License - see LICENSE file for details

---

**cts-MageCK v0.4.0 - Production-Ready Cell-Type-Specific CRISPR Screen Analysis** 🧬

*Fast demultiplexing + Gold standard MAGeCK RRA = Detect both depletion AND enrichment*

**Version:** 0.4.0
**Last Updated:** 2026-01-31
**Status:** Production Ready ✅
