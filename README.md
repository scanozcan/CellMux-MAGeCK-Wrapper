# ctm-MageCK: Cell-Type-Multiplexed CRISPR Screen Analysis

**Version:** 1.0.0
**Status:** Production Ready ✅
**License:** MIT

---

## 🎯 Overview

**ctm-MageCK** is a computational pipeline for analyzing **pooled CRISPR knockout screens with cell-type-specific barcoding**. It combines ultra-fast hash-table demultiplexing with the industry-standard MAGeCK statistical framework to identify genes essential for cell viability or conferring resistance across multiple cell types and experimental conditions simultaneously.

### Key Capabilities

- ✅ **Ultra-Fast Demultiplexing** - O(1) hash table algorithm processing 180-190K reads/sec
- ✅ **Dual Statistical Methods** - Both RRA (pairwise) and MLE (multi-condition) analyses
- ✅ **Multi-Condition Support** - Analyze multiple treatments/timepoints in a single experiment
- ✅ **Comprehensive Visualizations** - 62 publication-ready figures (volcano plots, heatmaps, QC metrics)
- ✅ **High Throughput** - Optimized for 10M+ reads/sample with <2% unassigned rate
- ✅ **Production-Grade Performance** - Complete workflow in 5-10 minutes

### Workflow

```
FASTQ Files (Paired-end)
    ↓
[1] Demultiplexing by cell-type barcode
    ↓
Cell-type-specific count matrices
    ↓ ↙             ↘
[2] MAGeCK RRA    [3] MAGeCK MLE
(pairwise)        (multi-condition)
    ↓ ↙             ↘
[4] Visualizations & Reports
    ↓
Publication-ready results
```

---

## 📋 Table of Contents

1. [Quick Start](#-quick-start)
2. [Installation](#-installation)
3. [Input Requirements](#-input-requirements)
4. [Running the Pipeline](#-running-the-pipeline)
5. [Output Files](#-output-files)
6. [Interpretation Guide](#-interpretation-guide)
7. [Example Results](#-example-results)
8. [Advanced Usage](#-advanced-usage)
9. [Citation](#-citation)

---

## 🚀 Quick Start

```bash
# 1. Install dependencies
conda create -n py37 python=3.7
conda activate py37
pip install mageck

conda create -n ctm-mageck python=3.12
conda activate ctm-mageck
pip install numpy pandas matplotlib seaborn scipy pyyaml matplotlib-venn

# 2. Prepare experiment configuration
cp config/experiment_template.yaml config/my_experiment.yaml
# Edit my_experiment.yaml with your sample names and cell types

# 3. Run complete pipeline
python scripts/run_mageck_pipeline_v2.py \
    --config config/my_experiment.yaml \
    --data-dir /path/to/fastq_files \
    --output results_rra

python scripts/run_mageck_mle_pipeline.py \
    --config config/my_experiment.yaml \
    --data-dir /path/to/fastq_files \
    --use-existing-counts results_rra/counts \
    --output results_mle

# 4. Generate visualizations
python scripts/post_analysis_visualization.py \
    --results-dir results_rra/results \
    --config config/my_experiment.yaml \
    --output visualizations_rra

python scripts/post_analysis_visualization_mle.py \
    --results-dir results_mle/results_mle \
    --config config/my_experiment.yaml \
    --output visualizations_mle
```

---

## 🔧 Installation

### Prerequisites

- Python 3.7+ (for MAGeCK compatibility)
- Python 3.12 (for visualization scripts)
- Conda or Miniconda
- 8+ GB RAM (16+ GB recommended for large datasets)
- 50+ GB disk space

### Setup

```bash
# Clone repository
git clone https://github.com/your-org/ctm-mageck.git
cd ctm-mageck

# Create conda environment for MAGeCK (Python 3.7)
conda create -n py37 python=3.7
conda activate py37
pip install mageck

# Verify MAGeCK installation
mageck test --version

# Create separate environment for pipeline scripts (Python 3.12)
conda create -n ctm-mageck python=3.12
conda activate ctm-mageck
pip install numpy pandas matplotlib seaborn scipy pyyaml matplotlib-venn

# Verify installation
python scripts/run_mageck_pipeline_v2.py --help
```

---

## 📥 Input Requirements

### 1. FASTQ Files (Paired-End)

**Read 1 (R1):** Contains sgRNA sequence
**Read 2 (R2):** Contains cell-type barcode

```
Expected structure:
[...primers...]NNNNNNNNNNNNNNNNNNNN[...]BBBBBBBB[...]
               └─ sgRNA (20bp) ─┘      └ Barcode (8bp)
                   Position: 12-31          Position: 22-29
```

**File naming convention:**
```
<SampleName>_R1.fastq
<SampleName>_R2.fastq
```

### 2. sgRNA Library File(s)

Tab-separated, no header:
```
AAAAACCGTATACAGGCCCC	KRT1
TGGAATCCTAACGGTAGGCA	KRT1
TATAAAAGGGCCAAGTACCC	CDH5
...
```

**Format:**
- Column 1: sgRNA sequence (20bp)
- Column 2: Gene symbol

**Options:**
- **Single library:** All cell types use same sgRNAs
- **Cell-type-specific libraries:** Different sgRNAs per cell type
  - Named: `{celltype}_library.txt` (e.g., `keratinocyte_library.txt`)

### 3. Cell-Type Barcode File

CSV format:
```csv
Celltype,Barcode
Keratinocyte,ATGCAGGG
Fibroblast,GTTGCAGC
Endothelial,ATAGCACG
```

### 4. Experiment Configuration (YAML)

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

---

## 🎬 Running the Pipeline

### Step 1: Demultiplexing + RRA Analysis

```bash
python scripts/run_mageck_pipeline_v2.py \
    --config config/experiment.yaml \
    --data-dir /path/to/fastq_files \
    --output results_rra \
    --conda-env py37
```

**Options:**
- `--skip-demux`: Use existing demux count files
- `--normalization {median,total,control}`: Normalization method (default: median)

**Output:** 9 pairwise comparisons (3 cell types × 3 conditions)

### Step 2: MLE Multi-Condition Analysis

```bash
python scripts/run_mageck_mle_pipeline.py \
    --config config/experiment.yaml \
    --data-dir /path/to/fastq_files \
    --use-existing-counts results_rra/counts \
    --output results_mle \
    --conda-env py37
```

**Output:** 3 joint analyses (one per cell type, all conditions vs control)

### Step 3: Generate Visualizations

```bash
# RRA visualizations (41 plots)
python scripts/post_analysis_visualization.py \
    --results-dir results_rra/results \
    --config config/experiment.yaml \
    --output visualizations_rra

# MLE visualizations (21 plots)
python scripts/post_analysis_visualization_mle.py \
    --results-dir results_mle/results_mle \
    --config config/experiment.yaml \
    --output visualizations_mle
```

---

## 📊 Output Files

### RRA Results (`results_rra/`)

```
results_rra/
├── demux/
│   └── {Sample}_count.txt             # Per-sample cell-type counts
├── counts/
│   ├── {CellType}_count_matrix.txt    # MAGeCK-format count matrix
│   └── all_samples_combined.txt       # Combined matrix (all cell types)
└── results/
    ├── Control_vs_{Condition}_{CellType}_gene_results.txt
    └── Control_vs_{Condition}_{CellType}.gene_summary.txt
```

**Gene results columns:**
- `Gene`: Gene symbol
- `logFC`: Log2 fold change (treatment vs control)
- `neg|score`: Depletion score (RRA)
- `neg|p-value`: P-value for depletion
- `neg|fdr`: FDR-corrected p-value (depletion)
- `pos|score`: Enrichment score (RRA)
- `pos|p-value`: P-value for enrichment
- `pos|fdr`: FDR-corrected p-value (enrichment)

### MLE Results (`results_mle/`)

```
results_mle/
└── results_mle/
    ├── {CellType}_MLE_gene_results.txt
    ├── {CellType}_MLE_design_matrix.txt
    └── {CellType}_MLE.gene_summary.txt
```

**Gene results columns (per condition):**
- `Gene`: Gene symbol
- `{Condition}|beta`: Effect size (log scale)
- `{Condition}|z`: Z-score
- `{Condition}|p-value`: P-value
- `{Condition}|fdr`: FDR-corrected p-value

### Visualizations

**RRA (41 plots):**
- 12 volcano plots (depleted & enriched)
- 9 QC plots (count distributions, replicate correlations)
- 12 cross-condition comparisons
- 5 heatmaps (top hits, cell-type specificity)
- 3 cell-type distribution plots

**MLE (21 plots):**
- 12 volcano plots (with dual FDR/p-value thresholds)
- 3 beta score heatmaps
- 3 dotplots (logFC vs FDR)
- 3 cross-cell-type comparison heatmaps

---

## 📖 Interpretation Guide

### Identifying Essential Genes

**RRA Analysis (Pairwise):**
```
Significantly depleted: FDR < 0.05, logFC < -0.5
Significantly enriched: FDR < 0.05, logFC > 0.5
```

**MLE Analysis (Multi-condition):**
```
High confidence: FDR < 0.05
Suggestive: p-value < 0.05 (but FDR > 0.05)
```

### Cell-Type Specificity

**Check heatmaps:**
- Gene depleted in ONE cell type only → Cell-type essential
- Gene depleted in ALL cell types → Common essential
- Gene depleted in one condition only → Condition-specific synthetic lethality

### Quality Control

**Good experiment indicators:**
- Unassigned reads: <5%
- Mean reads/sgRNA: >200
- Replicate correlation (R²): >0.8
- Zero-count sgRNAs: <1%

**Red flags:**
- High unassigned rate (>10%): Barcode mismatch or contamination
- Low coverage (<100 reads/sgRNA): Under-sequenced
- Poor replicate correlation (R²<0.6): Technical variability

---

## 📈 Example Results

The included dataset demonstrates analysis of a multi-condition CRISPR screen across three cell types (Keratinocyte, Fibroblast, Endothelial) with three experimental conditions.

### Detection Summary (FDR < 0.05)

| Cell Type | Condition1 | Condition2 | Condition3 |
|-----------|------------|------------|------------|
| **Keratinocyte** | 121 genes (11.8% depleted, 1.7% enriched) | 120 genes (11.8% depleted, 1.6% enriched) | 121 genes (11.8% depleted, 1.7% enriched) |
| **Fibroblast** | 101 genes (10.2% depleted, 1.0% enriched) | 99 genes (10.1% depleted, 0.9% enriched) | 101 genes (10.3% depleted, 0.9% enriched) |
| **Endothelial** | 109 genes (11.4% depleted, 0.8% enriched) | 109 genes (10.7% depleted, 1.5% enriched) | 109 genes (11.1% depleted, 1.1% enriched) |

### Quality Metrics

- **Coverage:** 1,000-1,800 reads/sgRNA across all cell types
- **Unassigned reads:** ~1.5% (excellent)
- **Zero count sgRNAs:** 0% (perfect)
- **Processing speed:** 180-190K reads/sec demultiplexing

### Top Hits Examples

**Keratinocyte Condition1 - Top Depleted:**
- NFKB1 (logFC=-1.96, FDR=2.6e-04)
- JUNB (logFC=-1.95, FDR=2.6e-04)
- MDM2 (logFC=-1.87, FDR=2.6e-04)

**Fibroblast Condition1 - Top Depleted:**
- MMP2 (logFC=-1.86, FDR=5.5e-04)
- CNN1 (logFC=-1.83, FDR=5.5e-04)
- FN1 (logFC=-1.83, FDR=5.5e-04)

**Endothelial Condition1 - Top Depleted:**
- THY1 (logFC=-2.16, FDR=5.5e-04)
- PECAM1 (logFC=-2.04, FDR=7.0e-04)
- ITGB5 (logFC=-1.92, FDR=5.5e-04)

---

## 🔬 Advanced Usage

### Custom Barcode Positions

Edit demuxer initialization in pipeline scripts:

```python
demuxer = FastDemultiplexer(
    barcode_csv=str(barcode_csv),
    barcode_start=22,       # Adjust if needed
    barcode_length=8,       # Adjust if needed
    grna_start=12,          # Adjust if needed
    grna_length=20,         # Standard
    max_barcode_mismatches=1,
    allow_grna_mismatch=True
)
```

### Alternative Normalization Methods

```bash
# Total count normalization
python scripts/run_mageck_pipeline_v2.py \
    --normalization total \
    ...

# Control-based normalization
python scripts/run_mageck_pipeline_v2.py \
    --normalization control \
    ...
```

### Analyzing Subset of Conditions

Edit `config/experiment.yaml` to include only desired conditions:

```yaml
conditions:
  - name: "Condition1"
    replicates: 3
  # Remove Condition2/3 if not needed
```

---

## 📈 Performance Benchmarks

| Dataset Size | Samples | Total Reads | Runtime | Memory |
|--------------|---------|-------------|---------|--------|
| Small | 6 | 30M | 2-3 min | 4 GB |
| Medium | 12 | 120M | 5-10 min | 8 GB |
| Large | 24 | 240M | 15-20 min | 16 GB |

**Hardware:** Intel Core i7, 16GB RAM, SSD storage

---

## 🆚 Why ctm-MageCK?

### vs. Standard MAGeCK

| Feature | Standard MAGeCK | ctm-MageCK |
|---------|-----------------|------------|
| Cell-type demultiplexing | ❌ Manual | ✅ Automated (O(1) speed) |
| Multi-cell-type support | ❌ Separate runs | ✅ Integrated workflow |
| Visualization suite | ❌ Basic | ✅ 62 publication-ready plots |
| Cross-condition comparison | ❌ Manual | ✅ Automated heatmaps/overlays |
| Cell-type specificity analysis | ❌ None | ✅ Built-in |

### vs. CRISPRCloud/CRISPR-Analyzer

- **Faster:** O(1) demultiplexing vs. O(n) string matching
- **More comprehensive:** Both RRA and MLE in single workflow
- **Better visualization:** Automated multi-panel figures
- **Production-ready:** Handles 10M+ reads/sample routinely

---

## 📚 Citation

If you use ctm-MageCK in your research, please cite:

**ctm-MageCK:**
```
[Your publication here]
```

**MAGeCK (underlying algorithm):**
```
Li W, Xu H, Xiao T, et al. MAGeCK enables robust identification of essential genes
from genome-scale CRISPR/Cas9 knockout screens. Genome Biol. 2014;15(12):554.
doi:10.1186/s13059-014-0554-4
```

**MAGeCK-MLE:**
```
Li W, Köster J, Xu H, et al. Quality control, modeling, and visualization of CRISPR
screens with MAGeCK-VISPR. Genome Biol. 2015;16:281.
doi:10.1186/s13059-015-0843-6
```

---

## 🐛 Troubleshooting

### Issue: High unassigned reads (>10%)

**Causes:**
- Barcode sequences don't match configuration
- Wrong barcode positions in R2
- Sequencing quality issues

**Solutions:**
1. Verify barcodes in `celltype_barcodes.csv`
2. Check barcode positions (R2, positions 22-29)
3. Try allowing 2 mismatches: `max_barcode_mismatches=2`

### Issue: No significant genes detected

**Causes:**
- Under-sequenced (<200 reads/sgRNA)
- No true selection pressure
- Wrong control samples
- FDR threshold too strict

**Solutions:**
1. Check coverage in QC report
2. Verify control samples are correct
3. Try RRA p-value threshold (<0.05) for exploratory analysis
4. Use MLE for more power with multi-condition data

### Issue: Poor replicate correlation

**Causes:**
- Technical variability
- Batch effects
- Sample swaps

**Solutions:**
1. Check sample labels
2. Review QC plots for outliers
3. Consider removing low-quality replicates
4. Increase sequencing depth for future experiments

---

## 🔗 Resources

- **MAGeCK Documentation:** https://sourceforge.net/p/mageck/wiki/Home/
- **MAGeCK GitHub:** https://github.com/liulab-dfci/MAGeCK
- **Paper (MAGeCK):** https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4

---

## 👥 Support

For questions, issues, or feature requests:
- **GitHub Issues:** https://github.com/your-org/ctm-mageck/issues
- **Email:** support@your-institution.edu

---

## 📝 License

MIT License - see LICENSE file for details.

---

**ctm-MageCK v1.0.0** - Production-Ready Cell-Type-Specific CRISPR Screen Analysis
*Ultra-fast demultiplexing + Gold-standard MAGeCK statistics*

**Last Updated:** 2026-02-03
**Maintained by:** [Your Lab/Organization]
