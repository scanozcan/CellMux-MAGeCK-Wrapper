#!/usr/bin/env python3
"""
cts-MageCK v0.4.0 - Pipeline with MAGeCK RRA wrapper.

Complete pipeline:
1. Demultiplex pooled FASTQ files by cell type (generates counts directly)
2. Combine counts into matrices per cell type
3. Run MAGeCK test (RRA) for each comparison
"""

import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cts_mageck.config.experiment import ExperimentConfig
from cts_mageck.demux.demux_hash_table import FastDemultiplexer
from cts_mageck.wrappers.mageck_wrapper import run_celltype_comparison


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='cts-MageCK v0.4.0 - MAGeCK RRA Pipeline'
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to experiment YAML configuration'
    )

    parser.add_argument(
        '--data-dir',
        type=str,
        required=True,
        help='Directory containing test data (library files + FASTQ)'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='pipeline_outputs_v4',
        help='Output directory (default: pipeline_outputs_v4)'
    )

    parser.add_argument(
        '--conda-env',
        type=str,
        default='py37',
        help='Conda environment with MAGeCK installed (default: py37)'
    )

    parser.add_argument(
        '--normalization',
        type=str,
        default='median',
        choices=['median', 'total', 'control'],
        help='MAGeCK normalization method (default: median)'
    )

    parser.add_argument(
        '--skip-demux',
        action='store_true',
        help='Skip demultiplexing step (use existing count files)'
    )

    return parser.parse_args()


def combine_sample_counts(
    demux_dir: Path,
    sample_names: list,
    cell_types: list,
    output_file: Path
) -> pd.DataFrame:
    """
    Combine demux count files into a single count matrix.

    Args:
        demux_dir: Directory with demux count files
        sample_names: List of sample names
        cell_types: List of cell type names
        output_file: Path to save combined matrix

    Returns:
        Combined count matrix DataFrame
    """
    # Load first sample to get sgRNA/Gene template
    first_sample = sample_names[0]
    first_file = demux_dir / f"{first_sample}_count.txt"

    if not first_file.exists():
        raise FileNotFoundError(f"Missing count file: {first_file}")

    template_df = pd.read_csv(first_file, sep='\t')

    # Start with sgRNA and Gene columns
    matrix_df = template_df[['sgRNA', 'Gene']].copy()
    matrix_df['in_library'] = 1  # All sgRNAs are in library

    # Add counts for each sample
    for sample_name in sample_names:
        count_file = demux_dir / f"{sample_name}_count.txt"

        if not count_file.exists():
            print(f"  Warning: Missing {count_file.name}, using zeros")
            for cell_type in cell_types:
                matrix_df[f"{sample_name}_{cell_type}"] = 0
            continue

        count_df = pd.read_csv(count_file, sep='\t')

        # Add columns for each cell type in this sample
        for cell_type in cell_types:
            if cell_type in count_df.columns:
                matrix_df[f"{sample_name}_{cell_type}"] = count_df[cell_type].values
            else:
                matrix_df[f"{sample_name}_{cell_type}"] = 0

    # Save combined matrix
    matrix_df.to_csv(output_file, sep='\t', index=False)

    print(f"  Combined {len(sample_names)} samples × {len(cell_types)} cell types")
    print(f"  Output: {output_file.name}")

    return matrix_df


def extract_celltype_matrix(
    combined_matrix: pd.DataFrame,
    cell_type: str,
    sample_names: list,
    output_file: Path,
    library_file: Path = None
) -> pd.DataFrame:
    """
    Extract counts for a specific cell type into MAGeCK format.

    IMPORTANT: Filters to only include sgRNAs from this cell type's library,
    removing sgRNAs from other cell types that have all-zero counts.

    Args:
        combined_matrix: Combined count matrix
        cell_type: Cell type to extract
        sample_names: List of sample names
        output_file: Path to save cell-type-specific matrix
        library_file: Optional path to cell-type-specific library file for filtering

    Returns:
        Cell-type-specific count matrix
    """
    # Start with sgRNA, Gene, in_library
    celltype_df = combined_matrix[['sgRNA', 'Gene', 'in_library']].copy()

    # Add sample columns for this cell type
    for sample_name in sample_names:
        col_name = f"{sample_name}_{cell_type}"
        if col_name in combined_matrix.columns:
            celltype_df[sample_name] = combined_matrix[col_name].values
        else:
            celltype_df[sample_name] = 0

    # Set in_library based on whether sgRNA is in this cell type's library
    if library_file is not None and library_file.exists():
        # Load cell-type-specific library (no header, first column is sgRNA)
        lib_df = pd.read_csv(library_file, sep='\t', header=None, names=['sgRNA', 'Gene'])
        lib_sgrnas = set(lib_df['sgRNA'].values)

        # Set in_library: 1 if in this library, 0 otherwise
        celltype_df['in_library'] = celltype_df['sgRNA'].isin(lib_sgrnas).astype(int)

        # Filter to only keep sgRNAs from this library
        celltype_df = celltype_df[celltype_df['in_library'] == 1].copy()

        n_in_library = (celltype_df['in_library'] == 1).sum()
        n_filtered = len(combined_matrix) - len(celltype_df)
        print(f"    Filtered to {n_in_library:,} sgRNAs from library (removed {n_filtered:,} from other cell types)")
    else:
        # Fallback: filter out rows where all counts are zero
        count_cols = [col for col in celltype_df.columns if col not in ['sgRNA', 'Gene', 'in_library']]
        row_sums = celltype_df[count_cols].sum(axis=1)
        celltype_df = celltype_df[row_sums > 0].copy()

        n_filtered = len(combined_matrix) - len(celltype_df)
        print(f"    Filtered to {len(celltype_df):,} sgRNAs with non-zero counts (removed {n_filtered:,} all-zero rows)")

    # Save
    celltype_df.to_csv(output_file, sep='\t', index=False)

    return celltype_df


def main():
    """Main execution."""
    args = parse_args()

    print("=" * 80)
    print("cts-MageCK v0.4.0 - MAGeCK RRA PIPELINE")
    print("=" * 80)

    # Load config
    print(f"\nLoading configuration: {args.config}")
    config = ExperimentConfig.from_yaml(args.config)

    # Setup directories
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True, parents=True)

    demux_dir = output_dir / 'demux'
    counts_dir = output_dir / 'counts'
    results_dir = output_dir / 'results'

    demux_dir.mkdir(exist_ok=True)
    counts_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)

    # Get all sample names
    all_samples = []
    sample_groups = [config.control] + config.conditions
    for sample_group in sample_groups:
        for rep_idx in range(1, sample_group.replicates + 1):
            all_samples.append(f"{sample_group.name}_Rep{rep_idx}")

    # ========================================================================
    # STEP 1: DEMULTIPLEXING (generates count files directly)
    # ========================================================================
    if not args.skip_demux:
        print("\n" + "=" * 80)
        print("STEP 1: DEMULTIPLEXING BY CELL TYPE")
        print("=" * 80)

        # Initialize demuxer
        barcode_csv = data_dir / "celltype_barcodes.csv"
        demuxer = FastDemultiplexer(
            barcode_csv=str(barcode_csv),
            barcode_start=22,
            barcode_length=8,
            grna_start=12,
            grna_length=20,
            max_barcode_mismatches=1,
            allow_grna_mismatch=True
        )

        # Process all samples
        for sample_name in all_samples:
            r1_file = data_dir / 'fastq_files' / f"{sample_name}_R1.fastq"
            r2_file = data_dir / 'fastq_files' / f"{sample_name}_R2.fastq"

            if not r1_file.exists() or not r2_file.exists():
                print(f"\nWarning: Missing FASTQ files for {sample_name}, skipping...")
                continue

            print(f"\nDemultiplexing {sample_name}...")

            # Use combined library (contains all cell types)
            library_file = data_dir / "combined_library.txt"
            output_prefix = demux_dir / sample_name

            demuxer.process(
                read1_file=str(r1_file),
                read2_file=str(r2_file),
                library_file=str(library_file),
                output_prefix=str(output_prefix)
            )

        print("\n✓ Demultiplexing complete!")
    else:
        print("\n⊳ Skipping demultiplexing (using existing files)")

    # ========================================================================
    # STEP 2: COMBINE COUNTS INTO CELL-TYPE MATRICES
    # ========================================================================
    print("\n" + "=" * 80)
    print("STEP 2: CREATING CELL-TYPE-SPECIFIC COUNT MATRICES")
    print("=" * 80)

    cell_type_names = [ct.name for ct in config.cell_types]

    # First, combine all samples into one big matrix
    print("\nCombining all sample counts...")
    combined_matrix_file = counts_dir / "all_samples_combined.txt"

    combined_matrix = combine_sample_counts(
        demux_dir=demux_dir,
        sample_names=all_samples,
        cell_types=cell_type_names,
        output_file=combined_matrix_file
    )

    # Extract per-cell-type matrices for MAGeCK
    print("\nExtracting cell-type-specific matrices...")
    for cell_type in config.cell_types:
        matrix_file = counts_dir / f"{cell_type.name}_count_matrix.txt"
        library_file = data_dir / f"{cell_type.name.lower()}_library.txt"

        print(f"\n  [{cell_type.name}]")
        extract_celltype_matrix(
            combined_matrix=combined_matrix,
            cell_type=cell_type.name,
            sample_names=all_samples,
            output_file=matrix_file,
            library_file=library_file
        )

    print("\n✓ Count matrix creation complete!")

    # ========================================================================
    # STEP 3: RUN MAGeCK TEST
    # ========================================================================
    print("\n" + "=" * 80)
    print("STEP 3: RUNNING MAGeCK RRA ANALYSIS")
    print("=" * 80)

    comparison_count = 0

    for cell_type in config.cell_types:
        print(f"\n[{cell_type.name}]")

        matrix_file = counts_dir / f"{cell_type.name}_count_matrix.txt"

        if not matrix_file.exists():
            print(f"  Error: Count matrix not found: {matrix_file}")
            continue

        # Run comparison for each condition
        for condition in config.conditions:
            comparison_name = f"Control_vs_{condition.name}"

            # Get control and treatment sample names
            control_samples = [
                f"Control_Rep{i}" for i in range(1, config.control.replicates + 1)
            ]
            treatment_samples = [
                f"{condition.name}_Rep{i}" for i in range(1, condition.replicates + 1)
            ]

            output_prefix = results_dir / f"{comparison_name}_{cell_type.name}"

            print(f"\n  Running {comparison_name}...")
            print(f"    Control: {', '.join(control_samples)}")
            print(f"    Treatment: {', '.join(treatment_samples)}")

            try:
                results = run_celltype_comparison(
                    count_matrix_file=str(matrix_file),
                    control_samples=control_samples,
                    treatment_samples=treatment_samples,
                    output_prefix=str(output_prefix),
                    conda_env=args.conda_env,
                    normalization=args.normalization
                )

                # Save parsed gene results in cts-MageCK format
                gene_output = results_dir / f"{comparison_name}_{cell_type.name}_gene_results.txt"
                results['gene_results'].to_csv(gene_output, sep='\t', index=False)

                print(f"    ✓ Gene results saved: {gene_output.name}")

                comparison_count += 1

            except Exception as e:
                print(f"    ERROR: {e}")
                import traceback
                traceback.print_exc()
                continue

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "=" * 80)
    print("✓ PIPELINE COMPLETE!")
    print("=" * 80)

    print(f"\nCompleted {comparison_count} comparisons")
    print(f"Output directory: {output_dir}")

    print("\nGenerated files:")
    print(f"  Demux: {len(list(demux_dir.glob('*_count.txt')))} count files")
    print(f"  Counts: {len(list(counts_dir.glob('*_count_matrix.txt')))} count matrices")
    print(f"  Results: {len(list(results_dir.glob('*.gene_summary.txt')))} MAGeCK gene summaries")
    print(f"  Results: {len(list(results_dir.glob('*_gene_results.txt')))} cts-MageCK gene results")

    print("\nNext steps:")
    print("  python scripts/post_analysis_visualization.py \\")
    print(f"    --results-dir {results_dir} \\")
    print(f"    --config {args.config} \\")
    print(f"    --output visualizations_v4")


if __name__ == '__main__':
    main()
