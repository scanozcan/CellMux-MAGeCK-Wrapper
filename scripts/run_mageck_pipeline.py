#!/usr/bin/env python3
"""
cts-MageCK v0.4.0 - Pipeline with MAGeCK RRA wrapper.

Complete pipeline:
1. Demultiplex pooled FASTQ files by cell type
2. Count sgRNAs per cell type
3. Run MAGeCK test (RRA) for each comparison
"""

import sys
import argparse
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cts_mageck.config.experiment import ExperimentConfig
from cts_mageck.demux.demux_hash_table import FastDemultiplexer
from cts_mageck.counting.count_sgrnas import count_sgrnas_from_fastq, combine_counts_to_matrix
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
        help='Skip demultiplexing step (use existing demux files)'
    )

    parser.add_argument(
        '--skip-count',
        action='store_true',
        help='Skip counting step (use existing count files)'
    )

    return parser.parse_args()


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

    # ========================================================================
    # STEP 1: DEMULTIPLEXING
    # ========================================================================
    if not args.skip_demux:
        print("\n" + "=" * 80)
        print("STEP 1: DEMULTIPLEXING BY CELL TYPE")
        print("=" * 80)

        # Load libraries
        library_files = {
            ct.name: data_dir / f"{ct.name.lower()}_library.txt"
            for ct in config.cell_types
        }

        # Initialize demuxer
        demuxer = FastDemultiplexer(
            celltype_barcodes={ct.name: ct.barcode for ct in config.cell_types},
            library_files=library_files
        )

        # Demultiplex all samples
        all_samples = [config.control] + config.conditions
        for sample_group in all_samples:
            for rep_idx in range(1, sample_group.replicates + 1):
                sample_name = f"{sample_group.name}_Rep{rep_idx}"

                r1_file = data_dir / 'fastq_files' / f"{sample_name}_R1.fastq"
                r2_file = data_dir / 'fastq_files' / f"{sample_name}_R2.fastq"

                print(f"\nDemultiplexing {sample_name}...")
                demuxer.demultiplex_paired_end(
                    r1_file=str(r1_file),
                    r2_file=str(r2_file),
                    output_dir=str(demux_dir),
                    sample_name=sample_name
                )

        print("\n✓ Demultiplexing complete!")
    else:
        print("\n⊳ Skipping demultiplexing (using existing files)")

    # ========================================================================
    # STEP 2: COUNTING sgRNAs
    # ========================================================================
    if not args.skip_count:
        print("\n" + "=" * 80)
        print("STEP 2: COUNTING sgRNAs PER CELL TYPE")
        print("=" * 80)

        for cell_type in config.cell_types:
            print(f"\n[{cell_type.name}]")

            library_file = data_dir / f"{cell_type.name.lower()}_library.txt"

            # Count for each sample
            count_files = {}
            all_samples = [config.control] + config.conditions

            for sample_group in all_samples:
                for rep_idx in range(1, sample_group.replicates + 1):
                    sample_name = f"{sample_group.name}_Rep{rep_idx}"
                    demux_fastq = demux_dir / f"{sample_name}_{cell_type.name}_R1.fastq"

                    if not demux_fastq.exists():
                        print(f"  Warning: Missing {demux_fastq.name}, skipping...")
                        continue

                    count_file = counts_dir / f"{sample_name}_{cell_type.name}_counts.txt"
                    print(f"  Counting {sample_name}...")

                    count_sgrnas_from_fastq(
                        fastq_file=str(demux_fastq),
                        library_file=str(library_file),
                        output_file=str(count_file)
                    )

                    count_files[sample_name] = str(count_file)

            # Combine into count matrix
            print(f"\n  Combining count matrix for {cell_type.name}...")
            matrix_file = counts_dir / f"{cell_type.name}_count_matrix.txt"

            combine_counts_to_matrix(
                count_files=count_files,
                library_file=str(library_file),
                output_file=str(matrix_file)
            )

        print("\n✓ Counting complete!")
    else:
        print("\n⊳ Skipping counting (using existing files)")

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
    print(f"  Demux: {len(list(demux_dir.glob('*.fastq')))} FASTQ files")
    print(f"  Counts: {len(list(counts_dir.glob('*_count_matrix.txt')))} count matrices")
    print(f"  Results: {len(list(results_dir.glob('*.gene_summary.txt')))} MAGeCK gene summaries")
    print(f"  Results: {len(list(results_dir.glob('*_gene_results.txt')))} cts-MageCK gene results")

    print("\nNext steps:")
    print("  1. Visualize results using post_analysis_visualization.py from v3")
    print("  2. Compare with v3 RRA implementation")


if __name__ == '__main__':
    main()
