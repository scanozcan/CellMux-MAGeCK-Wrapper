#!/usr/bin/env python3
"""
Wrapper for MAGeCK RRA and MLE analysis.

This wrapper interfaces with the official MAGeCK tool (Li et al., Genome Biology 2014)
to perform:
- RRA (Robust Rank Aggregation) - pairwise comparisons
- MLE (Maximum Likelihood Estimation) - multi-condition comparisons

MLE is particularly useful for comparing multiple conditions simultaneously
against a control group.
"""

import subprocess
import pandas as pd
from pathlib import Path
from typing import List, Optional, Dict


class MAGeCKWrapper:
    """
    Wrapper for running MAGeCK RRA analysis.

    Uses the official MAGeCK implementation via conda environment.
    """

    def __init__(self, conda_env: str = "py37"):
        """
        Initialize MAGeCK wrapper.

        Args:
            conda_env: Name of conda environment with MAGeCK installed
        """
        self.conda_env = conda_env

    def run_test(
        self,
        count_matrix: str,
        control_samples: List[str],
        treatment_samples: List[str],
        output_prefix: str,
        normalization: str = "median",
        control_sgrna: Optional[str] = None,
        fdr_threshold: float = 0.25
    ) -> dict:
        """
        Run MAGeCK test (RRA analysis).

        Args:
            count_matrix: Path to count matrix file (tab-separated)
            control_samples: List of control sample column names
            treatment_samples: List of treatment sample column names
            output_prefix: Output file prefix
            normalization: Normalization method ('median', 'total', 'control')
            control_sgrna: Optional path to control sgRNA list for normalization
            fdr_threshold: FDR threshold for calling significance (default: 0.25)

        Returns:
            Dictionary with paths to output files
        """
        # Build MAGeCK command
        cmd = [
            "conda", "run", "-n", self.conda_env,
            "mageck", "test",
            "-k", count_matrix,
            "-t", ",".join(treatment_samples),
            "-c", ",".join(control_samples),
            "-n", output_prefix,
            "--norm-method", normalization,
            "--adjust-method", "fdr"
        ]

        # Add control sgRNA file if provided
        if control_sgrna:
            cmd.extend(["--control-sgrna", control_sgrna])

        # Add FDR threshold
        cmd.extend(["--gene-test-fdr-threshold", str(fdr_threshold)])

        print(f"\nRunning MAGeCK test...")
        print(f"  Command: {' '.join(cmd)}")

        # Run MAGeCK
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            print(result.stdout)

            # Return paths to output files
            output_files = {
                'gene_summary': f"{output_prefix}.gene_summary.txt",
                'sgrna_summary': f"{output_prefix}.sgrna_summary.txt",
                'log': f"{output_prefix}.log"
            }

            return output_files

        except subprocess.CalledProcessError as e:
            print(f"ERROR running MAGeCK:")
            print(f"  Return code: {e.returncode}")
            print(f"  STDOUT: {e.stdout}")
            print(f"  STDERR: {e.stderr}")
            raise

    def parse_gene_results(
        self,
        gene_summary_file: str,
        direction: str = "both"
    ) -> pd.DataFrame:
        """
        Parse MAGeCK gene summary file.

        Args:
            gene_summary_file: Path to .gene_summary.txt file
            direction: Which hits to return ('neg', 'pos', or 'both')

        Returns:
            DataFrame with gene results
        """
        # Read gene summary
        gene_df = pd.read_csv(gene_summary_file, sep='\t')

        # Rename columns to match cts-MageCK v3 format
        gene_df = gene_df.rename(columns={
            'id': 'Gene',
            'num': 'Num_sgRNAs',
            'neg|score': 'RRA_Score_Neg',
            'neg|p-value': 'PValue_Neg',
            'neg|fdr': 'FDR_Neg',
            'neg|rank': 'Rank_Neg',
            'pos|score': 'RRA_Score_Pos',
            'pos|p-value': 'PValue_Pos',
            'pos|fdr': 'FDR_Pos',
            'pos|rank': 'Rank_Pos',
            'neg|lfc': 'LFC_Neg',
            'pos|lfc': 'LFC_Pos'
        })

        # Add unified columns for compatibility
        if direction == "neg":
            gene_df['RRA_Score'] = gene_df['RRA_Score_Neg']
            gene_df['FDR'] = gene_df['FDR_Neg']
            gene_df['logFC'] = gene_df['LFC_Neg']
            gene_df['Direction'] = 'depleted'
        elif direction == "pos":
            gene_df['RRA_Score'] = gene_df['RRA_Score_Pos']
            gene_df['FDR'] = gene_df['FDR_Pos']
            gene_df['logFC'] = gene_df['LFC_Pos']
            gene_df['Direction'] = 'enriched'
        else:  # both
            # Use minimum FDR and corresponding direction
            gene_df['FDR'] = gene_df[['FDR_Neg', 'FDR_Pos']].min(axis=1)
            gene_df['Direction'] = gene_df.apply(
                lambda row: 'depleted' if row['FDR_Neg'] < row['FDR_Pos'] else 'enriched',
                axis=1
            )
            gene_df['logFC'] = gene_df.apply(
                lambda row: row['LFC_Neg'] if row['Direction'] == 'depleted' else row['LFC_Pos'],
                axis=1
            )
            gene_df['RRA_Score'] = gene_df.apply(
                lambda row: row['RRA_Score_Neg'] if row['Direction'] == 'depleted' else row['RRA_Score_Pos'],
                axis=1
            )

        # Add significance column
        gene_df['Significance'] = gene_df['FDR'].apply(
            lambda x: 'significant' if x <= 0.05 else 'not_significant'
        )

        return gene_df

    def parse_sgrna_results(
        self,
        sgrna_summary_file: str
    ) -> pd.DataFrame:
        """
        Parse MAGeCK sgRNA summary file.

        Args:
            sgrna_summary_file: Path to .sgrna_summary.txt file

        Returns:
            DataFrame with sgRNA results
        """
        # Read sgRNA summary
        sgrna_df = pd.read_csv(sgrna_summary_file, sep='\t')

        # Rename columns
        sgrna_df = sgrna_df.rename(columns={
            'sgrna': 'sgRNA',
            'Gene': 'Gene',
            'control_count': 'Control_Mean',
            'treatment_count': 'Treatment_Mean',
            'LFC': 'logFC',
            'control_var': 'Control_Var',
            'treatment_var': 'Treatment_Var',
            'p.low': 'PValue_Low',
            'p.high': 'PValue_High',
            'p.twosided': 'PValue',
            'FDR': 'FDR'
        })

        return sgrna_df

    def generate_design_matrix(
        self,
        control_samples: List[str],
        condition_groups: Dict[str, List[str]],
        output_file: str
    ) -> str:
        """
        Generate design matrix for MAGeCK MLE.

        Args:
            control_samples: List of control sample names
            condition_groups: Dict mapping condition names to sample lists
                             e.g., {'Condition1': ['Condition1_Rep1', ...], ...}
            output_file: Path to save design matrix

        Returns:
            Path to design matrix file

        Example design matrix:
            Samples      baseline  Condition1  Condition2  Condition3
            Control_Rep1     1         0           0           0
            Control_Rep2     1         0           0           0
            Control_Rep3     1         0           0           0
            Cond1_Rep1       1         1           0           0
            Cond1_Rep2       1         1           0           0
            Cond1_Rep3       1         1           0           0
            Cond2_Rep1       1         0           1           0
            ...
        """
        # Collect all samples
        all_samples = control_samples.copy()
        for condition_samples in condition_groups.values():
            all_samples.extend(condition_samples)

        # Create design matrix
        design_data = []

        for sample in all_samples:
            row = {'Samples': sample, 'baseline': 1}

            # Check which condition this sample belongs to
            for condition_name, condition_samples in condition_groups.items():
                if sample in condition_samples:
                    row[condition_name] = 1
                else:
                    row[condition_name] = 0

            design_data.append(row)

        # Create DataFrame
        design_df = pd.DataFrame(design_data)

        # Reorder columns: Samples, baseline, then conditions in order
        column_order = ['Samples', 'baseline'] + list(condition_groups.keys())
        design_df = design_df[column_order]

        # Save
        design_df.to_csv(output_file, sep='\t', index=False)

        print(f"\n  Design matrix saved: {output_file}")
        print(f"  Samples: {len(all_samples)}")
        print(f"  Conditions: {len(condition_groups)} (+ baseline)")

        return output_file

    def run_mle(
        self,
        count_matrix: str,
        design_matrix: str,
        output_prefix: str,
        normalization: str = "median",
        control_sgrna: Optional[str] = None
    ) -> dict:
        """
        Run MAGeCK MLE (Maximum Likelihood Estimation).

        MLE estimates the selection effect of each condition compared to baseline,
        using a negative binomial model with all conditions simultaneously.

        Args:
            count_matrix: Path to count matrix file (tab-separated)
            design_matrix: Path to design matrix file
            output_prefix: Output file prefix
            normalization: Normalization method ('median', 'total', 'control')
            control_sgrna: Optional path to control sgRNA list for normalization

        Returns:
            Dictionary with paths to output files
        """
        # Build MAGeCK command
        cmd = [
            "conda", "run", "-n", self.conda_env,
            "mageck", "mle",
            "-k", count_matrix,
            "-d", design_matrix,
            "-n", output_prefix,
            "--norm-method", normalization
        ]

        # Add control sgRNA file if provided
        if control_sgrna:
            cmd.extend(["--control-sgrna", control_sgrna])

        print(f"\nRunning MAGeCK MLE...")
        print(f"  Command: {' '.join(cmd)}")

        # Run MAGeCK
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            print(result.stdout)

            # Return paths to output files
            output_files = {
                'gene_summary': f"{output_prefix}.gene_summary.txt",
                'sgrna_summary': f"{output_prefix}.sgrna_summary.txt",
                'log': f"{output_prefix}.log"
            }

            return output_files

        except subprocess.CalledProcessError as e:
            print(f"ERROR running MAGeCK MLE:")
            print(f"  Return code: {e.returncode}")
            print(f"  STDOUT: {e.stdout}")
            print(f"  STDERR: {e.stderr}")
            raise

    def parse_mle_gene_results(
        self,
        gene_summary_file: str,
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Parse MAGeCK MLE gene summary file.

        Args:
            gene_summary_file: Path to .gene_summary.txt file
            fdr_threshold: FDR threshold for significance

        Returns:
            DataFrame with gene results for each condition
        """
        # Read gene summary
        gene_df = pd.read_csv(gene_summary_file, sep='\t')

        # MLE output has columns like:
        # Gene, [condition]|beta, [condition]|z, [condition]|p-value, [condition]|fdr
        # for each condition in the design matrix

        # Rename 'Gene' column if it's called something else
        if 'gene' in gene_df.columns:
            gene_df = gene_df.rename(columns={'gene': 'Gene'})

        # Add significance columns for each condition
        fdr_cols = [col for col in gene_df.columns if col.endswith('|fdr')]
        for fdr_col in fdr_cols:
            condition = fdr_col.replace('|fdr', '')
            sig_col = f"{condition}|significant"
            gene_df[sig_col] = gene_df[fdr_col].apply(
                lambda x: 'significant' if x <= fdr_threshold else 'not_significant'
            )

        return gene_df

    def parse_mle_sgrna_results(
        self,
        sgrna_summary_file: str
    ) -> pd.DataFrame:
        """
        Parse MAGeCK MLE sgRNA summary file.

        Args:
            sgrna_summary_file: Path to .sgrna_summary.txt file

        Returns:
            DataFrame with sgRNA results
        """
        # Read sgRNA summary
        sgrna_df = pd.read_csv(sgrna_summary_file, sep='\t')

        return sgrna_df


def run_celltype_comparison(
    count_matrix_file: str,
    control_samples: List[str],
    treatment_samples: List[str],
    output_prefix: str,
    conda_env: str = "py37",
    normalization: str = "median"
) -> dict:
    """
    Convenience function to run MAGeCK test for a single comparison.

    Args:
        count_matrix_file: Path to count matrix
        control_samples: List of control sample names
        treatment_samples: List of treatment sample names
        output_prefix: Output file prefix
        conda_env: Conda environment name
        normalization: Normalization method

    Returns:
        Dictionary with parsed results
    """
    wrapper = MAGeCKWrapper(conda_env=conda_env)

    # Run MAGeCK test
    output_files = wrapper.run_test(
        count_matrix=count_matrix_file,
        control_samples=control_samples,
        treatment_samples=treatment_samples,
        output_prefix=output_prefix,
        normalization=normalization
    )

    # Parse results
    gene_results = wrapper.parse_gene_results(output_files['gene_summary'])
    sgrna_results = wrapper.parse_sgrna_results(output_files['sgrna_summary'])

    return {
        'gene_results': gene_results,
        'sgrna_results': sgrna_results,
        'output_files': output_files
    }


def run_mle_multi_condition(
    count_matrix_file: str,
    control_samples: List[str],
    condition_groups: Dict[str, List[str]],
    output_prefix: str,
    conda_env: str = "py37",
    normalization: str = "median"
) -> dict:
    """
    Run MAGeCK MLE for multiple conditions with auto-generated design matrix.

    This compares all conditions simultaneously against the control group,
    estimating individual effects for each condition.

    Args:
        count_matrix_file: Path to count matrix
        control_samples: List of control sample names
        condition_groups: Dict mapping condition names to sample lists
                         e.g., {'Condition1': ['Condition1_Rep1', ...],
                                'Condition2': ['Condition2_Rep1', ...]}
        output_prefix: Output file prefix
        conda_env: Conda environment name
        normalization: Normalization method

    Returns:
        Dictionary with parsed results and design matrix

    Example:
        >>> results = run_mle_multi_condition(
        ...     count_matrix_file="Keratinocyte_count_matrix.txt",
        ...     control_samples=["Control_Rep1", "Control_Rep2", "Control_Rep3"],
        ...     condition_groups={
        ...         "Condition1": ["Condition1_Rep1", "Condition1_Rep2", "Condition1_Rep3"],
        ...         "Condition2": ["Condition2_Rep1", "Condition2_Rep2", "Condition2_Rep3"],
        ...         "Condition3": ["Condition3_Rep1", "Condition3_Rep2", "Condition3_Rep3"]
        ...     },
        ...     output_prefix="results/Keratinocyte_MLE"
        ... )
    """
    wrapper = MAGeCKWrapper(conda_env=conda_env)

    # Generate design matrix
    design_matrix_file = f"{output_prefix}_design_matrix.txt"
    wrapper.generate_design_matrix(
        control_samples=control_samples,
        condition_groups=condition_groups,
        output_file=design_matrix_file
    )

    # Run MAGeCK MLE
    output_files = wrapper.run_mle(
        count_matrix=count_matrix_file,
        design_matrix=design_matrix_file,
        output_prefix=output_prefix,
        normalization=normalization
    )

    # Parse results
    gene_results = wrapper.parse_mle_gene_results(output_files['gene_summary'])
    sgrna_results = wrapper.parse_mle_sgrna_results(output_files['sgrna_summary'])

    return {
        'gene_results': gene_results,
        'sgrna_results': sgrna_results,
        'output_files': output_files,
        'design_matrix': design_matrix_file
    }
