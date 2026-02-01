#!/usr/bin/env python3
"""
Experiment configuration for cts-MageCK.

Provides YAML-based configuration for multi-condition experiments.
"""

import yaml
from dataclasses import dataclass, field
from typing import List, Dict, Optional
from pathlib import Path


@dataclass
class CellTypeConfig:
    """Configuration for a single cell type."""
    name: str
    barcode: str

    def validate(self) -> List[str]:
        """Validate cell type configuration."""
        errors = []
        if not self.name:
            errors.append("Cell type name cannot be empty")
        if not self.barcode or len(self.barcode) != 8:
            errors.append(f"Barcode must be 8bp, got: {self.barcode}")
        return errors


@dataclass
class ConditionConfig:
    """Configuration for a single experimental condition."""
    name: str
    replicates: int = 3

    def get_sample_names(self) -> List[str]:
        """Get sample names for this condition."""
        return [f"{self.name}_Rep{i+1}" for i in range(self.replicates)]

    def validate(self) -> List[str]:
        """Validate condition configuration."""
        errors = []
        if not self.name:
            errors.append("Condition name cannot be empty")
        if self.replicates < 1:
            errors.append(f"Replicates must be >= 1, got: {self.replicates}")
        return errors


@dataclass
class ControlConfig:
    """Configuration for control samples."""
    name: str = "Control"
    replicates: int = 3

    def get_sample_names(self) -> List[str]:
        """Get sample names for control."""
        return [f"{self.name}_Rep{i+1}" for i in range(self.replicates)]

    def validate(self) -> List[str]:
        """Validate control configuration."""
        errors = []
        if not self.name:
            errors.append("Control name cannot be empty")
        if self.replicates < 1:
            errors.append(f"Replicates must be >= 1, got: {self.replicates}")
        return errors


@dataclass
class Comparison:
    """A single control vs condition comparison."""
    control_name: str
    condition_name: str
    control_samples: List[str]
    condition_samples: List[str]

    @property
    def name(self) -> str:
        """Get comparison name."""
        return f"Control_vs_{self.condition_name}"


@dataclass
class ExperimentConfig:
    """Master experiment configuration."""
    name: str
    control: ControlConfig
    conditions: List[ConditionConfig]
    cell_types: List[CellTypeConfig]
    library_config: Optional[str] = 'default'  # Library preset name

    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'ExperimentConfig':
        """
        Load experiment configuration from YAML file.

        Args:
            yaml_path: Path to YAML configuration file

        Returns:
            ExperimentConfig instance

        Raises:
            FileNotFoundError: If YAML file doesn't exist
            yaml.YAMLError: If YAML is invalid
            ValueError: If configuration is invalid

        Examples:
            >>> config = ExperimentConfig.from_yaml('config/experiment.yaml')
            >>> print(config.name)
            Multi_Condition_Screen
        """
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Config file not found: {yaml_path}")

        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)

        # Parse experiment section
        exp_data = data.get('experiment', {})

        # Parse control
        control_data = exp_data.get('control', {})
        control = ControlConfig(
            name=control_data.get('name', 'Control'),
            replicates=control_data.get('replicates', 3)
        )

        # Parse conditions
        conditions = []
        for cond_data in exp_data.get('conditions', []):
            condition = ConditionConfig(
                name=cond_data.get('name'),
                replicates=cond_data.get('replicates', 3)
            )
            conditions.append(condition)

        # Parse cell types
        cell_types = []
        for ct_data in exp_data.get('cell_types', []):
            cell_type = CellTypeConfig(
                name=ct_data.get('name'),
                barcode=ct_data.get('barcode')
            )
            cell_types.append(cell_type)

        # Create config
        config = cls(
            name=exp_data.get('name', 'Experiment'),
            control=control,
            conditions=conditions,
            cell_types=cell_types,
            library_config=exp_data.get('library_config', 'default')
        )

        # Validate
        errors = config.validate()
        if errors:
            raise ValueError(f"Invalid configuration:\n" + "\n".join(f"  - {e}" for e in errors))

        return config

    def get_all_sample_names(self) -> List[str]:
        """Get all sample names (control + all conditions)."""
        samples = self.control.get_sample_names()
        for condition in self.conditions:
            samples.extend(condition.get_sample_names())
        return samples

    def get_comparisons(self) -> List[Comparison]:
        """
        Generate all Control vs Condition comparisons.

        Returns:
            List of Comparison objects
        """
        comparisons = []
        for condition in self.conditions:
            comparison = Comparison(
                control_name=self.control.name,
                condition_name=condition.name,
                control_samples=self.control.get_sample_names(),
                condition_samples=condition.get_sample_names()
            )
            comparisons.append(comparison)
        return comparisons

    def get_cell_type_barcodes(self) -> Dict[str, str]:
        """Get mapping of cell type to barcode."""
        return {ct.name: ct.barcode for ct in self.cell_types}

    def validate(self) -> List[str]:
        """
        Validate experiment configuration.

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        # Validate name
        if not self.name:
            errors.append("Experiment name cannot be empty")

        # Validate control
        errors.extend([f"Control: {e}" for e in self.control.validate()])

        # Validate conditions
        if not self.conditions:
            errors.append("At least one condition must be specified")

        for i, condition in enumerate(self.conditions):
            cond_errors = condition.validate()
            errors.extend([f"Condition {i+1}: {e}" for e in cond_errors])

        # Check for duplicate condition names
        condition_names = [c.name for c in self.conditions]
        if len(condition_names) != len(set(condition_names)):
            errors.append("Duplicate condition names detected")

        # Validate cell types
        if not self.cell_types:
            errors.append("At least one cell type must be specified")

        for i, cell_type in enumerate(self.cell_types):
            ct_errors = cell_type.validate()
            errors.extend([f"Cell type {i+1}: {e}" for e in ct_errors])

        # Check for duplicate cell type names or barcodes
        ct_names = [ct.name for ct in self.cell_types]
        if len(ct_names) != len(set(ct_names)):
            errors.append("Duplicate cell type names detected")

        ct_barcodes = [ct.barcode for ct in self.cell_types]
        if len(ct_barcodes) != len(set(ct_barcodes)):
            errors.append("Duplicate cell type barcodes detected")

        return errors

    def get_summary(self) -> str:
        """
        Get human-readable summary of experiment configuration.

        Returns:
            Formatted string with experiment details
        """
        summary = []
        summary.append("=" * 60)
        summary.append(f"EXPERIMENT: {self.name}")
        summary.append("=" * 60)
        summary.append("")
        summary.append(f"Control: {self.control.name} ({self.control.replicates} replicates)")
        summary.append("")
        summary.append(f"Conditions ({len(self.conditions)}):")
        for condition in self.conditions:
            summary.append(f"  - {condition.name} ({condition.replicates} replicates)")
        summary.append("")
        summary.append(f"Cell Types ({len(self.cell_types)}):")
        for cell_type in self.cell_types:
            summary.append(f"  - {cell_type.name}: {cell_type.barcode}")
        summary.append("")
        summary.append(f"Total Samples: {len(self.get_all_sample_names())}")
        summary.append(f"Total Comparisons: {len(self.get_comparisons())}")
        summary.append("")
        summary.append("Comparisons:")
        for comparison in self.get_comparisons():
            summary.append(f"  - {comparison.name}")
        summary.append("=" * 60)

        return "\n".join(summary)


def create_template_yaml(output_path: str):
    """
    Create a template experiment YAML file.

    Args:
        output_path: Path where template should be created
    """
    template = {
        'experiment': {
            'name': 'Multi_Condition_Screen',
            'control': {
                'name': 'Control',
                'replicates': 3
            },
            'conditions': [
                {'name': 'Condition1', 'replicates': 3},
                {'name': 'Condition2', 'replicates': 3},
                {'name': 'Condition3', 'replicates': 3}
            ],
            'cell_types': [
                {'name': 'Keratinocyte', 'barcode': 'ATGCAGGG'},
                {'name': 'Fibroblast', 'barcode': 'GTTGCAGC'},
                {'name': 'Endothelial', 'barcode': 'ATAGCACG'}
            ],
            'library_config': 'default'
        }
    }

    with open(output_path, 'w') as f:
        yaml.dump(template, f, default_flow_style=False, sort_keys=False)

    print(f"✓ Created template configuration: {output_path}")


if __name__ == '__main__':
    """Test experiment configuration."""
    import tempfile

    print("Testing Experiment Configuration Module")
    print("=" * 60)

    # Create temporary template
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        temp_path = f.name
        create_template_yaml(temp_path)

    # Load and validate
    print("\nLoading configuration...")
    config = ExperimentConfig.from_yaml(temp_path)

    # Show summary
    print("\n" + config.get_summary())

    # Validate
    print("\nValidation:")
    errors = config.validate()
    if errors:
        print("  Errors found:")
        for error in errors:
            print(f"    ✗ {error}")
    else:
        print("  ✓ Configuration is valid!")

    # Clean up
    import os
    os.unlink(temp_path)

    print("\n" + "=" * 60)
    print("✓ Experiment configuration module working correctly!")
