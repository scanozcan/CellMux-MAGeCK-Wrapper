#!/usr/bin/env python3
"""
Test MAGeCK installation and conda environment.
"""

import subprocess
import sys


def test_mageck_conda(conda_env="py37"):
    """Test if MAGeCK is accessible via conda environment."""
    print("=" * 70)
    print("TESTING MAGeCK INSTALLATION")
    print("=" * 70)

    print(f"\nConda environment: {conda_env}")

    # Test 1: Check if conda is available
    print("\n[1/3] Checking conda installation...")
    try:
        result = subprocess.run(
            ["conda", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"  ✓ Conda installed: {result.stdout.strip()}")
    except FileNotFoundError:
        print("  ✗ ERROR: conda command not found")
        print("    Please install Anaconda or Miniconda")
        return False
    except subprocess.CalledProcessError as e:
        print(f"  ✗ ERROR: {e}")
        return False

    # Test 2: Check if conda environment exists
    print(f"\n[2/3] Checking if conda environment '{conda_env}' exists...")
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True,
            check=True
        )
        if conda_env in result.stdout:
            print(f"  ✓ Environment '{conda_env}' exists")
        else:
            print(f"  ✗ ERROR: Environment '{conda_env}' not found")
            print(f"\n  Create it with:")
            print(f"    conda create -n {conda_env} python=3.7")
            print(f"    conda activate {conda_env}")
            print(f"    pip install mageck")
            return False
    except subprocess.CalledProcessError as e:
        print(f"  ✗ ERROR: {e}")
        return False

    # Test 3: Check if MAGeCK is installed in the environment
    print(f"\n[3/3] Checking MAGeCK installation in '{conda_env}'...")
    try:
        result = subprocess.run(
            ["conda", "run", "-n", conda_env, "mageck", "-h"],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"  ✓ MAGeCK is installed and accessible")

        # Extract MAGeCK version
        for line in result.stdout.split('\n'):
            if 'MAGeCK' in line or 'Version' in line:
                print(f"    {line.strip()}")
                break

        return True

    except FileNotFoundError:
        print(f"  ✗ ERROR: MAGeCK not found in '{conda_env}'")
        print(f"\n  Install it with:")
        print(f"    conda activate {conda_env}")
        print(f"    pip install mageck")
        return False
    except subprocess.CalledProcessError as e:
        print(f"  ✗ ERROR: {e.stderr}")
        return False


def main():
    """Run tests."""
    conda_env = "py37"

    if len(sys.argv) > 1:
        conda_env = sys.argv[1]

    success = test_mageck_conda(conda_env)

    print("\n" + "=" * 70)
    if success:
        print("✓ ALL TESTS PASSED")
        print("=" * 70)
        print("\nMAGeCK is ready to use!")
        print(f"\nYou can now run the pipeline:")
        print(f"  python scripts/run_mageck_pipeline.py \\")
        print(f"    --config config/experiment_template.yaml \\")
        print(f"    --data-dir test_data \\")
        print(f"    --output pipeline_outputs_v4 \\")
        print(f"    --conda-env {conda_env}")
    else:
        print("✗ TESTS FAILED")
        print("=" * 70)
        print("\nPlease fix the issues above before running the pipeline.")

    print()
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
