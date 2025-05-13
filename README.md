# Molecular Structure Optimiser

A Python-based workflow for converting InChI strings and XYZ files into DFT-optimised 3D molecular structures.

## Description

This tool provides a streamlined command-line interface for converting chemical identifiers (InChI strings) or existing molecular geometries (XYZ files) into quantum-mechanically optimised 3D structures. It integrates RDKit for initial structure generation, DFTD4 for dispersion corrections, and GPAW for DFT calculations within the Atomic Simulation Environment (ASE) framework.

The workflow demonstrates a practical pipeline for computational chemistry tasks:

1. Molecular structure generation from identifiers or existing files
2. Initial structure optimisation using force fields
3. DFT-level geometry optimisation with dispersion corrections
4. Output of optimised structures for further computational workflows

## Features

- Direct conversion of InChI strings to 3D structures
- Support for existing XYZ file input
- Initial structure generation and force field pre-optimisation via RDKit
- DFT optimisation using GPAW (PBE functional)
- Dispersion correction via DFTD4
- Standardised output format for integration with other computational workflows

## Dependencies

- Python 3.9+
- ASE (Atomic Simulation Environment)
- GPAW
- DFTD4
- RDKit

## Usage

### Converting from InChI String

```bash
python main.py --inchi "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
```

### Starting from XYZ File

```bash
python main.py --xyz molecule.xyz
```

### Output Files

For each calculation, the following files are generated:

- `[molecule_id]_optimised.xyz`: Optimised structure in XYZ format
- `[molecule_id].log`: Optimisation log
- `[molecule_id].trj`: Trajectory file for visualisation

## Implementation Details

### Structure Generation (RDKit)

The tool uses RDKit to generate 3D structures from InChI strings using the ETKDG v3 algorithm with a fixed random seed for reproducibility. Force field optimisation with UFF is applied as a pre-optimisation step before quantum calculations.

### Quantum Mechanical Calculations (ASE/GPAW/DFTD4)

- DFT calculations use GPAW's LCAO mode with a dzp basis set and PBE functional
- Dispersion corrections are applied using the DFTD4 method
- Calculations are combined using ASE's SumCalculator
- Geometry optimisation uses the BFGS algorithm with a force convergence of 0.05 eV/Å

## Potential Enhancements

### Configuration Flexibility

- Currently uses hard-coded calculator parameters
- Future versions could include command-line options or configuration files for:
  - Functional selection (beyond just PBE)
  - Basis set options
  - Convergence criteria customisation
  - Cell size parameters

### Extended Features

- Vibrational frequency analysis
- Additional property calculations (charges, dipoles, etc.)
- Support for more file formats (both input and output)
- Parallelisation options for larger systems

## Limitations

- Limited to single molecules (no periodic systems)
- Fixed DFT parameters (PBE functional, dzp basis)
- No solvent effects
- Basic error handling for structure generation failures
- Limited to geometry optimisation (no excited states, reaction paths, etc.)

This tool demonstrates foundational computational chemistry workflow concepts and serves as a starting point for more complex implementations.
