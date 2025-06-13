#!/usr/bin/env python3

import argparse
import os
import sys

from ase import Atoms
from ase.calculators.mixing import SumCalculator
from ase.io import read, write
from ase.optimize import BFGS
from dftd4.ase import DFTD4
from gpaw import GPAW
from rdkit import Chem
from rdkit.Chem import AllChem


def get_gpaw_calc():
    default_params = {
        "mode": "lcao",
        "basis": "dzp",
        "xc": "PBE",
        "txt": None,
    }
    return GPAW(**default_params)


def get_dftd4_calc():
    default_params = {"method": "pbe"}
    return DFTD4(**default_params)


def get_mol_from_inchi(inchi_str):
    try:
        mol = Chem.MolFromInchi(
            inchi_str, sanitize=True, removeHs=False, treatWarningAsError=True
        )
        if mol is None:
            raise ValueError(f"Error: RDKit could not parse `{inchi_str}'.")

        mol = Chem.rdmolops.AddHs(mol, explicitOnly=False)

        inchi_key = Chem.InchiToInchiKey(inchi_str)
        if not inchi_key:
            inchi_key = "unknown_inchi_key"  # Fallback

        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D  # For reproducibility
        if AllChem.EmbedMolecule(mol, params) == -1:
            # NOTE:
            # This will likely lead to a junk starting point, but use it for now
            # as a fallback
            AllChem.EmbedMolecule(mol, useRandomCoords=True)

        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception:
            # The UFF optimisation failed so we use initial coords as a fallback
            pass

        # RDKit mol to ASE mol
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        try:
            positions = mol.GetConformer().GetPositions()
        except ValueError:
            raise ValueError("Error: RDKit molecule has no conformer after embedding.")

        atoms = Atoms(symbols=symbols, positions=positions)
        atoms.center(vacuum=10)  # Create a cell with 10 Å of space around molecule

        return atoms, inchi_key

    except Exception as e:
        raise ValueError(
            f"Error generating molecule from `{inchi_str}' with RDKit: {e}"
        )


def optimise_mol(atoms, logfile, trajfile):
    gpaw_calc = get_gpaw_calc()
    dftd4_calc = get_dftd4_calc()

    if atoms.get_cell().volume < 1e-9:
        atoms.center(vacuum=10)  # Create a cell with 10 Å of space around molecule

    atoms.calc = SumCalculator([dftd4_calc, gpaw_calc])

    optimiser = BFGS(atoms, logfile=logfile, trajectory=trajfile)
    optimiser.run(fmax=0.05)

    energy = atoms.get_potential_energy()

    return atoms, energy


def main():
    parser = argparse.ArgumentParser(
        description="Command-line tool for ASE calculations using DFTD4 and GPAW."
    )

    mol_input_group = parser.add_mutually_exclusive_group(required=True)
    mol_input_group.add_argument(
        "--inchi",
        type=str,
        help="Input molecule via an InChI string.",
    )
    mol_input_group.add_argument(
        "--xyz",
        type=str,
        help="Input a molecule from an XYZ file.",
    )

    args = parser.parse_args()

    mol_id = "unknown_molecule"
    mol_structure = None

    if args.inchi:
        try:
            mol_structure, mol_id = get_mol_from_inchi(args.inchi)
        except ValueError as e:
            print(e, file=sys.stderr)
            sys.exit(1)
    elif args.xyz:
        try:
            mol_structure = read(args.xyz)
            mol_id = os.path.splitext(os.path.basename(args.xyz))[0]
        except Exception as e:
            print(f"Error reading XYZ file `{args.xyz}': {e}", file=sys.stderr)
            sys.exit(1)

    optimised_mol, calculated_energy = optimise_mol(
        mol_structure, f"{mol_id}.log", f"{mol_id}.trj"
    )

    if optimised_mol is not None:
        final_structure_filename = f"{mol_id}_optimised.xyz"
        try:
            write(final_structure_filename, optimised_mol, format="xyz")
        except Exception as e:
            print(f"Error writing final structure: {e}", file=sys.stderr)
            sys.exit(1)

    if calculated_energy is not None:
        print(f"Energy: {calculated_energy:.5f} eV", file=sys.stdout)
    else:
        print("Error: Final energy could not be determined.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
