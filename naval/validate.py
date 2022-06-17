import os
import sys

from Bio.PDB import PDBParser

from naval.nucleotide_geometry import NUCLEOTIDE_RES_NAMES
from naval.nucleotide_geometry import NucleotideGeometry
from naval.bases_validator import BasesValidator


def read_structure(pdb_file_path):
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


def iterate_struct(structure):
    lines = []
    pdbcode = structure.id
    print(f"# PDB id: {pdbcode}")
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                if res_name in NUCLEOTIDE_RES_NAMES:
                    geometry = NucleotideGeometry(pdbcode, model, chain, residue)
                    geometry.calculate_conformation()
                    geometry.prepare_report_torsion()

                    bases = BasesValidator(geometry)
                    lines.extend(bases.validate())
    return lines


def main(filepath):
    sructure = read_structure(filepath)
    iterate_struct(sructure)
    return 0


if __name__ == "__main__":
    pdb_filepath = sys.argv[1]
    main(pdb_filepath)
