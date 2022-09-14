import os
import sys

from Bio.PDB import PDBParser

from naval.nucleotide_geometry import NUCLEOTIDE_RES_NAMES
from naval.nucleotide_geometry import NucleotideGeometry
from naval.validators.bases_validator import BasesValidator
from naval.validators.po4_validator import Po4Validator
from naval.printer import CsvPrinter


def read_structure(pdb_file_path):
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


def validate_structure(structure):
    residue_cache = []
    records = []
    pdbcode = structure.id
    print(f"# PDB id: {pdbcode}")
    for model in structure:
        for chain in model:
            for residue in chain:
                # ResidueCacheEntry
                residue_cache.append(NucleotideGeometry(pdbcode, model, chain, residue))

    for prev_residue, current_residue in zip(residue_cache[1:], residue_cache[:-1]):
        current_residue.prev_res = prev_residue
        prev_residue.next_res = current_residue

    for geometry in residue_cache:
        if geometry.residue.get_resname() in NUCLEOTIDE_RES_NAMES:
            geometry.calculate_conformation()
            geometry.prepare_report_torsion()

            bases = BasesValidator(geometry)
            records.extend(bases.validate())

            po4 = Po4Validator(geometry)
            records.extend(po4.validate())

    return records


def main(filepath):
    sructure = read_structure(filepath)
    validation_records = validate_structure(sructure)
    printer = CsvPrinter()
    lines = printer.print(validation_records)
    # print or save to file
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    PDB_FILEPATH = sys.argv[1]
    main(PDB_FILEPATH)
