import os
import sys
from typing import List

from Bio.PDB import MMCIFParser, PDBParser

from naval.nucleotide_geometry import NucleotideGeometry
from naval.printer import CsvPrinter
from naval.residue_cache_entry import ResidueCacheEntry
from naval.validators.bases_validator import BasesValidator
from naval.validators.po4_validator import Po4Validator
from naval.validators.sugar_validator import BasicSugarValidator


def read_structure(pdb_file_path):
    if pdb_file_path.endswith("pdb"):
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
    elif pdb_file_path.endswith("cif"):
        parser = MMCIFParser(QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


def fill_residue_cache(structure, pdbcode) -> List[ResidueCacheEntry]:
    residue_cache = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_cache.append(ResidueCacheEntry(pdbcode, model, chain, residue))
    return residue_cache


def link_residues(residue_cache) -> List[ResidueCacheEntry]:
    for prev_residue, current_residue in zip(residue_cache[:-1], residue_cache[1:]):
        # TODO: add only when same chain and is not hetatm
        # TODO: move to residue cache entry
        # same chain or next seqid or the same with insetion code
        if current_residue.chain == prev_residue.chain and (
            abs(current_residue.resseq - prev_residue.resseq) == 1
            or (
                current_residue.resseq == prev_residue.resseq
                and (current_residue.inscode != " " or prev_residue.inscode != " ")
            )
        ):
            current_residue.prev_res = prev_residue
            prev_residue.next_res = current_residue
    return residue_cache


def calculate_geometry(residue_cache):
    for residue_entry in residue_cache:
        if residue_entry.is_nucleotide():
            geometry = NucleotideGeometry(residue_entry)
            geometry.calculate_conformation()
            geometry.prepare_report_torsion()
            residue_entry.geometry = geometry
    return residue_cache


def validate_structure(structure):
    pdbcode = structure.id
    print(f"# PDB id: {pdbcode}")

    residue_cache = fill_residue_cache(structure, pdbcode)
    residue_cache = link_residues(residue_cache)
    residue_cache = calculate_geometry(residue_cache)

    records = []
    for residue_entry in residue_cache:
        if residue_entry.is_nucleotide():
            geometry = residue_entry.geometry

            bases = BasesValidator(geometry)
            records.extend(bases.validate())

            po4 = Po4Validator(geometry)
            records.extend(po4.validate())

            sugar = BasicSugarValidator(geometry)
            records.extend(sugar.validate())

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
