import os
import sys
from typing import List, Tuple, Union

from Bio.PDB import MMCIFParser, PDBParser, Structure

from naval.nucleotide_geometry import NucleotideGeometry
from naval.printer import AnglesCsvPrinter, BondsCsvPrinter, GeometryCsvPrinter
from naval.residue_cache_entry import ResidueCacheEntry
from naval.validation_record import TorsionRecord, ValidationRecord
from naval.validators.bases_validator import BasesValidator
from naval.validators.geometry_validator import GeometryValidator
from naval.validators.po4_validator import Po4Validator
from naval.validators.sugar_pucker_validator import SugarPuckerBasedSugarValidator

MAX_RESIDUE_DISTANCE = 2.0


def read_structure(pdb_file_path: str) -> Structure:
    """
    Read and parse pdb/mm-cif nucleotide structure
    """
    if pdb_file_path.endswith("pdb"):
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
    elif pdb_file_path.endswith("cif"):
        parser = MMCIFParser(QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


def fill_residue_cache(structure: Structure, pdbcode: str) -> List[ResidueCacheEntry]:
    """
    Prepare a list of ResidueCacheEntry elements with all residues in the structire
    """
    residue_cache = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_cache.append(ResidueCacheEntry(pdbcode, model, chain, residue))
    return residue_cache


def calc_res_pair_dist(atom_name1, res1, atom_name2, res2):
    "Return minimal distance for a pair of atoms, taking into account the alternative conformations"
    try:
        atom_group1 = res1.residue[atom_name1]
        atoms1 = atom_group1.disordered_get_list() if atom_group1.is_disordered() else [atom_group1]

        atom_group2 = res2.residue[atom_name2]
        atoms2 = atom_group2.disordered_get_list() if atom_group2.is_disordered() else [atom_group2]
    except KeyError:
        return float("inf")

    dists = []
    for atom1 in atoms1:
        for atom2 in atoms2:
            if atom1.get_altloc() == atom2.get_altloc() or atom1.get_altloc() == " " or atom2.get_altloc() == " ":
                dists.append(round(atom2 - atom1, 3))
    if dists:
        return min(dists)
    return float("inf")


def link_residues(residue_cache: List[ResidueCacheEntry]) -> List[ResidueCacheEntry]:
    """
    Link all residures so that it is possible to easily select previous or next residue
    """
    for prev_residue, current_residue in zip(residue_cache[:-1], residue_cache[1:]):
        # TODO: add only when same chain and is not hetatm
        # TODO: move to residue cache entry
        # same chain or next seqid or the same with insetion code
        if current_residue.chain == prev_residue.chain and (
            abs(current_residue.resseq - prev_residue.resseq) == 1
            or (current_residue.resseq == prev_residue.resseq and (current_residue.inscode != " " or prev_residue.inscode != " "))
        ):
            # O3'(prev) is close enough to P(next) or O3'(next) is close enough to P(prev)
            if (
                min(
                    calc_res_pair_dist("O3'", prev_residue, "P", current_residue),
                    calc_res_pair_dist("P", prev_residue, "O3'", current_residue),
                )
                < MAX_RESIDUE_DISTANCE
            ):
                current_residue.prev_res = prev_residue
                prev_residue.next_res = current_residue
    return residue_cache


def calculate_geometry(residue_cache: List[ResidueCacheEntry]) -> List[ResidueCacheEntry]:
    """
    Iterate over all residues and caclulate required torsion angles and pseudorotation for all nucleotides
    """
    for residue_entry in residue_cache:
        if residue_entry.is_nucleotide():
            geometry = NucleotideGeometry(residue_entry)
            geometry.calculate_conformation()
            # geometry.prepare_report_torsion()
            residue_entry.geometry = geometry
    return residue_cache


def validate_structure(structure) -> Tuple[List[ValidationRecord], List[TorsionRecord]]:
    """
    Calculates torsion angles and pass residues through validators.
    """
    pdbcode = structure.id
    print(f"# PDB id: {pdbcode}")

    residue_cache = fill_residue_cache(structure, pdbcode)
    residue_cache = link_residues(residue_cache)
    residue_cache = calculate_geometry(residue_cache)

    validation_records = []
    geometry_records = []
    for residue_entry in residue_cache:
        if residue_entry.is_nucleotide():
            geometry = residue_entry.geometry

            if geometry:
                geometry_validator = GeometryValidator(geometry)
                geometry_records.extend(geometry_validator.validate())

                bases = BasesValidator(geometry)
                validation_records.extend(bases.validate())

                po4 = Po4Validator(geometry)
                validation_records.extend(po4.validate())

                sugar = SugarPuckerBasedSugarValidator(geometry)
                validation_records.extend(sugar.validate())

    return validation_records, geometry_records


def print_records(
    printer: Union[AnglesCsvPrinter, BondsCsvPrinter, GeometryCsvPrinter],
    validation_records: Union[List[ValidationRecord], List[TorsionRecord]],
    out_filename: str,
):
    """
    Save validation records to a file.
    """
    lines = printer.print(validation_records)  # type: ignore
    # save to file
    with open(out_filename, "w", encoding="utf-8") as out_file:
        out_file.write("\n".join(lines))
        out_file.write("\n")


def main(structure_filepath: str, bonds_out_filepath: str, angles_out_filepath: str, geometry_out_path: str):
    sructure = read_structure(structure_filepath)
    validation_records, geometry_records = validate_structure(sructure)

    bonds_printer = BondsCsvPrinter()
    print_records(bonds_printer, validation_records, bonds_out_filepath)

    angles_printer = AnglesCsvPrinter()
    print_records(angles_printer, validation_records, angles_out_filepath)

    geometry_printer = GeometryCsvPrinter()
    print_records(geometry_printer, geometry_records, geometry_out_path)
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage:\n\tvalidate.py struct_in.pdb|struct_in.cif bonds_out.csv angles_out.csv torsion_out.csv")
        sys.exit(-1)
    PDB_OR_MMCIF_FILEPATH = sys.argv[1]
    BONDS_OUT_FILE = sys.argv[2]
    AGLES_OUT_FILE = sys.argv[3]
    GEOMETRY_OUT_PATH = sys.argv[4]
    main(PDB_OR_MMCIF_FILEPATH, BONDS_OUT_FILE, AGLES_OUT_FILE, GEOMETRY_OUT_PATH)
