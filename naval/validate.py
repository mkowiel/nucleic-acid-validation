import os
import sys
import numpy as np

from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

NUCLEOTIDE_RES_NAMES = ["A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU"]


def read_structure(pdb_file_path):
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


def pick_atoms(atom_group):
    return atom_group.disordered_get_list() if atom_group.is_disordered() else [atom_group]


def calculate_torsions(chain, resseq, atoms_names):
    torsions = []

    res_obj = chain[resseq]

    atoms1 = pick_atoms(res_obj[atoms_names[0]])
    atoms2 = pick_atoms(res_obj[atoms_names[1]])
    atoms3 = pick_atoms(res_obj[atoms_names[2]])
    atoms4 = pick_atoms(res_obj[atoms_names[3]])

    for atom1 in atoms1:
        for atom2 in atoms2:
            for atom3 in atoms3:
                for atom4 in atoms4:
                    torsion = calc_dihedral(
                        atom1.get_vector(), atom2.get_vector(), atom3.get_vector(), atom4.get_vector()
                    )
                    torsion = round(np.rad2deg(torsion), 1)
                    torsions.append(
                        (
                            chain.get_id(),
                            resseq,
                            atom1.get_name(),
                            atom1.get_altloc(),
                            atom2.get_name(),
                            atom2.get_altloc(),
                            atom3.get_name(),
                            atom3.get_altloc(),
                            atom4.get_name(),
                            atom4.get_altloc(),
                            torsion,
                        )
                    )
    return torsions


def calculate_ksi(chain, resseq, res_name):
    atoms_names = ["O4'", "C1'", "N1", "C2"]
    if res_name in ["A", "G", "DA", "DG"]:
        atoms_names = ["O4'", "C1'", "N9", "C4"]
    return calculate_torsions(chain, resseq, atoms_names)


def iterate_struct(structure):
    lines = []
    pdbcode = structure.id

    for model in structure:
        model_id = model.get_id()
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_name = residue.get_resname()
                if res_name in NUCLEOTIDE_RES_NAMES:
                    res_full_id = residue.get_id()

                    resseq = res_full_id[1]
                    inscode = res_full_id[2]

                    line = ",".join(
                        str(_)
                        for _ in (
                            pdbcode,
                            model_id,
                            chain_id,
                            res_name,
                            resseq,
                            inscode,
                        )
                    )
                    print(line)
                    print(calculate_ksi(chain, resseq, res_name))
                    lines.append(line)
    return lines


def main(filepath):
    read_structure(filepath)
    return 0


if __name__ == "__main__":
    pdb_filepath = sys.argv[1]
    main(pdb_filepath)
