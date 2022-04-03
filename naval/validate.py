import os
import sys
from typing import Dict, List
import numpy as np

from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

NUCLEOTIDE_RES_NAMES = ("A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU")
PURINES_RES_NAMES = ("A", "G", "DA", "DG")


def read_structure(pdb_file_path):
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


class NucleotideGeometry:
    """
    Class to keep cache torsion angles for given residue
    """
    
    # pylint: disable=too-many-instance-attributes
    def __init__(self, model, chain, residue) -> None:
        self.model = model
        self.chain = chain
        self.residue = residue
        self.res_name = self.residue.get_resname()
        self.res_full_id = self.residue.get_id()
        self.resseq = self.res_full_id[1]

        self.alpha = None  # O3'(i-1)-P-O5'-C5'
        self.beta = None  # P-O5'-C5'-C4'
        self.gamma = None  # O5'-C5'-C4'-C3'
        self.delta = None  # C5'-C4'-C3'-O3'
        self.epsilon = None  # C4'-C3'-O3'-P(i+1)
        self.zeta = None  # C3'-O3'-P(i+1)-O5'(i+1)
        self.chi = None  # O4'-C1'-N1-C2 or O4'-C1'-N9-C4
        self.tau_max = None  # sugar pucker amplitude
        self.pseudorotation = None  # the phase angle of pseudorotation

    def get_residue_object(self, relative_position: int = 0):
        return self.chain[self.resseq + relative_position]

    def pick_atoms(self, atom_name: str, relative_position: int):
        relative_residue = self.get_residue_object(relative_position)
        atom_group = relative_residue[atom_name]
        return atom_group.disordered_get_list() if atom_group.is_disordered() else [atom_group]

    @staticmethod
    def _round_torsion(atom1, atom2, atom3, atom4):
        torsion = calc_dihedral(atom1.get_vector(), atom2.get_vector(), atom3.get_vector(), atom4.get_vector())
        return round(np.rad2deg(torsion), 1)

    def _calculate_disordered_torsions(
        self, atom_names: List[str], atom_relative_positions: List[int]
    ) -> Dict[str, float]:
        torsions = {}

        atoms1 = self.pick_atoms(atom_names[0], atom_relative_positions[0])
        atoms2 = self.pick_atoms(atom_names[1], atom_relative_positions[1])
        atoms3 = self.pick_atoms(atom_names[2], atom_relative_positions[2])
        atoms4 = self.pick_atoms(atom_names[3], atom_relative_positions[3])

        for atom1 in atoms1:
            for atom2 in atoms2:
                for atom3 in atoms3:
                    for atom4 in atoms4:
                        alt_locs = set(
                            [
                                atom1.get_altloc().strip(),
                                atom2.get_altloc().strip(),
                                atom3.get_altloc().strip(),
                                atom4.get_altloc().strip(),
                            ]
                        )
                        alt_locs.discard("")
                        # if not mixed (for example only "", or only one alternative fonformation "" and "A")
                        if len(alt_locs) <= 1:
                            torsion = self._round_torsion(atom1, atom2, atom3, atom4)
                            alt_loc = ""
                            if len(alt_locs) == 1:
                                alt_loc = alt_locs.pop()
                            torsions[alt_loc] = torsion
        return torsions

    def calculate_torsions(self, atom_names: List[str], atom_relative_positions: List[int]) -> Dict[str, float]:
        if self.residue.is_disordered() == 0:
            atom1 = self.pick_atoms(atom_names[0], atom_relative_positions[0])[0]
            atom2 = self.pick_atoms(atom_names[1], atom_relative_positions[1])[0]
            atom3 = self.pick_atoms(atom_names[2], atom_relative_positions[2])[0]
            atom4 = self.pick_atoms(atom_names[3], atom_relative_positions[3])[0]
            return {"": self._round_torsion(atom1, atom2, atom3, atom4)}

        return self._calculate_disordered_torsions(atom_names, atom_relative_positions)

    def calculate_chi(self):
        atoms_names = ["O4'", "C1'", "N1", "C2"]
        if self.res_name in PURINES_RES_NAMES:
            atoms_names = ["O4'", "C1'", "N9", "C4"]
        self.chi = self.calculate_torsions(atoms_names, (0, 0, 0, 0))

    def calculate_conformation(self):
        self.calculate_chi()


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
                    geometry = NucleotideGeometry(model, chain, residue)
                    geometry.calculate_conformation()

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
                    print(line, geometry.chi)
                    lines.append(line)
    return lines


def main(filepath):
    read_structure(filepath)
    return 0


if __name__ == "__main__":
    pdb_filepath = sys.argv[1]
    main(pdb_filepath)
