import numpy as np

from Bio.PDB.vectors import calc_angle

from naval.validation_record import ValidationRecord
from naval.nucleotide_geometry import NucleotideGeometry


class Validator:
    """
    Base validator class
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        self.geometry = geometry
        self.csd_sig = csd_sig

        self.bonds_definition: dict = {}
        self.angles_definition: dict = {}

    def _validate_bonds(self, res_name, resseq, chain):
        # pylint: disable=too-many-locals

        records = []
        for definition in self.bonds_definition[res_name]:
            atom_group1 = chain[resseq][definition.atom1]
            atom_group2 = chain[resseq][definition.atom2]

            atoms1 = atom_group1.disordered_get_list() if atom_group1.is_disordered() else [atom_group1]
            atoms2 = atom_group2.disordered_get_list() if atom_group2.is_disordered() else [atom_group2]

            for atom1 in atoms1:
                for atom2 in atoms2:
                    if (
                        atom1.get_altloc() == atom2.get_altloc()
                        or atom1.get_altloc() == " "
                        or atom2.get_altloc() == " "
                    ):

                        dist = round(atom2 - atom1, 3)
                        records.append(
                            ValidationRecord(
                                "bond",
                                definition.name,
                                self.geometry,
                                atom1,
                                atom2,
                                None,
                                dist,
                                definition.csd_target,
                                definition.csd_std,
                                definition.pdb_3low,
                                definition.pdb_3high,
                                definition.pdb_4low,
                                definition.pdb_4high,
                            )
                        )
        return records

    def _validate_angles(self, res_name, resseq, chain):
        # pylint: disable=too-many-locals

        records = []
        for definition in self.angles_definition[res_name]:
            atom_group1 = chain[resseq][definition.atom1]
            atom_group2 = chain[resseq][definition.atom2]
            atom_group3 = chain[resseq][definition.atom3]

            atoms1 = atom_group1.disordered_get_list() if atom_group1.is_disordered() else [atom_group1]
            atoms2 = atom_group2.disordered_get_list() if atom_group2.is_disordered() else [atom_group2]
            atoms3 = atom_group3.disordered_get_list() if atom_group3.is_disordered() else [atom_group3]

            for atom1 in atoms1:
                for atom2 in atoms2:
                    for atom3 in atoms3:
                        if (
                            # pylint: disable=too-many-boolean-expressions
                            (atom1.get_altloc() == atom2.get_altloc() and atom2.get_altloc() == atom3.get_altloc())
                            or (atom1.get_altloc() == " " and atom2.get_altloc() == atom3.get_altloc())
                            or (atom2.get_altloc() == " " and atom1.get_altloc() == atom3.get_altloc())
                            or (atom3.get_altloc() == " " and atom1.get_altloc() == atom2.get_altloc())
                            or (atom1.get_altloc() != " " and atom2.get_altloc() == " " and atom3.get_altloc() == " ")
                            or (atom1.get_altloc() == " " and atom2.get_altloc() != " " and atom3.get_altloc() == " ")
                            or (atom1.get_altloc() == " " and atom2.get_altloc() == " " and atom3.get_altloc() != " ")
                        ):
                            angle_value = calc_angle(
                                atom1.get_vector(),
                                atom2.get_vector(),
                                atom3.get_vector(),
                            )
                            angle_value = np.round(np.rad2deg(angle_value), 1)

                            records.append(
                                ValidationRecord(
                                    "angle",
                                    definition.name,
                                    self.geometry,
                                    atom1,
                                    atom2,
                                    atom3,
                                    angle_value,
                                    definition.csd_target,
                                    definition.csd_std,
                                    definition.pdb_3low,
                                    definition.pdb_3high,
                                    definition.pdb_4low,
                                    definition.pdb_4high,
                                )
                            )

        return records

    def validate(self):
        res_name = self.geometry.res_name
        resseq = self.geometry.resseq
        chain = self.geometry.chain

        records = []
        records.extend(self._validate_bonds(res_name, resseq, chain))
        records.extend(self._validate_angles(res_name, resseq, chain))
        return records
