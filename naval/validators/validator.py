from typing import List

import numpy as np

from Bio.PDB.vectors import calc_angle
from Bio.PDB import Chain

from naval.validation_record import ValidationRecord
from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition
from naval.restraint_definition import BondDefinition


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

    def _atom_names_bonds(self, res_name: str) -> List[BondDefinition]:
        # TODO return list of (d.atom1, d.atom2)
        return self.bonds_definition[res_name]

    def _atom_names_angles(self, res_name: str) -> List[AngleDefinition]:
        return self.angles_definition[res_name]

    def _find_bond_definitions(self, res_name: str, altloc: str) -> List[BondDefinition]:
        return self.bonds_definition[res_name]

    def _find_anlge_definitions(self, res_name: str, altloc: str) -> List[AngleDefinition]:
        return self.angles_definition[res_name]

    def _select_bond_definition(self, definitions, atom1: str, atom2: str) -> BondDefinition:
        for definition in definitions:
            if definition.atom1 == atom1 and definition.atom2 == atom2:
                return definition
        return None

    def _select_angle_definition(self, definitions, atom1: str, atom2: str, atom3: str) -> AngleDefinition:
        for definition in definitions:
            if definition.atom1 == atom1 and definition.atom2 == atom2 and definition.atom3 == atom3:
                return definition
        return None

    def _validate_bonds(self, res_name: str, resseq: str, chain: Chain) -> List[ValidationRecord]:
        # pylint: disable=too-many-locals

        records = []
        for atom_definition in self._atom_names_bonds(res_name):
            if atom_definition.atom1 in chain[resseq] and atom_definition.atom2 in chain[resseq]:
                atom_group1 = chain[resseq][atom_definition.atom1]
                atom_group2 = chain[resseq][atom_definition.atom2]

                atoms1 = atom_group1.disordered_get_list() if atom_group1.is_disordered() else [atom_group1]
                atoms2 = atom_group2.disordered_get_list() if atom_group2.is_disordered() else [atom_group2]

                for atom1 in atoms1:
                    for atom2 in atoms2:
                        if (
                            atom1.get_altloc() == atom2.get_altloc()
                            or atom1.get_altloc() == " "
                            or atom2.get_altloc() == " "
                        ):
                            altloc_set = set([atom1.get_altloc(), atom2.get_altloc()])
                            altloc_set.discard(" ")
                            altloc = "" if len(altloc_set) == 0 else altloc_set.pop()

                            definitions = self._find_bond_definitions(res_name, altloc)
                            definition = self._select_bond_definition(definitions, atom1.name, atom2.name)

                            dist = round(atom2 - atom1, 3)

                            if definition:
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

    def _validate_angles(self, res_name: str, resseq: str, chain: Chain) -> List[ValidationRecord]:
        # pylint: disable=too-many-locals

        records = []
        for atom_definition in self._atom_names_angles(res_name):
            if (
                atom_definition.atom1 in chain[resseq]
                and atom_definition.atom2 in chain[resseq]
                and atom_definition.atom3 in chain[resseq]
            ):
                atom_group1 = chain[resseq][atom_definition.atom1]
                atom_group2 = chain[resseq][atom_definition.atom2]
                atom_group3 = chain[resseq][atom_definition.atom3]

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
                                or (
                                    atom1.get_altloc() != " "
                                    and atom2.get_altloc() == " "
                                    and atom3.get_altloc() == " "
                                )
                                or (
                                    atom1.get_altloc() == " "
                                    and atom2.get_altloc() != " "
                                    and atom3.get_altloc() == " "
                                )
                                or (
                                    atom1.get_altloc() == " "
                                    and atom2.get_altloc() == " "
                                    and atom3.get_altloc() != " "
                                )
                            ):
                                altloc_set = set([atom1.get_altloc(), atom2.get_altloc(), atom3.get_altloc()])
                                altloc_set.discard(" ")
                                altloc = "" if len(altloc_set) == 0 else altloc_set.pop()

                                definitions = self._find_anlge_definitions(res_name, altloc)
                                definition = self._select_angle_definition(
                                    definitions, atom1.name, atom2.name, atom3.name
                                )

                                angle_value = calc_angle(
                                    atom1.get_vector(),
                                    atom2.get_vector(),
                                    atom3.get_vector(),
                                )
                                angle_value = np.round(np.rad2deg(angle_value), 1)

                                if definition:
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

    def validate(self) -> List[ValidationRecord]:
        res_name = self.geometry.res_name
        resseq = self.geometry.resseq
        chain = self.geometry.chain

        records = []
        records.extend(self._validate_bonds(res_name, resseq, chain))
        records.extend(self._validate_angles(res_name, resseq, chain))
        return records
