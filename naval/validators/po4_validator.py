from typing import List

from naval.validators.validator import Validator
from naval.restraint_definition import AngleDefinition
from naval.restraint_definition import BondDefinition
from naval.nucleotide_geometry import NucleotideGeometry


PO4_BONDS = {
    "PO4==AA_0": [
        BondDefinition("PO4==AA_0", "OP1", "P", 1.487, 0.01, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "OP2", "P", 1.483, 0.01, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O3'", "P", 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O5'", "P", 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O3'", "C3'", 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O5'", "C5'", 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AA_1": [
        BondDefinition("PO4==AA_1", "OP1", "P", 1.483, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "OP2", "P", 1.487, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O3'", "P", 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O5'", "P", 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O3'", "C3'", 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O5'", "C5'", 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AA_2": [
        BondDefinition("PO4==AA_2", "OP1", "P", 1.487, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "OP2", "P", 1.483, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O3'", "P", 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O5'", "P", 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O3'", "C3'", 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O5'", "C5'", 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AA_3": [
        BondDefinition("PO4==AA_3", "OP1", "P", 1.483, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "OP2", "P", 1.487, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O3'", "P", 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O5'", "P", 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O3'", "C3'", 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O5'", "C5'", 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_0": [
        BondDefinition("PO4==AS_0", "OP1", "P", 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "OP2", "P", 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O3'", "P", 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O5'", "P", 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O3'", "C3'", 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O5'", "C5'", 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_1": [
        BondDefinition("PO4==AS_1", "OP1", "P", 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "OP2", "P", 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O3'", "P", 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O5'", "P", 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O3'", "C3'", 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O5'", "C5'", 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_2": [
        BondDefinition("PO4==AS_2", "OP1", "P", 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "OP2", "P", 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O3'", "P", 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O5'", "P", 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O3'", "C3'", 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O5'", "C5'", 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_3": [
        BondDefinition("PO4==AS_3", "OP1", "P", 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "OP2", "P", 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O3'", "P", 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O5'", "P", 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O3'", "C3'", 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O5'", "C5'", 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
}


PO4_ANGLES = {
    "PO4==AA_0": [
        AngleDefinition("PO4==AA_0", "OP1", "P", "OP2", 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O3'", 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O5'", 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O3'", 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O5'", 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "O3'", "P", "O5'", 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "P", "O3'", "C3'", 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "P", "O5'", "C5'", 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AA_1": [
        AngleDefinition("PO4==AA_1", "OP1", "P", "OP2", 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O3'", 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O5'", 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O3'", 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O5'", 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "O3'", "P", "O5'", 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "P", "O3'", "C3'", 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "P", "O5'", "C5'", 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AA_2": [
        AngleDefinition("PO4==AA_2", "OP1", "P", "OP2", 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O3'", 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O5'", 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O3'", 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O5'", 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "O3'", "P", "O5'", 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "P", "O3'", "C3'", 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "P", "O5'", "C5'", 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AA_3": [
        AngleDefinition("PO4==AA_3", "OP1", "P", "OP2", 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O3'", 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O5'", 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O3'", 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O5'", 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "O3'", "P", "O5'", 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "P", "O3'", "C3'", 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "P", "O5'", "C5'", 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_0": [
        AngleDefinition("PO4==AS_0", "OP1", "P", "OP2", 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O3'", 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O5'", 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O3'", 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O5'", 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "O3'", "P", "O5'", 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "P", "O3'", "C3'", 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "P", "O5'", "C5'", 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_1": [
        AngleDefinition("PO4==AS_1", "OP1", "P", "OP2", 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O3'", 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O5'", 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O3'", 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O5'", 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "O3'", "P", "O5'", 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "P", "O3'", "C3'", 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "P", "O5'", "C5'", 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_2": [
        AngleDefinition("PO4==AS_2", "OP1", "P", "OP2", 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O3'", 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O5'", 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O3'", 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O5'", 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "O3'", "P", "O5'", 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "P", "O3'", "C3'", 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "P", "O5'", "C5'", 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_3": [
        AngleDefinition("PO4==AS_3", "OP1", "P", "OP2", 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O3'", 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O5'", 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O3'", 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O5'", 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "O3'", "P", "O5'", 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "P", "O3'", "C3'", 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "P", "O5'", "C5'", 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
    ],
}


class Po4Validator(Validator):
    """
    Validator for nucleotde basees
    """

    # pylint: disable=too-few-public-methods
    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        super().__init__(geometry, csd_sig)

        self.bonds_definition = PO4_BONDS
        self.angles_definition = PO4_ANGLES

    def _atom_names_bonds(self, res_name: str) -> List[BondDefinition]:
        return self.bonds_definition["PO4==AS_0"]

    def _atom_names_angles(self, res_name: str) -> List[AngleDefinition]:
        return self.angles_definition["PO4==AS_0"]

    def _find_bond_definitions(self, res_name: str, altloc: str) -> List[BondDefinition]:
        # pylint: disable=too-many-return-statements
        zeta = self.geometry.zeta_conformation.get(altloc, self.geometry.zeta_conformation.get("", None))
        alpha = self.geometry.alpha_conformation.get(altloc, self.geometry.alpha_conformation.get("", None))
        print(altloc, zeta, alpha, self.geometry.zeta, self.geometry.alpha)
        if zeta == "sc-" and alpha == "sc-":
            return self.bonds_definition["PO4==AS_1"]
        if zeta == "sc+" and alpha == "sc+":
            return self.bonds_definition["PO4==AS_3"]
        if zeta == "sc-" and alpha == "ap":
            return self.bonds_definition["PO4==AA_0"]
        if zeta == "ap" and alpha == "sc-":
            return self.bonds_definition["PO4==AA_1"]
        if zeta == "ap" and alpha == "sc+":
            return self.bonds_definition["PO4==AA_2"]
        if zeta == "sc+" and alpha == "ap":
            return self.bonds_definition["PO4==AA_3"]
        return []

    def _find_anlge_definitions(self, res_name: str, altloc: str) -> List[AngleDefinition]:
        # pylint: disable=too-many-return-statements
        zeta = self.geometry.zeta_conformation.get(altloc, self.geometry.zeta_conformation.get("", None))
        alpha = self.geometry.alpha_conformation.get(altloc, self.geometry.alpha_conformation.get("", None))

        if zeta == "sc-" and alpha == "sc-":
            return self.angles_definition["PO4==AS_1"]
        if zeta == "sc+" and alpha == "sc+":
            return self.angles_definition["PO4==AS_3"]
        if zeta == "sc-" and alpha == "ap":
            return self.angles_definition["PO4==AA_0"]
        if zeta == "ap" and alpha == "sc-":
            return self.angles_definition["PO4==AA_1"]
        if zeta == "ap" and alpha == "sc+":
            return self.angles_definition["PO4==AA_2"]
        if zeta == "sc+" and alpha == "ap":
            return self.angles_definition["PO4==AA_3"]
        return []
