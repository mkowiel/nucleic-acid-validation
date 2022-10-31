from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.validator import Validator

BASIC_SUGAR_BONDS = {
    "sugar_basic==A_G": [
        BondDefinition("sugar_basic==A_G", "C1'", "C2'", 0, 0, 1.530, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C2'", "C3'", 0, 0, 1.528, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C3'", "C4'", 0, 0, 1.525, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C4'", "O4'", 0, 0, 1.450, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C1'", "O4'", 0, 0, 1.413, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C4'", "C5'", 0, 0, 1.508, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C2'", "O2'", 0, 0, 1.412, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==A_G", "C1'", "N9", 0, 0, 1.460, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "sugar_basic==U_T_C": [
        BondDefinition("sugar_basic==U_T_C", "C1'", "C2'", 0, 0, 1.532, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C2'", "C3'", 0, 0, 1.526, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C3'", "C4'", 0, 0, 1.523, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C4'", "O4'", 0, 0, 1.450, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C1'", "O4'", 0, 0, 1.411, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C4'", "C5'", 0, 0, 1.507, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C2'", "O2'", 0, 0, 1.412, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==U_T_C", "C1'", "N1", 0, 0, 1.480, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "sugar_basic==DA_DG": [
        BondDefinition("sugar_basic==DA_DG", "C1'", "C2'", 0, 0, 1.520, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DA_DG", "C2'", "C3'", 0, 0, 1.521, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DA_DG", "C3'", "C4'", 0, 0, 1.528, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DA_DG", "C4'", "O4'", 0, 0, 1.444, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DA_DG", "C1'", "O4'", 0, 0, 1.422, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DA_DG", "C4'", "C5'", 0, 0, 1.509, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DA_DG", "C1'", "N9", 0, 0, 1.456, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "sugar_basic==DU_DT_DC": [
        BondDefinition("sugar_basic==DU_DT_DC", "C1'", "C2'", 0, 0, 1.519, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DU_DT_DC", "C2'", "C3'", 0, 0, 1.518, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DU_DT_DC", "C3'", "C4'", 0, 0, 1.525, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DU_DT_DC", "C4'", "O4'", 0, 0, 1.446, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DU_DT_DC", "C1'", "O4'", 0, 0, 1.417, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DU_DT_DC", "C4'", "C5'", 0, 0, 1.511, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("sugar_basic==DU_DT_DC", "C1'", "N1", 0, 0, 1.476, 0.014, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
}


BASIC_SUGAR_ANGLES = {
    "sugar_basic==A_G": [
        AngleDefinition("sugar_basic==A_G", "C1'", "C2'", "C3'", 0, 0, 0, 101.4, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C2'", "C3'", "C4'", 0, 0, 0, 102.6, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C3'", "C4'", "O4'", 0, 0, 0, 105.6, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C1'", "O4'", "C4'", 0, 0, 0, 109.7, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C2'", "C1'", "O4'", 0, 0, 0, 106.4, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C1'", "C2'", "O2'", 0, 0, 0, 111.2, 2.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C3'", "C2'", "O2'", 0, 0, 0, 113.0, 2.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C2'", "C3'", "O3'", 0, 0, 0, 111.1, 2.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C4'", "C3'", "O3'", 0, 0, 0, 110.3, 2.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C3'", "C4'", "C5'", 0, 0, 0, 115.6, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C5'", "C4'", "O4'", 0, 0, 0, 108.9, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "N9", "C1'", "O4'", 0, 0, 0, 108.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "N9", "C1'", "C2'", 0, 0, 0, 113.9, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C4'", "C5'", "O5'", 0, 0, 0, 111.6, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C1'", "N9", "C4", 0, 0, 0, 126.8, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==A_G", "C1'", "N9", "C8", 0, 0, 0, 127.1, 1.8, 0, 180, 180, 0, 360, 0, 360),
    ],
    "sugar_basic==U_T_C": [
        AngleDefinition("sugar_basic==U_T_C", "C1'", "C2'", "C3'", 0, 0, 0, 101.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C2'", "C3'", "C4'", 0, 0, 0, 102.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C3'", "C4'", "O4'", 0, 0, 0, 104.9, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C1'", "O4'", "C4'", 0, 0, 0, 109.8, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C2'", "C1'", "O4'", 0, 0, 0, 106.6, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C1'", "C2'", "O2'", 0, 0, 0, 110.4, 2.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C3'", "C2'", "O2'", 0, 0, 0, 111.9, 2.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C2'", "C3'", "O3'", 0, 0, 0, 111.3, 3.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C4'", "C3'", "O3'", 0, 0, 0, 111.7, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C3'", "C4'", "C5'", 0, 0, 0, 116.0, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C5'", "C4'", "O4'", 0, 0, 0, 109.6, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "N1", "C1'", "O4'", 0, 0, 0, 108.5, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "N1", "C1'", "C2'", 0, 0, 0, 112.9, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C4'", "C5'", "O5'", 0, 0, 0, 110.7, 2.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C1'", "N1", "C2", 0, 0, 0, 117.5, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==U_T_C", "C1'", "N1", "C6", 0, 0, 0, 121.3, 1.5, 0, 180, 180, 0, 360, 0, 360),
    ],
    "sugar_basic==DA_DG": [
        AngleDefinition("sugar_basic==DA_DG", "C1'", "C2'", "C3'", 0, 0, 0, 102.4, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C3'", "C4'", "O4'", 0, 0, 0, 105.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C1'", "O4'", "C4'", 0, 0, 0, 109.1, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C2'", "C1'", "O4'", 0, 0, 0, 106.0, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C2'", "C3'", "O3'", 0, 0, 0, 110.6, 2.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C4'", "C3'", "O3'", 0, 0, 0, 109.6, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C3'", "C4'", "C5'", 0, 0, 0, 114.7, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C5'", "C4'", "O4'", 0, 0, 0, 108.9, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "N9 ", "C1'", "O4'", 0, 0, 0, 107.9, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "N9 ", "C1'", "C2'", 0, 0, 0, 115.0, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C4'", "C5'", "O5'", 0, 0, 0, 111.0, 2.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C1'", "N9", "C4", 0, 0, 0, 126.9, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DA_DG", "C1'", "N9", "C8", 0, 0, 0, 126.8, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "sugar_basic==DU_DT_DC": [
        AngleDefinition("sugar_basic==DU_DT_DC", "C1'", "C2'", "C3'", 0, 0, 0, 102.6, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C3'", "C4'", "O4'", 0, 0, 0, 105.8, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C1'", "O4'", "C4'", 0, 0, 0, 110.0, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C2'", "C1'", "O4'", 0, 0, 0, 106.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C2'", "C3'", "O3'", 0, 0, 0, 110.1, 2.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C4'", "C3'", "O3'", 0, 0, 0, 110.2, 2.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C3'", "C4'", "C5'", 0, 0, 0, 115.0, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C5'", "C4'", "O4'", 0, 0, 0, 109.6, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "N1", "C1'", "O4'", 0, 0, 0, 107.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "N1", "C1'", "C2'", 0, 0, 0, 113.9, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C4'", "C5'", "O5'", 0, 0, 0, 110.2, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C1'", "N1", "C2", 0, 0, 0, 118.1, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("sugar_basic==DU_DT_DC", "C1'", "N1", "C6", 0, 0, 0, 120.5, 1.3, 0, 180, 180, 0, 360, 0, 360),
    ],
}


class BasicSugarValidator(Validator):
    """
    Validator for nucleotde basees
    """

    # pylint: disable=too-few-public-methods
    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        super().__init__(geometry, csd_sig)

        self.bonds_definition = BASIC_SUGAR_BONDS
        self.angles_definition = BASIC_SUGAR_ANGLES

    def _atom_names_bonds(self, res_name: str) -> List[BondDefinition]:
        if res_name in ("A", "G"):
            return self.bonds_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.bonds_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.bonds_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.bonds_definition["sugar_basic==DU_DT_DC"]
        raise Exception("Non-standard residue")

    def _atom_names_angles(self, res_name: str) -> List[AngleDefinition]:
        if res_name in ("A", "G"):
            return self.angles_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.angles_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.angles_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.angles_definition["sugar_basic==DU_DT_DC"]
        raise Exception("Non-standard residue")

    def _find_bond_definitions(self, res_name: str, altloc: str, atom1_name: str, atom2_name: str) -> List[BondDefinition]:
        if res_name in ("A", "G"):
            return self.bonds_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.bonds_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.bonds_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.bonds_definition["sugar_basic==DU_DT_DC"]
        raise Exception("Non-standard residue")

    def _find_anlge_definitions(
        self, res_name: str, altloc: str, atom1_name: str, atom2_name: str, atom3_name: str
    ) -> List[AngleDefinition]:
        # pylint: disable=too-many-arguments
        if res_name in ("A", "G"):
            return self.angles_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.angles_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.angles_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.angles_definition["sugar_basic==DU_DT_DC"]
        raise Exception("Non-standard residue")
