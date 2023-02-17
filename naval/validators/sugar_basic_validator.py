from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.validator import NonStandardResidueException, Validator

BASIC_SUGAR_BONDS = {
    "sugar_basic==A_G": [
        BondDefinition("sugar_basic==A_G", "C1'", "C2'", 0, 0, 1.530, 0.011, 1494018, 1.527, 0.006, 1.499, 1.546, 1.476, 1.564),
        BondDefinition("sugar_basic==A_G", "C2'", "C3'", 0, 0, 1.528, 0.009, 1494012, 1.522, 0.006, 1.494, 1.547, 1.464, 1.568),
        BondDefinition("sugar_basic==A_G", "C3'", "C4'", 0, 0, 1.525, 0.008, 1494012, 1.519, 0.007, 1.489, 1.544, 1.467, 1.564),
        BondDefinition("sugar_basic==A_G", "C4'", "O4'", 0, 0, 1.450, 0.009, 1494018, 1.450, 0.005, 1.424, 1.468, 1.405, 1.494),
        BondDefinition("sugar_basic==A_G", "C1'", "O4'", 0, 0, 1.413, 0.008, 1494020, 1.412, 0.006, 1.386, 1.433, 1.366, 1.472),
        BondDefinition("sugar_basic==A_G", "C4'", "C5'", 0, 0, 1.508, 0.010, 1493990, 1.507, 0.007, 1.477, 1.532, 1.458, 1.560),
        BondDefinition("sugar_basic==A_G", "C2'", "O2'", 0, 0, 1.412, 0.009, 1493809, 1.417, 0.005, 1.393, 1.434, 1.375, 1.464),
        BondDefinition("sugar_basic==A_G", "C1'", "N9", 0, 0, 1.460, 0.012, 1493933, 1.471, 0.008, 1.438, 1.497, 1.393, 1.525),
    ],
    "sugar_basic==U_T_C": [
        BondDefinition("sugar_basic==U_T_C", "C1'", "C2'", 0, 0, 1.532, 0.009, 1127586, 1.528, 0.005, 1.501, 1.546, 1.481, 1.564),
        BondDefinition("sugar_basic==U_T_C", "C2'", "C3'", 0, 0, 1.526, 0.009, 1127585, 1.522, 0.006, 1.494, 1.546, 1.472, 1.566),
        BondDefinition("sugar_basic==U_T_C", "C3'", "C4'", 0, 0, 1.523, 0.011, 1127590, 1.519, 0.007, 1.490, 1.543, 1.464, 1.566),
        BondDefinition("sugar_basic==U_T_C", "C4'", "O4'", 0, 0, 1.450, 0.008, 1127581, 1.450, 0.005, 1.424, 1.468, 1.403, 1.484),
        BondDefinition("sugar_basic==U_T_C", "C1'", "O4'", 0, 0, 1.411, 0.008, 1127587, 1.412, 0.005, 1.387, 1.434, 1.373, 1.451),
        BondDefinition("sugar_basic==U_T_C", "C4'", "C5'", 0, 0, 1.507, 0.008, 1127572, 1.507, 0.007, 1.477, 1.531, 1.457, 1.558),
        BondDefinition("sugar_basic==U_T_C", "C2'", "O2'", 0, 0, 1.412, 0.009, 1127523, 1.417, 0.005, 1.395, 1.434, 1.377, 1.455),
        BondDefinition("sugar_basic==U_T_C", "C1'", "N1", 0, 0, 1.480, 0.013, 1127476, 1.479, 0.011, 1.439, 1.549, 1.416, 1.587),
    ],
    "sugar_basic==DA_DG": [
        BondDefinition("sugar_basic==DA_DG", "C1'", "C2'", 0, 0, 1.520, 0.010, 101019, 1.522, 0.007, 1.491, 1.553, 1.459, 1.589),
        BondDefinition("sugar_basic==DA_DG", "C2'", "C3'", 0, 0, 1.521, 0.010, 101013, 1.526, 0.008, 1.495, 1.557, 1.446, 1.606),
        BondDefinition("sugar_basic==DA_DG", "C3'", "C4'", 0, 0, 1.528, 0.009, 101006, 1.526, 0.007, 1.492, 1.555, 1.462, 1.577),
        BondDefinition("sugar_basic==DA_DG", "C4'", "O4'", 0, 0, 1.444, 0.008, 101026, 1.449, 0.007, 1.418, 1.477, 1.388, 1.521),
        BondDefinition("sugar_basic==DA_DG", "C1'", "O4'", 0, 0, 1.422, 0.011, 101035, 1.414, 0.009, 1.382, 1.456, 1.352, 1.487),
        BondDefinition("sugar_basic==DA_DG", "C4'", "C5'", 0, 0, 1.509, 0.008, 100940, 1.515, 0.007, 1.483, 1.545, 1.434, 1.581),
        BondDefinition("sugar_basic==DA_DG", "C1'", "N9", 0, 0, 1.456, 0.009, 101009, 1.461, 0.007, 1.426, 1.493, 1.389, 1.528),
    ],
    "sugar_basic==DU_DT_DC": [
        BondDefinition("sugar_basic==DU_DT_DC", "C1'", "C2'", 0, 0, 1.519, 0.010, 100654, 1.522, 0.008, 1.491, 1.554, 1.464, 1.625),
        BondDefinition("sugar_basic==DU_DT_DC", "C2'", "C3'", 0, 0, 1.518, 0.010, 100645, 1.527, 0.009, 1.494, 1.559, 1.226, 1.609),
        BondDefinition("sugar_basic==DU_DT_DC", "C3'", "C4'", 0, 0, 1.525, 0.011, 100648, 1.527, 0.008, 1.492, 1.557, 1.464, 1.644),
        BondDefinition("sugar_basic==DU_DT_DC", "C4'", "O4'", 0, 0, 1.446, 0.010, 100659, 1.447, 0.007, 1.415, 1.477, 1.288, 1.522),
        BondDefinition("sugar_basic==DU_DT_DC", "C1'", "O4'", 0, 0, 1.417, 0.011, 100665, 1.412, 0.009, 1.380, 1.454, 1.286, 1.508),
        BondDefinition("sugar_basic==DU_DT_DC", "C4'", "C5'", 0, 0, 1.511, 0.011, 100565, 1.514, 0.007, 1.482, 1.545, 1.415, 1.583),
        BondDefinition("sugar_basic==DU_DT_DC", "C1'", "N1", 0, 0, 1.476, 0.014, 100552, 1.488, 0.017, 1.429, 1.558, 1.345, 1.623),
    ],
}


BASIC_SUGAR_ANGLES = {
    "sugar_basic==A_G": [
        AngleDefinition("sugar_basic==A_G", "C1'", "C2'", "C3'", 0, 0, 0, 101.4, 1.1, 1494011, 101.4, 0.8, 98.5, 105.4, 94.9, 106.8),
        AngleDefinition("sugar_basic==A_G", "C2'", "C3'", "C4'", 0, 0, 0, 102.6, 0.8, 1494009, 102.5, 0.8, 99.0, 106.2, 95.9, 108.3),
        AngleDefinition("sugar_basic==A_G", "C3'", "C4'", "O4'", 0, 0, 0, 105.6, 1.1, 1494010, 104.3, 1.0, 101.2, 107.6, 97.6, 109.1),
        AngleDefinition("sugar_basic==A_G", "C1'", "O4'", "C4'", 0, 0, 0, 109.7, 0.9, 1494017, 109.8, 0.5, 107.0, 111.8, 101.0, 112.8),
        AngleDefinition("sugar_basic==A_G", "C2'", "C1'", "O4'", 0, 0, 0, 106.4, 1.0, 1494017, 107.3, 0.8, 103.6, 109.5, 100.3, 111.3),
        AngleDefinition("sugar_basic==A_G", "C1'", "C2'", "O2'", 0, 0, 0, 111.2, 2.7, 1493809, 109.0, 1.5, 104.4, 115.3, 96.8, 120.6),
        AngleDefinition("sugar_basic==A_G", "C3'", "C2'", "O2'", 0, 0, 0, 113.0, 2.5, 1493802, 111.4, 1.7, 106.4, 118.2, 101.4, 123.9),
        AngleDefinition("sugar_basic==A_G", "C2'", "C3'", "O3'", 0, 0, 0, 111.1, 2.8, 1493985, 112.9, 1.9, 105.0, 119.7, 97.6, 129.5),
        AngleDefinition("sugar_basic==A_G", "C4'", "C3'", "O3'", 0, 0, 0, 110.3, 2.4, 1493984, 112.1, 1.8, 104.6, 117.8, 97.9, 124.3),
        AngleDefinition("sugar_basic==A_G", "C3'", "C4'", "C5'", 0, 0, 0, 115.6, 1.3, 1493983, 115.8, 1.1, 110.3, 120.1, 105.5, 124.7),
        AngleDefinition("sugar_basic==A_G", "C5'", "C4'", "O4'", 0, 0, 0, 108.9, 1.1, 1493989, 109.7, 0.9, 105.5, 113.3, 102.3, 116.8),
        AngleDefinition("sugar_basic==A_G", "N9", "C1'", "O4'", 0, 0, 0, 108.2, 1.0, 1493930, 108.5, 1.1, 104.2, 113.7, 101.2, 118.5),
        AngleDefinition("sugar_basic==A_G", "N9", "C1'", "C2'", 0, 0, 0, 113.9, 1.4, 1493928, 112.5, 1.7, 107.0, 119.9, 103.0, 126.2),
        AngleDefinition("sugar_basic==A_G", "C4'", "C5'", "O5'", 0, 0, 0, 111.6, 1.9, 1493962, 111.4, 1.1, 107.0, 115.4, 103.7, 118.5),
        AngleDefinition("sugar_basic==A_G", "C1'", "N9", "C4", 0, 0, 0, 126.8, 1.8, 1493908, 126.7, 2.0, 119.3, 135.0, 112.9, 142.9),
        AngleDefinition("sugar_basic==A_G", "C1'", "N9", "C8", 0, 0, 0, 127.1, 1.8, 1493885, 127.1, 1.9, 119.4, 134.0, 112.4, 142.1),
    ],
    "sugar_basic==U_T_C": [
        AngleDefinition("sugar_basic==U_T_C", "C1'", "C2'", "C3'", 0, 0, 0, 101.2, 0.8, 1127585, 101.4, 0.7, 98.8, 105.2, 96.1, 106.6),
        AngleDefinition("sugar_basic==U_T_C", "C2'", "C3'", "C4'", 0, 0, 0, 102.2, 1.0, 1127578, 102.5, 0.8, 99.2, 106.0, 96.1, 107.8),
        AngleDefinition(
            "sugar_basic==U_T_C", "C3'", "C4'", "O4'", 0, 0, 0, 104.9, 1.3, 1127579, 104.2, 1.0, 101.3, 107.5, 98.0, 108.7
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C1'", "O4'", "C4'", 0, 0, 0, 109.8, 0.7, 1127580, 109.8, 0.4, 107.4, 111.6, 103.8, 112.7
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C2'", "C1'", "O4'", 0, 0, 0, 106.6, 0.9, 1127585, 107.4, 0.8, 104.0, 109.3, 101.3, 110.8
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C1'", "C2'", "O2'", 0, 0, 0, 110.4, 2.7, 1127522, 108.9, 1.4, 104.4, 114.9, 97.3, 119.8
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C3'", "C2'", "O2'", 0, 0, 0, 111.9, 2.9, 1127522, 111.3, 1.6, 106.5, 118.3, 102.0, 123.2
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C2'", "C3'", "O3'", 0, 0, 0, 111.3, 3.3, 1127567, 113.0, 1.8, 105.5, 119.1, 99.5, 128.0
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C4'", "C3'", "O3'", 0, 0, 0, 111.7, 2.3, 1127572, 112.4, 1.7, 105.2, 117.7, 99.4, 124.1
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C3'", "C4'", "C5'", 0, 0, 0, 116.0, 1.4, 1127568, 115.8, 1.0, 110.4, 119.8, 104.4, 124.0
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "C5'", "C4'", "O4'", 0, 0, 0, 109.6, 0.9, 1127558, 109.8, 0.8, 105.8, 113.0, 102.2, 116.6
        ),
        AngleDefinition(
            "sugar_basic==U_T_C", "N1", "C1'", "O4'", 0, 0, 0, 108.5, 0.8, 1127447, 108.8, 1.0, 104.8, 113.5, 100.7, 117.8
        ),
        AngleDefinition("sugar_basic==U_T_C", "N1", "C1'", "C2'", 0, 0, 0, 112.9, 1.3, 1127446, 112.6, 1.5, 106.1, 119.7, 96.1, 125.5),
        AngleDefinition(
            "sugar_basic==U_T_C", "C4'", "C5'", "O5'", 0, 0, 0, 110.7, 2.0, 1127493, 111.4, 1.0, 107.2, 115.3, 104.1, 118.1
        ),
        AngleDefinition("sugar_basic==U_T_C", "C1'", "N1", "C2", 0, 0, 0, 117.5, 1.5, 1127454, 118.5, 2.0, 111.3, 127.1, 105.0, 132.0),
        AngleDefinition("sugar_basic==U_T_C", "C1'", "N1", "C6", 0, 0, 0, 121.3, 1.5, 1127451, 121.1, 1.7, 113.8, 127.5, 108.4, 133.0),
    ],
    "sugar_basic==DA_DG": [
        AngleDefinition("sugar_basic==DA_DG", "C1'", "C2'", "C3'", 0, 0, 0, 102.4, 1.4, 101012, 101.9, 1.6, 95.9, 106.6, 91.2, 109.8),
        AngleDefinition("sugar_basic==DA_DG", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 1.1, 101005, 103.3, 1.1, 98.2, 106.8, 93.4, 110.3),
        AngleDefinition("sugar_basic==DA_DG", "C3'", "C4'", "O4'", 0, 0, 0, 105.6, 1.2, 101006, 105.7, 0.9, 100.7, 108.6, 97.4, 113.1),
        AngleDefinition("sugar_basic==DA_DG", "C1'", "O4'", "C4'", 0, 0, 0, 109.1, 1.5, 101024, 109.4, 1.3, 103.2, 112.7, 96.6, 115.0),
        AngleDefinition("sugar_basic==DA_DG", "C2'", "C1'", "O4'", 0, 0, 0, 106.0, 1.0, 101019, 106.0, 1.2, 101.4, 109.6, 97.2, 112.3),
        AngleDefinition("sugar_basic==DA_DG", "C2'", "C3'", "O3'", 0, 0, 0, 110.6, 2.5, 100887, 111.8, 2.0, 101.5, 121.4, 93.5, 128.0),
        AngleDefinition("sugar_basic==DA_DG", "C4'", "C3'", "O3'", 0, 0, 0, 109.6, 1.8, 100881, 109.2, 2.0, 100.5, 118.8, 92.2, 126.3),
        AngleDefinition("sugar_basic==DA_DG", "C3'", "C4'", "C5'", 0, 0, 0, 114.7, 1.3, 100922, 114.9, 1.3, 108.1, 120.7, 98.5, 127.3),
        AngleDefinition("sugar_basic==DA_DG", "C5'", "C4'", "O4'", 0, 0, 0, 108.9, 1.2, 100940, 109.7, 1.1, 103.9, 116.3, 99.8, 121.2),
        AngleDefinition("sugar_basic==DA_DG", "N9 ", "C1'", "O4'", 0, 0, 0, 107.9, 1.1, 101005, 108.3, 1.4, 101.3, 115.3, 95.0, 122.6),
        AngleDefinition("sugar_basic==DA_DG", "N9 ", "C1'", "C2'", 0, 0, 0, 115.0, 1.0, 100989, 114.5, 1.7, 106.7, 121.1, 98.1, 126.9),
        AngleDefinition("sugar_basic==DA_DG", "C4'", "C5'", "O5'", 0, 0, 0, 111.0, 2.2, 100772, 110.5, 1.3, 104.1, 116.4, 98.5, 124.8),
        AngleDefinition("sugar_basic==DA_DG", "C1'", "N9", "C4", 0, 0, 0, 126.9, 1.9, 100988, 126.5, 1.0, 121.8, 130.4, 116.1, 137.2),
        AngleDefinition("sugar_basic==DA_DG", "C1'", "N9", "C8", 0, 0, 0, 126.8, 1.9, 100988, 127.4, 0.9, 123.7, 132.0, 116.5, 136.6),
    ],
    "sugar_basic==DU_DT_DC": [
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C1'", "C2'", "C3'", 0, 0, 0, 102.6, 1.1, 100641, 102.1, 1.7, 96.2, 106.8, 92.5, 108.4
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 1.0, 100642, 103.4, 1.2, 98.2, 107.1, 94.6, 109.4
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C3'", "C4'", "O4'", 0, 0, 0, 105.8, 0.9, 100645, 105.7, 1.0, 100.5, 108.5, 95.2, 110.4
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C1'", "O4'", "C4'", 0, 0, 0, 110.0, 0.7, 100654, 109.0, 1.4, 102.9, 112.4, 97.9, 115.1
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C2'", "C1'", "O4'", 0, 0, 0, 106.3, 0.8, 100651, 105.8, 1.2, 101.0, 109.6, 96.9, 111.6
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C2'", "C3'", "O3'", 0, 0, 0, 110.1, 2.6, 100512, 111.8, 1.9, 101.8, 121.0, 92.9, 127.5
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C4'", "C3'", "O3'", 0, 0, 0, 110.2, 2.1, 100514, 109.3, 1.9, 100.9, 118.7, 95.5, 127.9
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C3'", "C4'", "C5'", 0, 0, 0, 115.0, 1.4, 100549, 114.8, 1.3, 108.3, 120.9, 102.5, 126.4
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C5'", "C4'", "O4'", 0, 0, 0, 109.6, 1.1, 100556, 109.7, 1.2, 104.1, 115.7, 98.6, 122.2
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "N1", "C1'", "O4'", 0, 0, 0, 107.9, 0.7, 100546, 108.5, 1.4, 101.6, 114.7, 94.6, 120.0
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "N1", "C1'", "C2'", 0, 0, 0, 113.9, 1.1, 100535, 114.7, 1.8, 107.3, 121.8, 101.6, 130.5
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C4'", "C5'", "O5'", 0, 0, 0, 110.2, 1.9, 100398, 110.6, 1.3, 104.6, 116.9, 99.1, 123.8
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C1'", "N1", "C2", 0, 0, 0, 118.1, 1.5, 100521, 119.0, 1.0, 114.9, 123.2, 108.6, 130.3
        ),
        AngleDefinition(
            "sugar_basic==DU_DT_DC", "C1'", "N1", "C6", 0, 0, 0, 120.5, 1.3, 100521, 120.1, 0.9, 116.1, 124.0, 108.7, 133.5
        ),
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
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

    def _atom_names_angles(self, res_name: str) -> List[AngleDefinition]:
        if res_name in ("A", "G"):
            return self.angles_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.angles_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.angles_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.angles_definition["sugar_basic==DU_DT_DC"]
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

    def _find_bond_definitions(self, res_name: str, altloc: str, atom1_name: str, atom2_name: str) -> List[BondDefinition]:
        if res_name in ("A", "G"):
            return self.bonds_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.bonds_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.bonds_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.bonds_definition["sugar_basic==DU_DT_DC"]
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

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
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")
