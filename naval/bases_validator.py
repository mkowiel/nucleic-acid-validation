import numpy as np

from Bio.PDB.vectors import calc_angle

from naval.nucleotide_geometry import NucleotideGeometry
from naval.validation_record import ValidationRecord


class BondDefinition:
    """
    Simple container class for bond definitions
    """

    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-instance-attributes

    __slots__ = (
        "atom1",
        "atom2",
        "csd_target",
        "csd_std",
        "pdb_count",
        "pdb_mean",
        "pdb_std",
        "pdb_3low",
        "pdb_3high",
        "pdb_4low",
        "pdb_4high",
    )

    def __init__(
        self,
        atom1: str,
        atom2: str,
        csd_target: float,
        csd_std: float,
        pdb_count: int,
        pdb_mean: float,
        pdb_std: float,
        pdb_3low: float,
        pdb_3high: float,
        pdb_4low: float,
        pdb_4high: float,
    ):
        # pylint: disable=too-many-arguments
        self.atom1 = atom1
        self.atom2 = atom2
        self.csd_target = csd_target
        self.csd_std = csd_std
        self.pdb_count = pdb_count
        self.pdb_mean = pdb_mean
        self.pdb_std = pdb_std
        self.pdb_3low = pdb_3low
        self.pdb_3high = pdb_3high
        self.pdb_4low = pdb_4low
        self.pdb_4high = pdb_4high


class AngleDefinition:
    """
    Simple container class for angle definitions
    """

    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-instance-attributes

    __slots__ = (
        "atom1",
        "atom2",
        "atom3",
        "csd_target",
        "csd_std",
        "pdb_count",
        "pdb_mean",
        "pdb_std",
        "pdb_3low",
        "pdb_3high",
        "pdb_4low",
        "pdb_4high",
    )

    def __init__(
        self,
        atom1: str,
        atom2: str,
        atom3: str,
        csd_target: float,
        csd_std: float,
        pdb_count: int,
        pdb_mean: float,
        pdb_std: float,
        pdb_3low: float,
        pdb_3high: float,
        pdb_4low: float,
        pdb_4high: float,
    ):
        # pylint: disable=too-many-arguments
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.csd_target = csd_target
        self.csd_std = csd_std
        self.pdb_count = pdb_count
        self.pdb_mean = pdb_mean
        self.pdb_std = pdb_std
        self.pdb_3low = pdb_3low
        self.pdb_3high = pdb_3high
        self.pdb_4low = pdb_4low
        self.pdb_4high = pdb_4high


BASES_BONDS = {
    "A": [
        BondDefinition("N1", "C2", 1.339, 0.007, 544647, 1.338, 0.006, 1.307, 1.367, 1.278, 1.410),
        BondDefinition("C2", "N3", 1.330, 0.007, 544637, 1.331, 0.006, 1.301, 1.355, 1.274, 1.411),
        BondDefinition("N3", "C4", 1.346, 0.006, 544636, 1.343, 0.008, 1.299, 1.372, 1.264, 1.410),
        BondDefinition("C4", "C5", 1.382, 0.008, 544646, 1.382, 0.008, 1.344, 1.414, 1.310, 1.466),
        BondDefinition("C5", "C6", 1.406, 0.008, 544647, 1.404, 0.010, 1.351, 1.438, 1.302, 1.483),
        BondDefinition("C6", "N1", 1.353, 0.007, 544647, 1.349, 0.007, 1.310, 1.378, 1.279, 1.414),
        BondDefinition("C5", "N7", 1.388, 0.007, 544647, 1.385, 0.007, 1.348, 1.412, 1.316, 1.503),
        BondDefinition("N7", "C8", 1.311, 0.007, 544631, 1.310, 0.005, 1.285, 1.333, 1.264, 1.380),
        BondDefinition("C8", "N9", 1.370, 0.008, 544630, 1.371, 0.007, 1.334, 1.405, 1.285, 1.454),
        BondDefinition("N9", "C4", 1.374, 0.007, 544645, 1.373, 0.009, 1.321, 1.407, 1.273, 1.454),
        BondDefinition("C6", "N6", 1.334, 0.007, 544645, 1.334, 0.006, 1.301, 1.360, 1.270, 1.483),
    ],
    "G": [
        BondDefinition("N1", "C2", 1.372, 0.006, 782271, 1.372, 0.007, 1.334, 1.400, 1.302, 1.433),
        BondDefinition("C2", "N3", 1.327, 0.005, 782270, 1.322, 0.007, 1.288, 1.353, 1.260, 1.393),
        BondDefinition("N3", "C4", 1.352, 0.006, 782270, 1.349, 0.007, 1.311, 1.377, 1.278, 1.402),
        BondDefinition("C4", "C5", 1.379, 0.006, 782273, 1.377, 0.007, 1.342, 1.404, 1.312, 1.447),
        BondDefinition("C5", "C6", 1.418, 0.008, 782271, 1.417, 0.009, 1.373, 1.451, 1.333, 1.488),
        BondDefinition("C6", "N1", 1.392, 0.006, 782271, 1.390, 0.007, 1.353, 1.418, 1.314, 1.450),
        BondDefinition("C5", "N7", 1.388, 0.006, 782274, 1.387, 0.006, 1.353, 1.411, 1.323, 1.447),
        BondDefinition("N7", "C8", 1.308, 0.006, 782274, 1.304, 0.005, 1.280, 1.326, 1.257, 1.388),
        BondDefinition("C8", "N9", 1.375, 0.006, 782272, 1.372, 0.007, 1.337, 1.399, 1.305, 1.444),
        BondDefinition("N9", "C4", 1.374, 0.006, 782271, 1.374, 0.009, 1.333, 1.407, 1.287, 1.438),
        BondDefinition("C6", "O6", 1.238, 0.007, 782242, 1.237, 0.007, 1.205, 1.271, 1.177, 1.321),
        BondDefinition("C2", "N2", 1.338, 0.007, 782269, 1.339, 0.006, 1.308, 1.365, 1.280, 1.398),
    ],
    "U": [
        BondDefinition("N1", "C2", 1.381, 0.009, 387309, 1.381, 0.011, 1.332, 1.426, 1.292, 1.478),
        BondDefinition("C2", "N3", 1.373, 0.008, 387309, 1.372, 0.008, 1.332, 1.406, 1.299, 1.454),
        BondDefinition("N3", "C4", 1.381, 0.008, 387310, 1.378, 0.009, 1.335, 1.412, 1.301, 1.455),
        BondDefinition("C4", "C5", 1.432, 0.008, 387310, 1.430, 0.008, 1.394, 1.462, 1.362, 1.510),
        BondDefinition("C5", "C6", 1.337, 0.008, 387310, 1.336, 0.006, 1.306, 1.360, 1.283, 1.473),
        BondDefinition("C6", "N1", 1.374, 0.008, 387310, 1.374, 0.009, 1.330, 1.405, 1.290, 1.457),
        BondDefinition("C2", "O2", 1.219, 0.008, 387309, 1.218, 0.007, 1.178, 1.250, 1.146, 1.287),
        BondDefinition("C4", "O4", 1.231, 0.008, 387306, 1.231, 0.007, 1.195, 1.267, 1.172, 1.339),
    ],
    "T": [
        BondDefinition("N1", "C2", 1.376, 0.008, 41088, 1.378, 0.007, 1.353, 1.421, 1.336, 1.488),
        BondDefinition("C2", "N3", 1.372, 0.007, 41088, 1.371, 0.005, 1.339, 1.394, 1.305, 1.465),
        BondDefinition("N3", "C4", 1.382, 0.008, 41088, 1.382, 0.005, 1.347, 1.405, 1.312, 1.437),
        BondDefinition("C4", "C5", 1.446, 0.008, 41088, 1.445, 0.006, 1.418, 1.484, 1.402, 1.525),
        BondDefinition("C5", "C6", 1.340, 0.007, 41088, 1.342, 0.005, 1.319, 1.379, 1.279, 1.463),
        BondDefinition("C6", "N1", 1.381, 0.007, 41088, 1.381, 0.005, 1.353, 1.407, 1.315, 1.450),
        BondDefinition("C2", "O2", 1.222, 0.008, 41086, 1.219, 0.005, 1.192, 1.248, 1.152, 1.274),
        BondDefinition("C4", "O4", 1.229, 0.008, 41086, 1.228, 0.005, 1.203, 1.255, 1.163, 1.300),
        BondDefinition("C7", "C5", 1.498, 0.006, 40843, 1.499, 0.004, 1.472, 1.524, 1.446, 1.559),
    ],
    "C": [
        BondDefinition("N1", "C2", 1.395, 0.009, 596005, 1.398, 0.010, 1.354, 1.436, 1.313, 1.485),
        BondDefinition("C2", "N3", 1.353, 0.007, 596005, 1.353, 0.007, 1.320, 1.381, 1.289, 1.417),
        BondDefinition("N3", "C4", 1.337, 0.008, 596005, 1.332, 0.007, 1.295, 1.360, 1.261, 1.391),
        BondDefinition("C4", "C5", 1.424, 0.010, 596004, 1.423, 0.007, 1.388, 1.449, 1.364, 1.477),
        BondDefinition("C5", "C6", 1.338, 0.008, 596003, 1.338, 0.005, 1.311, 1.358, 1.291, 1.387),
        BondDefinition("C6", "N1", 1.365, 0.007, 596005, 1.366, 0.007, 1.327, 1.392, 1.292, 1.418),
        BondDefinition("C2", "O2", 1.240, 0.008, 596004, 1.240, 0.007, 1.205, 1.270, 1.177, 1.302),
        BondDefinition("C4", "N4", 1.330, 0.008, 596001, 1.334, 0.006, 1.300, 1.357, 1.273, 1.388),
    ],
}
BASES_BONDS["DA"] = BASES_BONDS["A"]
BASES_BONDS["DG"] = BASES_BONDS["G"]
BASES_BONDS["DU"] = BASES_BONDS["U"]
BASES_BONDS["DT"] = BASES_BONDS["T"]
BASES_BONDS["DC"] = BASES_BONDS["C"]


BASES_ANGLES = {
    "A": [
        AngleDefinition("C6", "N1", "C2", 118.6, 0.6, 544647, 118.5, 0.7, 114.2, 121.8, 109.1, 126.0),
        AngleDefinition("N1", "C2", "N3", 129.4, 0.7, 544637, 129.3, 0.7, 126.4, 133.8, 118.2, 137.6),
        AngleDefinition("C2", "N3", "C4", 110.5, 0.6, 544636, 110.6, 0.9, 104.4, 114.1, 98.9, 121.0),
        AngleDefinition("N3", "C4", "C5", 126.9, 0.6, 544636, 126.8, 0.9, 123.1, 131.6, 117.6, 138.4),
        AngleDefinition("C4", "C5", "C6", 117.1, 0.5, 544646, 117.0, 0.7, 114.5, 120.4, 111.9, 124.9),
        AngleDefinition("C5", "C6", "N1", 117.5, 0.5, 544647, 117.7, 0.7, 113.8, 121.2, 109.4, 125.9),
        AngleDefinition("N3", "C4", "N9", 127.2, 0.7, 544635, 127.4, 0.9, 122.8, 131.1, 115.0, 136.4),
        AngleDefinition("C6", "C5", "N7", 132.2, 0.6, 544647, 132.1, 0.9, 126.6, 135.4, 119.3, 139.2),
        AngleDefinition("C5", "C4", "N9", 105.9, 0.4, 544645, 105.8, 0.6, 103.0, 108.6, 100.6, 111.5),
        AngleDefinition("C4", "N9", "C8", 105.7, 0.4, 544628, 105.8, 0.8, 102.0, 109.2, 98.9, 117.9),
        AngleDefinition("N9", "C8", "N7", 113.9, 0.5, 544630, 113.9, 0.8, 110.5, 118.5, 105.8, 122.5),
        AngleDefinition("C8", "N7", "C5", 103.8, 0.4, 544630, 103.7, 0.7, 98.6, 106.7, 93.5, 109.6),
        AngleDefinition("N7", "C5", "C4", 110.6, 0.5, 544646, 110.8, 0.7, 107.9, 115.0, 104.1, 119.3),
        AngleDefinition("N6", "C6", "N1", 118.6, 0.7, 544645, 118.7, 1.2, 113.5, 125.0, 108.2, 131.9),
        AngleDefinition("N6", "C6", "C5", 123.9, 0.7, 544645, 123.6, 1.1, 117.6, 128.5, 111.8, 133.5),
    ],
    "G": [
        AngleDefinition("C6", "N1", "C2", 125.4, 0.5, 782270, 125.0, 0.6, 121.5, 127.7, 117.2, 131.6),
        AngleDefinition("N1", "C2", "N3", 123.6, 0.5, 782270, 124.1, 0.7, 121.3, 127.8, 118.7, 132.1),
        AngleDefinition("C2", "N3", "C4", 112.0, 0.4, 782270, 111.9, 0.7, 108.0, 114.9, 102.2, 118.9),
        AngleDefinition("N3", "C4", "C5", 128.6, 0.5, 782270, 128.5, 0.8, 125.1, 132.1, 120.7, 139.1),
        AngleDefinition("C4", "C5", "C6", 118.9, 0.4, 782270, 119.0, 0.7, 116.3, 122.3, 112.9, 125.7),
        AngleDefinition("C5", "C6", "N1", 111.5, 0.5, 782271, 111.5, 0.7, 107.9, 114.9, 102.7, 118.7),
        AngleDefinition("N3", "C4", "N9", 125.8, 0.7, 782268, 126.1, 1.0, 121.9, 130.1, 112.5, 133.6),
        AngleDefinition("C6", "C5", "N7", 130.3, 0.5, 782271, 130.2, 0.9, 125.7, 133.5, 120.8, 136.3),
        AngleDefinition("C5", "C4", "N9", 105.6, 0.5, 782271, 105.4, 0.6, 102.6, 108.1, 100.5, 110.8),
        AngleDefinition("C4", "N9", "C8", 106.2, 0.4, 782271, 106.3, 0.7, 103.0, 109.7, 99.7, 112.5),
        AngleDefinition("N9", "C8", "N7", 113.2, 0.4, 782272, 113.2, 0.7, 110.0, 116.7, 107.5, 122.4),
        AngleDefinition("C8", "N7", "C5", 104.2, 0.4, 782274, 104.3, 0.6, 100.6, 107.1, 93.1, 110.3),
        AngleDefinition("N7", "C5", "C4", 110.8, 0.4, 782273, 110.8, 0.6, 108.1, 114.0, 104.9, 118.4),
        AngleDefinition("O6", "C6", "N1", 120.1, 0.5, 782242, 120.0, 1.1, 114.9, 125.8, 110.2, 132.4),
        AngleDefinition("O6", "C6", "C5", 128.4, 0.6, 782242, 128.5, 1.0, 123.2, 133.1, 118.0, 138.6),
        AngleDefinition("N2", "C2", "N1", 116.5, 0.6, 782269, 116.2, 1.1, 110.6, 120.9, 102.4, 125.2),
        AngleDefinition("N2", "C2", "N3", 119.9, 0.6, 782269, 119.8, 1.0, 114.8, 124.4, 108.9, 132.2),
    ],
    "U": [
        AngleDefinition("C6", "N1", "C2", 121.1, 0.5, 387309, 120.8, 1.0, 116.4, 125.1, 111.8, 128.5),
        AngleDefinition("N1", "C2", "N3", 114.9, 0.6, 387309, 115.2, 0.9, 111.7, 120.0, 108.0, 125.0),
        AngleDefinition("C2", "N3", "C4", 127.0, 0.5, 387309, 126.9, 0.8, 122.5, 130.1, 118.3, 134.4),
        AngleDefinition("N3", "C4", "C5", 114.5, 0.6, 387310, 114.6, 0.8, 110.2, 118.5, 103.9, 122.1),
        AngleDefinition("C4", "C5", "C6", 119.7, 0.6, 387310, 119.8, 0.8, 116.7, 124.3, 113.1, 128.7),
        AngleDefinition("C5", "C6", "N1", 122.7, 0.5, 387310, 122.7, 0.9, 118.1, 126.4, 113.3, 130.4),
        AngleDefinition("O2", "C2", "N1", 122.8, 0.7, 387309, 123.0, 1.3, 117.3, 128.5, 111.8, 132.3),
        AngleDefinition("O2", "C2", "N3", 122.3, 0.6, 387309, 121.8, 1.3, 115.4, 127.0, 110.7, 131.0),
        AngleDefinition("O4", "C4", "C5", 126.0, 0.7, 387305, 126.0, 1.0, 121.3, 131.2, 117.3, 136.3),
        AngleDefinition("O4", "C4", "N3", 119.5, 0.7, 387305, 119.5, 1.1, 114.1, 124.7, 109.3, 130.1),
    ],
    "T": [
        AngleDefinition("C6", "N1", "C2", 121.2, 0.5, 41088, 121.2, 0.3, 119.0, 123.0, 116.8, 124.7),
        AngleDefinition("N1", "C2", "N3", 114.7, 0.6, 41088, 114.7, 0.4, 113.0, 119.9, 111.7, 123.5),
        AngleDefinition("C2", "N3", "C4", 127.1, 0.5, 41088, 127.2, 0.4, 121.9, 129.0, 119.0, 130.7),
        AngleDefinition("N3", "C4", "C5", 115.2, 0.5, 41088, 115.3, 0.3, 113.7, 118.7, 111.1, 122.1),
        AngleDefinition("C4", "C5", "C6", 118.1, 0.5, 41088, 118.1, 0.3, 116.5, 120.3, 113.7, 124.0),
        AngleDefinition("C5", "C6", "N1", 123.6, 0.5, 41088, 123.5, 0.4, 120.4, 125.1, 117.3, 127.8),
        AngleDefinition("O2", "C2", "N1", 123.0, 0.7, 41086, 123.3, 0.5, 120.1, 126.7, 116.3, 132.0),
        AngleDefinition("O2", "C2", "N3", 122.3, 0.6, 41086, 122.0, 0.6, 117.5, 124.2, 113.2, 126.7),
        AngleDefinition("O4", "C4", "C5", 125.0, 0.7, 41086, 123.3, 1.3, 119.9, 127.3, 117.4, 130.0),
        AngleDefinition("O4", "C4", "N3", 119.8, 0.6, 41086, 121.4, 1.2, 117.2, 124.3, 115.0, 126.9),
        AngleDefinition("C7", "C5", "C4", 118.7, 0.6, 40843, 119.9, 0.7, 116.9, 123.1, 114.1, 125.4),
        AngleDefinition("C7", "C5", "C6", 123.2, 0.6, 40843, 122.1, 0.7, 118.5, 124.5, 114.9, 126.9),
    ],
    "C": [
        AngleDefinition("C6", "N1", "C2", 120.3, 0.5, 596004, 120.1, 0.9, 116.5, 124.3, 112.4, 128.4),
        AngleDefinition("N1", "C2", "N3", 119.1, 0.6, 596005, 119.3, 0.8, 116.0, 123.2, 112.4, 127.4),
        AngleDefinition("C2", "N3", "C4", 120.1, 0.5, 596005, 119.9, 0.7, 116.1, 122.8, 112.4, 127.8),
        AngleDefinition("N3", "C4", "C5", 121.6, 0.6, 596004, 121.9, 0.6, 118.9, 125.1, 114.3, 128.0),
        AngleDefinition("C4", "C5", "C6", 117.5, 0.5, 596003, 117.5, 0.6, 114.8, 120.8, 111.9, 123.9),
        AngleDefinition("C5", "C6", "N1", 121.2, 0.6, 596003, 121.2, 0.9, 117.0, 124.7, 113.5, 127.4),
        AngleDefinition("O2", "C2", "N1", 118.8, 0.8, 596004, 119.1, 1.2, 113.8, 124.2, 108.9, 131.2),
        AngleDefinition("O2", "C2", "N3", 122.0, 0.6, 596004, 121.6, 1.1, 116.3, 126.5, 109.3, 130.8),
        AngleDefinition("N4", "C4", "C5", 120.3, 0.7, 596000, 120.1, 0.9, 115.4, 124.4, 111.4, 129.4),
        AngleDefinition("N4", "C4", "N3", 118.1, 0.6, 596001, 118.0, 1.0, 113.4, 122.6, 109.2, 126.6),
    ],
}
BASES_ANGLES["DA"] = BASES_ANGLES["A"]
BASES_ANGLES["DG"] = BASES_ANGLES["G"]
BASES_ANGLES["DU"] = BASES_ANGLES["U"]
BASES_ANGLES["DT"] = BASES_ANGLES["T"]
BASES_ANGLES["DC"] = BASES_ANGLES["C"]


class Validator:
    """
    Base validator class
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        self.geometry = geometry
        self.csd_sig = csd_sig

    def validate(self):
        pass


class BasesValidator(Validator):
    """
    Validator for nucleotde basees
    """

    # pylint: disable=too-few-public-methods

    def _validate_bonds(self, res_name, resseq, chain):
        # pylint: disable=too-many-locals

        records = []
        for definition in BASES_BONDS[res_name]:
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
        for definition in BASES_ANGLES[res_name]:
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
