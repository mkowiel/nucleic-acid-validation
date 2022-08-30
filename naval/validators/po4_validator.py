from naval.validators.validator import Validator
from naval.restraint_definition import AngleDefinition
from naval.restraint_definition import BondDefinition
from naval.nucleotide_geometry import NucleotideGeometry


PO4_BONDS = {
    "PO4==AA_0": [
        BondDefinition("PO4==AA_0", "OP1", "P", 1.487, 0.01, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_0", "OP2", "P", 1.483, 0.01, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_0", "O3'", "P", 1.580, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_0", "O5'", "P", 1.603, 0.011, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_0", "O3'", "C3'", 1.422, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_0", "O5'", "C5'", 1.428, 0.013, None, None, None, None, None, None, None),
    ],
    "PO4==AA_1": [
        BondDefinition("PO4==AA_1", "OP1", "P", 1.483, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_1", "OP2", "P", 1.487, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_1", "O3'", "P", 1.603, 0.011, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_1", "O5'", "P", 1.580, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_1", "O3'", "C3'", 1.422, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_1", "O5'", "C5'", 1.428, 0.013, None, None, None, None, None, None, None),
    ],
    "PO4==AA_2": [
        BondDefinition("PO4==AA_2", "OP1", "P", 1.487, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_2", "OP2", "P", 1.483, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_2", "O3'", "P", 1.603, 0.011, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_2", "O5'", "P", 1.580, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_2", "O3'", "C3'", 1.422, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_2", "O5'", "C5'", 1.428, 0.013, None, None, None, None, None, None, None),
    ],
    "PO4==AA_3": [
        BondDefinition("PO4==AA_3", "OP1", "P", 1.483, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_3", "OP2", "P", 1.487, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_3", "O3'", "P", 1.580, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_3", "O5'", "P", 1.603, 0.011, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_3", "O3'", "C3'", 1.422, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AA_3", "O5'", "C5'", 1.428, 0.013, None, None, None, None, None, None, None),
    ],
    "PO4==AS_0": [
        BondDefinition("PO4==AS_0", "OP1", "P", 1.484, 0.012, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_0", "OP2", "P", 1.478, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_0", "O3'", "P", 1.599, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_0", "O5'", "P", 1.601, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_0", "O3'", "C3'", 1.438, 0.007, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_0", "O5'", "C5'", 1.437, 0.017, None, None, None, None, None, None, None),
    ],
    "PO4==AS_1": [
        BondDefinition("PO4==AS_1", "OP1", "P", 1.478, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_1", "OP2", "P", 1.484, 0.012, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_1", "O3'", "P", 1.601, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_1", "O5'", "P", 1.599, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_1", "O3'", "C3'", 1.438, 0.007, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_1", "O5'", "C5'", 1.437, 0.017, None, None, None, None, None, None, None),
    ],
    "PO4==AS_2": [
        BondDefinition("PO4==AS_2", "OP1", "P", 1.484, 0.012, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_2", "OP2", "P", 1.478, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_2", "O3'", "P", 1.601, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_2", "O5'", "P", 1.599, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_2", "O3'", "C3'", 1.438, 0.007, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_2", "O5'", "C5'", 1.437, 0.017, None, None, None, None, None, None, None),
    ],
    "PO4==AS_3": [
        BondDefinition("PO4==AS_3", "OP1", "P", 1.478, 0.010, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_3", "OP2", "P", 1.484, 0.012, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_3", "O3'", "P", 1.599, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_3", "O5'", "P", 1.601, 0.016, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_3", "O3'", "C3'", 1.438, 0.007, None, None, None, None, None, None, None),
        BondDefinition("PO4==AS_3", "O5'", "C5'", 1.437, 0.017, None, None, None, None, None, None, None),
    ],
}


PO4_ANGLES = {
    "PO4==AA_0": [
        AngleDefinition("PO4==AA_0", "OP1", "P", "OP2", 117.6, 1.2, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O3'", 106.2, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O5'", 110.2, 1.3, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O3'", 112.2, 1.0, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O5'", 109.3, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "O3'", "P", "O5'", 99.9, 0.7, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "P", "O3'", "C3'", 120.8, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_0", "P", "O5'", "C5'", 120.3, 1.9, None, None, None, None, None, None, None),
    ],
    "PO4==AA_1": [
        AngleDefinition("PO4==AA_1", "OP1", "P", "OP2", 117.6, 1.2, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O3'", 109.3, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O5'", 112.2, 1.0, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O3'", 110.2, 1.3, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O5'", 106.2, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "O3'", "P", "O5'", 99.9, 0.7, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "P", "O3'", "C3'", 120.3, 1.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_1", "P", "O5'", "C5'", 120.8, 1.1, None, None, None, None, None, None, None),
    ],
    "PO4==AA_2": [
        AngleDefinition("PO4==AA_2", "OP1", "P", "OP2", 117.6, 1.2, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O3'", 110.2, 1.3, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O5'", 106.2, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O3'", 109.3, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O5'", 112.2, 1.0, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "O3'", "P", "O5'", 99.9, 0.7, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "P", "O3'", "C3'", 120.3, 1.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_2", "P", "O5'", "C5'", 120.8, 1.1, None, None, None, None, None, None, None),
    ],
    "PO4==AA_3": [
        AngleDefinition("PO4==AA_3", "OP1", "P", "OP2", 117.6, 1.2, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O3'", 112.2, 1.0, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O5'", 109.3, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O3'", 106.2, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O5'", 110.2, 1.3, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "O3'", "P", "O5'", 99.9, 0.7, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "P", "O3'", "C3'", 120.8, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AA_3", "P", "O5'", "C5'", 120.3, 1.9, None, None, None, None, None, None, None),
    ],
    "PO4==AS_0": [
        AngleDefinition("PO4==AS_0", "OP1", "P", "OP2", 119.9, 1.6, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O3'", 104.5, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O5'", 110.3, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O3'", 111.5, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O5'", 105.2, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "O3'", "P", "O5'", 104.2, 1.5, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "P", "O3'", "C3'", 121.5, 3.0, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_0", "P", "O5'", "C5'", 121.6, 2.8, None, None, None, None, None, None, None),
    ],
    "PO4==AS_1": [
        AngleDefinition("PO4==AS_1", "OP1", "P", "OP2", 119.9, 1.6, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O3'", 105.2, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O5'", 111.5, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O3'", 110.3, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O5'", 104.5, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "O3'", "P", "O5'", 104.2, 1.5, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "P", "O3'", "C3'", 121.6, 2.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_1", "P", "O5'", "C5'", 121.5, 3.0, None, None, None, None, None, None, None),
    ],
    "PO4==AS_2": [
        AngleDefinition("PO4==AS_2", "OP1", "P", "OP2", 119.9, 1.6, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O3'", 110.3, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O5'", 104.5, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O3'", 105.2, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O5'", 111.5, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "O3'", "P", "O5'", 104.2, 1.5, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "P", "O3'", "C3'", 121.6, 2.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_2", "P", "O5'", "C5'", 121.5, 3.0, None, None, None, None, None, None, None),
    ],
    "PO4==AS_3": [
        AngleDefinition("PO4==AS_3", "OP1", "P", "OP2", 119.9, 1.6, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O3'", 111.5, 1.1, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O5'", 105.2, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O3'", 104.5, 0.9, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O5'", 110.3, 0.8, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "O3'", "P", "O5'", 104.2, 1.5, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "P", "O3'", "C3'", 121.5, 3.0, None, None, None, None, None, None, None),
        AngleDefinition("PO4==AS_3", "P", "O5'", "C5'", 121.6, 2.8, None, None, None, None, None, None, None),
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
