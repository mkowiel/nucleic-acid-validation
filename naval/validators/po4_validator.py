from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.validator import Validator

# TODO: maybe  "O3'", "C3'" relative positions should be -1, -1

PO4_BONDS = {
    "PO4==AA_0": [
        BondDefinition("PO4==AA_0", "OP1", "P", 0, 0, 1.487, 0.01, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "OP2", "P", 0, 0, 1.483, 0.01, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O3'", "P", -1, 0, 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O5'", "P", 0, 0, 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O3'", "C3'", 0, 0, 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_0", "O5'", "C5'", 0, 0, 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AA_1": [
        BondDefinition("PO4==AA_1", "OP1", "P", 0, 0, 1.483, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "OP2", "P", 0, 0, 1.487, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O3'", "P", -1, 0, 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O5'", "P", 0, 0, 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O3'", "C3'", 0, 0, 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_1", "O5'", "C5'", 0, 0, 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AA_2": [
        BondDefinition("PO4==AA_2", "OP1", "P", 0, 0, 1.487, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "OP2", "P", 0, 0, 1.483, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O3'", "P", -1, 0, 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O5'", "P", 0, 0, 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O3'", "C3'", 0, 0, 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_2", "O5'", "C5'", 0, 0, 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AA_3": [
        BondDefinition("PO4==AA_3", "OP1", "P", 0, 0, 1.483, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "OP2", "P", 0, 0, 1.487, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O3'", "P", -1, 0, 1.580, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O5'", "P", 0, 0, 1.603, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O3'", "C3'", 0, 0, 1.422, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AA_3", "O5'", "C5'", 0, 0, 1.428, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_0": [
        BondDefinition("PO4==AS_0", "OP1", "P", 0, 0, 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "OP2", "P", 0, 0, 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O3'", "P", -1, 0, 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O5'", "P", 0, 0, 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O3'", "C3'", 0, 0, 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_0", "O5'", "C5'", 0, 0, 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_1": [
        BondDefinition("PO4==AS_1", "OP1", "P", 0, 0, 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "OP2", "P", 0, 0, 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O3'", "P", -1, 0, 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O5'", "P", 0, 0, 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O3'", "C3'", 0, 0, 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_1", "O5'", "C5'", 0, 0, 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_2": [
        BondDefinition("PO4==AS_2", "OP1", "P", 0, 0, 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "OP2", "P", 0, 0, 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O3'", "P", -1, 0, 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O5'", "P", 0, 0, 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O3'", "C3'", 0, 0, 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_2", "O5'", "C5'", 0, 0, 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "PO4==AS_3": [
        BondDefinition("PO4==AS_3", "OP1", "P", 0, 0, 1.478, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "OP2", "P", 0, 0, 1.484, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O3'", "P", -1, 0, 1.599, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O5'", "P", 0, 0, 1.601, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O3'", "C3'", 0, 0, 1.438, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("PO4==AS_3", "O5'", "C5'", 0, 0, 1.437, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "other==A_G": [
        BondDefinition("other==A_G", "OP1", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==A_G", "OP2", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==A_G", "O3'", "P", -1, 0, 1.607, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==A_G", "O5'", "P", 0, 0, 1.593, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==A_G", "O3'", "C3'", 0, 0, 1.420, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==A_G", "O5'", "C5'", 0, 0, 1.426, 0.014, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "other==DA_DG": [
        BondDefinition("other==DA_DG", "OP1", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DA_DG", "OP2", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DA_DG", "O3'", "P", -1, 0, 1.607, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DA_DG", "O5'", "P", 0, 0, 1.593, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DA_DG", "O3'", "C3'", 0, 0, 1.427, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DA_DG", "O5'", "C5'", 0, 0, 1.424, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "other==U_T_C": [
        BondDefinition("other==U_T_C", "OP1", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==U_T_C", "OP2", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==U_T_C", "O3'", "P", -1, 0, 1.607, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==U_T_C", "O5'", "P", 0, 0, 1.593, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==U_T_C", "O3'", "C3'", 0, 0, 1.420, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==U_T_C", "O5'", "C5'", 0, 0, 1.429, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "other==DU_DT_DC": [
        BondDefinition("other==DU_DT_DC", "OP1", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DU_DT_DC", "OP2", "P", 0, 0, 1.485, 0.017, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DU_DT_DC", "O3'", "P", -1, 0, 1.607, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DU_DT_DC", "O5'", "P", 0, 0, 1.593, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DU_DT_DC", "O3'", "C3'", 0, 0, 1.430, 0.015, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("other==DU_DT_DC", "O5'", "C5'", 0, 0, 1.425, 0.014, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
}


PO4_ANGLES = {
    "PO4==AA_0": [
        AngleDefinition("PO4==AA_0", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O3'", 0, 0, -1, 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O5'", 0, 0, 0, 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O3'", 0, 0, -1, 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O5'", 0, 0, 0, 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "P", "O3'", "C3'", 0, -1, -1, 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_0", "P", "O5'", "C5'", 0, 0, 0, 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AA_1": [
        AngleDefinition("PO4==AA_1", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O3'", 0, 0, -1, 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O5'", 0, 0, 0, 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O3'", 0, 0, -1, 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O5'", 0, 0, 0, 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "P", "O3'", "C3'", 0, -1, -1, 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_1", "P", "O5'", "C5'", 0, 0, 0, 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AA_2": [
        AngleDefinition("PO4==AA_2", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O3'", 0, 0, -1, 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O5'", 0, 0, 0, 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O3'", 0, 0, -1, 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O5'", 0, 0, 0, 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "P", "O3'", "C3'", 0, -1, -1, 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_2", "P", "O5'", "C5'", 0, 0, 0, 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AA_3": [
        AngleDefinition("PO4==AA_3", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O3'", 0, 0, -1, 112.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O5'", 0, 0, 0, 109.3, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O3'", 0, 0, -1, 106.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O5'", 0, 0, 0, 110.2, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "P", "O3'", "C3'", 0, -1, -1, 120.8, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AA_3", "P", "O5'", "C5'", 0, 0, 0, 120.3, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_0": [
        AngleDefinition("PO4==AS_0", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O3'", 0, 0, -1, 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O5'", 0, 0, 0, 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O3'", 0, 0, -1, 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O5'", 0, 0, 0, 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "P", "O3'", "C3'", 0, -1, -1, 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_0", "P", "O5'", "C5'", 0, 0, 0, 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_1": [
        AngleDefinition("PO4==AS_1", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O3'", 0, 0, -1, 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O5'", 0, 0, 0, 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O3'", 0, 0, -1, 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O5'", 0, 0, 0, 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "P", "O3'", "C3'", 0, -1, -1, 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_1", "P", "O5'", "C5'", 0, 0, 0, 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_2": [
        AngleDefinition("PO4==AS_2", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O3'", 0, 0, -1, 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O5'", 0, 0, 0, 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O3'", 0, 0, -1, 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O5'", 0, 0, 0, 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "P", "O3'", "C3'", 0, -1, -1, 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_2", "P", "O5'", "C5'", 0, 0, 0, 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
    ],
    "PO4==AS_3": [
        AngleDefinition("PO4==AS_3", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O3'", 0, 0, -1, 111.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O5'", 0, 0, 0, 105.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O3'", 0, 0, -1, 104.5, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O5'", 0, 0, 0, 110.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "P", "O3'", "C3'", 0, -1, -1, 121.5, 3.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("PO4==AS_3", "P", "O5'", "C5'", 0, 0, 0, 121.6, 2.8, 0, 180, 180, 0, 360, 0, 360),
    ],
    "other": [
        AngleDefinition("other", "OP1", "P", "OP2", 0, 0, 0, 119.6, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "OP1", "P", "O3'", 0, 0, -1, 107.7, 3.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "OP1", "P", "O5'", 0, 0, 0, 108.1, 2.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "OP2", "P", "O3'", 0, 0, -1, 108.3, 3.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "OP2", "P", "O5'", 0, 0, 0, 108.3, 2.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "O3'", "P", "O5'", -1, 0, 0, 104.0, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "P", "O3'", "C3'", 0, -1, -1, 119.7, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("other", "P", "O5'", "C5'", 0, 0, 0, 120.9, 1.6, 0, 180, 180, 0, 360, 0, 360),
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

    def _find_bond_definitions(
        self, res_name: str, altloc: str, atom1_name: str, atom2_name: str
    ) -> List[BondDefinition]:
        # pylint: disable=too-many-return-statements
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-branches
        # TODO: fix zeta next and zeta prev for C3'-O3' and C5'-O5' and for angles containing C3' O3' and C5' and O5'
        if "O3'" in atom1_name and "C3'" == atom2_name:
            alpha = None
            zeta = self.geometry.zeta_conformation.get(altloc, self.geometry.zeta_conformation.get("", None))
            if self.geometry.residue_entry.next_res is not None and self.geometry.residue_entry.next_res.geometry:
                next_geometry = self.geometry.residue_entry.next_res.geometry
                alpha = next_geometry.alpha_conformation.get(altloc, next_geometry.alpha_conformation.get("", None))
        else:
            alpha = self.geometry.alpha_conformation.get(altloc, self.geometry.alpha_conformation.get("", None))
            zeta = None
            if self.geometry.residue_entry.prev_res and self.geometry.residue_entry.prev_res.geometry:
                prev_geometry = self.geometry.residue_entry.prev_res.geometry
                zeta = prev_geometry.zeta_conformation.get(altloc, prev_geometry.zeta_conformation.get("", None))

        # print(altloc, atom1_name, atom2_name, zeta, alpha)
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
        if res_name in ("A", "G"):
            return self.bonds_definition["other==A_G"]
        if res_name in ("U", "T", "C"):
            return self.bonds_definition["other==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.bonds_definition["other==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.bonds_definition["other==DU_DT_DC"]
        raise Exception("Non-standard residue")

    def _find_anlge_definitions(
        self, res_name: str, altloc: str, atom1_name: str, atom2_name: str, atom3_name: str
    ) -> List[AngleDefinition]:
        # pylint: disable=too-many-return-statements
        # pylint: disable=too-many-arguments
        # TODO: fix zeta next and zeta prev for C3'-O3' and C5'-O5' and for angles containing C3' O3' and C5' and O5'
        alpha = self.geometry.alpha_conformation.get(altloc, self.geometry.alpha_conformation.get("", None))
        zeta = None
        if self.geometry.residue_entry.prev_res and self.geometry.residue_entry.prev_res.geometry:
            prev_geometry = self.geometry.residue_entry.prev_res.geometry
            zeta = prev_geometry.zeta_conformation.get(altloc, prev_geometry.zeta_conformation.get("", None))

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
        return self.angles_definition["other"]
