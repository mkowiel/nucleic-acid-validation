from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.sugar_basic_validator import BASIC_SUGAR_ANGLES, BASIC_SUGAR_BONDS
from naval.validators.validator import Validator

# pylint: disable=too-many-lines

SUGAR_PUCER_BASED_SUGAR_BONDS = {
    "pucker==A_G_C2p_endo": [
        BondDefinition("pucker==A_G_C2p_endo", "C1'", "C2'", 0, 0, 1.528, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C2'", "C3'", 0, 0, 1.528, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C3'", "C4'", 0, 0, 1.526, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C4'", "O4'", 0, 0, 1.453, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C1'", "O4'", 0, 0, 1.413, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C4'", "C5'", 0, 0, 1.510, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C2'", "O2'", 0, 0, 1.410, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C2p_endo", "C1'", "N9", 0, 0, 1.457, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==A_G_C3p_endo": [
        BondDefinition("pucker==A_G_C3p_endo", "C1'", "C2'", 0, 0, 1.533, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C2'", "C3'", 0, 0, 1.529, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C3'", "C4'", 0, 0, 1.522, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C4'", "O4'", 0, 0, 1.445, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C1'", "O4'", 0, 0, 1.413, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C4'", "C5'", 0, 0, 1.504, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C2'", "O2'", 0, 0, 1.416, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_C3p_endo", "C1'", "N9", 0, 0, 1.470, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==A_G_other": [
        BondDefinition("pucker==A_G_other", "C1'", "C2'", 0, 0, 1.541, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C2'", "C3'", 0, 0, 1.537, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C3'", "C4'", 0, 0, 1.529, 0.005, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C4'", "O4'", 0, 0, 1.443, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C1'", "O4'", 0, 0, 1.421, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C4'", "C5'", 0, 0, 1.506, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C2'", "O2'", 0, 0, 1.413, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==A_G_other", "C1'", "N9", 0, 0, 1.454, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==U_T_C_C2p_endo": [
        BondDefinition("pucker==U_T_C_C2p_endo", "C1'", "C2'", 0, 0, 1.531, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C2'", "C3'", 0, 0, 1.526, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C3'", "C4'", 0, 0, 1.527, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C4'", "O4'", 0, 0, 1.452, 0.006, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C1'", "O4'", 0, 0, 1.411, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C4'", "C5'", 0, 0, 1.508, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C2'", "O2'", 0, 0, 1.408, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C2p_endo", "C1'", "N1", 0, 0, 1.471, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==U_T_C_C3p_endo": [
        BondDefinition("pucker==U_T_C_C3p_endo", "C1'", "C2'", 0, 0, 1.533, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C2'", "C3'", 0, 0, 1.526, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C3'", "C4'", 0, 0, 1.520, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C4'", "O4'", 0, 0, 1.449, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C1'", "O4'", 0, 0, 1.411, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C4'", "C5'", 0, 0, 1.507, 0.007, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C2'", "O2'", 0, 0, 1.416, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==U_T_C_C3p_endo", "C1'", "N1", 0, 0, 1.488, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==DA_DG_C2p_endo": [
        BondDefinition("pucker==DA_DG_C2p_endo", "C1'", "C2'", 0, 0, 1.517, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C2p_endo", "C2'", "C3'", 0, 0, 1.518, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C2p_endo", "C3'", "C4'", 0, 0, 1.528, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C2p_endo", "C4'", "O4'", 0, 0, 1.446, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C2p_endo", "C1'", "O4'", 0, 0, 1.421, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C2p_endo", "C4'", "C5'", 0, 0, 1.511, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C2p_endo", "C1'", "N9", 0, 0, 1.456, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==DA_DG_C3p_endo": [
        BondDefinition("pucker==DA_DG_C3p_endo", "C1'", "C2'", 0, 0, 1.528, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C3p_endo", "C2'", "C3'", 0, 0, 1.516, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C3p_endo", "C3'", "C4'", 0, 0, 1.520, 0.012, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C3p_endo", "C4'", "O4'", 0, 0, 1.443, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C3p_endo", "C1'", "O4'", 0, 0, 1.420, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C3p_endo", "C4'", "C5'", 0, 0, 1.505, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_C3p_endo", "C1'", "N9", 0, 0, 1.467, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==DA_DG_other": [
        BondDefinition("pucker==DA_DG_other", "C1'", "C2'", 0, 0, 1.525, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_other", "C2'", "C3'", 0, 0, 1.529, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_other", "C3'", "C4'", 0, 0, 1.531, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_other", "C4'", "O4'", 0, 0, 1.438, 0.006, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_other", "C1'", "O4'", 0, 0, 1.427, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_other", "C4'", "C5'", 0, 0, 1.504, 0.006, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DA_DG_other", "C1'", "N9", 0, 0, 1.454, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==DU_DT_DC_C2p_endo": [
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "C2'", 0, 0, 1.519, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C2'", "C3'", 0, 0, 1.518, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C3'", "C4'", 0, 0, 1.526, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C4'", "O4'", 0, 0, 1.447, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "O4'", 0, 0, 1.419, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C4'", "C5'", 0, 0, 1.511, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "N1", 0, 0, 1.474, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==DU_DT_DC_C3p_endo": [
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "C2'", 0, 0, 1.518, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", 0, 0, 1.518, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", 0, 0, 1.518, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "O4'", 0, 0, 1.446, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "O4'", 0, 0, 1.413, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "C5'", 0, 0, 1.509, 0.011, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "N1", 0, 0, 1.492, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
    "pucker==DU_DT_DC_other": [
        BondDefinition("pucker==DU_DT_DC_other", "C1'", "C2'", 0, 0, 1.526, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_other", "C2'", "C3'", 0, 0, 1.516, 0.013, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_other", "C3'", "C4'", 0, 0, 1.533, 0.009, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_other", "C4'", "O4'", 0, 0, 1.438, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_other", "C1'", "O4'", 0, 0, 1.413, 0.016, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_other", "C4'", "C5'", 0, 0, 1.509, 0.010, 0, 1.5, 0.01, 1, 2, 1, 2),
        BondDefinition("pucker==DU_DT_DC_other", "C1'", "N1", 0, 0, 1.472, 0.008, 0, 1.5, 0.01, 1, 2, 1, 2),
    ],
}


SUGAR_PUCER_BASED_SUGAR_ANGLES = {
    "pucker==A_G_C2p_endo": [
        AngleDefinition("pucker==A_G_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.7, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 106.0, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.7, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 105.9, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 112.1, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 114.0, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 109.9, 2.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.5, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 115.6, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 108.8, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "N9", "C1'", "O4'", 0, 0, 0, 108.2, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "N9", "C1'", "C2'", 0, 0, 0, 114.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.9, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C1'", "N9", "C4", 0, 0, 0, 127.3, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C2p_endo", "C1'", "N9", "C8", 0, 0, 0, 126.7, 1.7, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==A_G_C3p_endo": [
        AngleDefinition("pucker==A_G_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.4, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.4, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.4, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.1, 0.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 107.5, 0.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 108.9, 2.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 110.1, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 114.0, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 112.5, 2.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 115.8, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.2, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "N9", "C1'", "O4'", 0, 0, 0, 108.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "N9", "C1'", "C2'", 0, 0, 0, 112.9, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 110.5, 2.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C1'", "N9", "C4", 0, 0, 0, 125.8, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_C3p_endo", "C1'", "N9", "C8", 0, 0, 0, 128.1, 1.6, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==A_G_other": [
        AngleDefinition("pucker==A_G_other", "C1'", "C2'", "C3'", 0, 0, 0, 103.3, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C2'", "C3'", "C4'", 0, 0, 0, 103.0, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C3'", "C4'", "O4'", 0, 0, 0, 104.7, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C1'", "O4'", "C4'", 0, 0, 0, 108.2, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C2'", "C1'", "O4'", 0, 0, 0, 107.0, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C1'", "C2'", "O2'", 0, 0, 0, 110.7, 3.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C3'", "C2'", "O2'", 0, 0, 0, 112.7, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C2'", "C3'", "O3'", 0, 0, 0, 114.1, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C4'", "C3'", "O3'", 0, 0, 0, 111.7, 2.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C3'", "C4'", "C5'", 0, 0, 0, 115.1, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C5'", "C4'", "O4'", 0, 0, 0, 108.7, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "N9", "C1'", "O4'", 0, 0, 0, 107.6, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "N9", "C1'", "C2'", 0, 0, 0, 115.2, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C4'", "C5'", "O5'", 0, 0, 0, 111.7, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C1'", "N9", "C4", 0, 0, 0, 126.3, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==A_G_other", "C1'", "N9", "C8", 0, 0, 0, 127.5, 2.1, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==U_T_C_C2p_endo": [
        AngleDefinition("pucker==U_T_C_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.7, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 106.0, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.8, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 105.8, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 112.4, 1.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 113.3, 2.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 108.7, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 110.6, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 115.4, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "N1", "C1'", "O4'", 0, 0, 0, 108.0, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "N1", "C1'", "C2'", 0, 0, 0, 113.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 110.0, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C1'", "N1", "C2", 0, 0, 0, 118.6, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C2p_endo", "C1'", "N1", "C6", 0, 0, 0, 120.1, 1.2, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==U_T_C_C3p_endo": [
        AngleDefinition("pucker==U_T_C_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.2, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 101.7, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.0, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.7, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 107.3, 0.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 108.7, 2.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 110.6, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 113.5, 2.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 112.6, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 116.4, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.8, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "N1", "C1'", "O4'", 0, 0, 0, 108.8, 0.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "N1", "C1'", "C2'", 0, 0, 0, 112.3, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.2, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C1'", "N1", "C2", 0, 0, 0, 116.6, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==U_T_C_C3p_endo", "C1'", "N1", "C6", 0, 0, 0, 122.4, 0.8, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==DA_DG_C2p_endo": [
        AngleDefinition("pucker==DA_DG_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.9, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 103.0, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 106.0, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.6, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 106.0, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 109.6, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.4, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 114.4, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 108.9, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "N9 ", "C1'", "O4'", 0, 0, 0, 108.0, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "N9 ", "C1'", "C2'", 0, 0, 0, 115.0, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.0, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C1'", "N9", "C4", 0, 0, 0, 126.8, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C2p_endo", "C1'", "N9", "C8", 0, 0, 0, 126.8, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==DA_DG_C3p_endo": [
        AngleDefinition("pucker==DA_DG_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 103.5, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 101.8, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.9, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.0, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 106.7, 0.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 113.5, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.1, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 116.3, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.6, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "N9 ", "C1'", "O4'", 0, 0, 0, 107.5, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "N9 ", "C1'", "C2'", 0, 0, 0, 113.6, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.9, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C1'", "N9", "C4", 0, 0, 0, 125.7, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_C3p_endo", "C1'", "N9", "C8", 0, 0, 0, 128.1, 1.3, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==DA_DG_other": [
        AngleDefinition("pucker==DA_DG_other", "C1'", "C2'", "C3'", 0, 0, 0, 103.8, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C3'", "C4'", "O4'", 0, 0, 0, 104.8, 1.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C1'", "O4'", "C4'", 0, 0, 0, 107.3, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C2'", "C1'", "O4'", 0, 0, 0, 105.7, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C2'", "C3'", "O3'", 0, 0, 0, 112.3, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C4'", "C3'", "O3'", 0, 0, 0, 110.2, 1.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C3'", "C4'", "C5'", 0, 0, 0, 114.9, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C5'", "C4'", "O4'", 0, 0, 0, 108.7, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "N9 ", "C1'", "O4'", 0, 0, 0, 107.7, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "N9 ", "C1'", "C2'", 0, 0, 0, 115.4, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C4'", "C5'", "O5'", 0, 0, 0, 110.9, 2.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C1'", "N9", "C4", 0, 0, 0, 127.6, 1.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DA_DG_other", "C1'", "N9", "C8", 0, 0, 0, 126.2, 1.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==DU_DT_DC_C2p_endo": [
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 102.4, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 105.9, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.0, 0.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 106.2, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 109.5, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.8, 2.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 114.7, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.7, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "N1", "C1'", "O4'", 0, 0, 0, 107.8, 0.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "N1", "C1'", "C2'", 0, 0, 0, 114.0, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 110.1, 1.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "N1", "C2", 0, 0, 0, 118.3, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "N1", "C6", 0, 0, 0, 120.3, 1.1, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==DU_DT_DC_C3p_endo": [
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 103.0, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.2, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.8, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.1, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 107.1, 0.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 112.6, 2.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 112.1, 2.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 116.1, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.5, 0.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "N1", "C1'", "O4'", 0, 0, 0, 108.1, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "N1", "C1'", "C2'", 0, 0, 0, 113.0, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.1, 2.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "N1", "C2", 0, 0, 0, 116.6, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "N1", "C6", 0, 0, 0, 122.0, 0.9, 0, 180, 180, 0, 360, 0, 360),
    ],
    "pucker==DU_DT_DC_other": [
        AngleDefinition("pucker==DU_DT_DC_other", "C1'", "C2'", "C3'", 0, 0, 0, 103.9, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C2'", "C3'", "C4'", 0, 0, 0, 103.5, 1.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C3'", "C4'", "O4'", 0, 0, 0, 106.1, 1.5, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C1'", "O4'", "C4'", 0, 0, 0, 109.1, 0.9, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C2'", "C1'", "O4'", 0, 0, 0, 106.3, 1.1, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C2'", "C3'", "O3'", 0, 0, 0, 111.6, 2.6, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C4'", "C3'", "O3'", 0, 0, 0, 111.5, 1.4, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C3'", "C4'", "C5'", 0, 0, 0, 115.1, 1.2, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C5'", "C4'", "O4'", 0, 0, 0, 108.7, 1.0, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "N1", "C1'", "O4'", 0, 0, 0, 107.7, 0.8, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "N1", "C1'", "C2'", 0, 0, 0, 114.5, 1.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C4'", "C5'", "O5'", 0, 0, 0, 109.7, 2.3, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C1'", "N1", "C2", 0, 0, 0, 119.0, 1.7, 0, 180, 180, 0, 360, 0, 360),
        AngleDefinition("pucker==DU_DT_DC_other", "C1'", "N1", "C6", 0, 0, 0, 119.4, 1.6, 0, 180, 180, 0, 360, 0, 360),
    ],
}


class SugarPuckerBasedSugarValidator(Validator):
    """
    Validator for nucleotde basees
    """

    # pylint: disable=too-few-public-methods
    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        super().__init__(geometry, csd_sig)

        self.basic_bonds_definition = BASIC_SUGAR_BONDS
        self.basic_angles_definition = BASIC_SUGAR_ANGLES

        self.bonds_definition = SUGAR_PUCER_BASED_SUGAR_BONDS
        self.angles_definition = SUGAR_PUCER_BASED_SUGAR_ANGLES

    def _atom_names_bonds(self, res_name: str) -> List[BondDefinition]:
        if res_name in ("A", "G"):
            return self.bonds_definition["pucker==A_G_C2p_endo"]
        if res_name in ("U", "T", "C"):
            return self.bonds_definition["pucker==U_T_C_C2p_endo"]
        if res_name in ("DA", "DG"):
            return self.bonds_definition["pucker==DA_DG_C2p_endo"]
        if res_name in ("DU", "DT", "DC"):
            return self.bonds_definition["pucker==DU_DT_DC_C2p_endo"]
        raise Exception("Non-standard residue")

    def _atom_names_angles(self, res_name: str) -> List[AngleDefinition]:
        if res_name in ("A", "G"):
            return self.angles_definition["pucker==A_G_C2p_endo"]
        if res_name in ("U", "T", "C"):
            return self.angles_definition["pucker==U_T_C_C2p_endo"]
        if res_name in ("DA", "DG"):
            return self.angles_definition["pucker==DA_DG_C2p_endo"]
        if res_name in ("DU", "DT", "DC"):
            return self.angles_definition["pucker==DU_DT_DC_C2p_endo"]
        raise Exception("Non-standard residue")

    def _find_bond_definitions(self, res_name: str, altloc: str, atom1_name: str, atom2_name: str) -> List[BondDefinition]:
        # pylint: disable=too-many-return-statements
        # pylint: disable=too-many-branches
        sugar_conformation = self.geometry.sugar_conformation.get(altloc, self.geometry.sugar_conformation.get("", None))

        if sugar_conformation == "C2'-endo":
            if res_name in ("A", "G"):
                return self.bonds_definition["pucker==A_G_C2p_endo"]
            if res_name in ("U", "T", "C"):
                return self.bonds_definition["pucker==U_T_C_C2p_endo"]
            if res_name in ("DA", "DG"):
                return self.bonds_definition["pucker==DA_DG_C2p_endo"]
            if res_name in ("DU", "DT", "DC"):
                return self.bonds_definition["pucker==DU_DT_DC_C2p_endo"]
            raise Exception("Non-standard residue")

        if sugar_conformation == "C3'-endo":
            if res_name in ("A", "G"):
                return self.bonds_definition["pucker==A_G_C3p_endo"]
            if res_name in ("U", "T", "C"):
                return self.bonds_definition["pucker==U_T_C_C3p_endo"]
            if res_name in ("DA", "DG"):
                return self.bonds_definition["pucker==DA_DG_C3p_endo"]
            if res_name in ("DU", "DT", "DC"):
                return self.bonds_definition["pucker==DU_DT_DC_C3p_endo"]
            raise Exception("Non-standard residue")

        if sugar_conformation == "other":
            if res_name in ("A", "G"):
                return self.bonds_definition["pucker==A_G_other"]
            if res_name in ("U", "T", "C"):
                # other is missing
                return self.basic_bonds_definition["sugar_basic==U_T_C"]
            if res_name in ("DA", "DG"):
                return self.bonds_definition["pucker==DA_DG_other"]
            if res_name in ("DU", "DT", "DC"):
                return self.bonds_definition["pucker==DU_DT_DC_other"]
            raise Exception("Non-standard residue")

        if res_name in ("A", "G"):
            return self.basic_bonds_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.basic_bonds_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.basic_bonds_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.basic_bonds_definition["sugar_basic==DU_DT_DC"]

        raise Exception("Non-standard residue")

    def _find_anlge_definitions(
        self, res_name: str, altloc: str, atom1_name: str, atom2_name: str, atom3_name: str
    ) -> List[AngleDefinition]:
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-return-statements
        # pylint: disable=too-many-branches

        sugar_conformation = self.geometry.sugar_conformation.get(altloc, self.geometry.sugar_conformation.get("", None))

        if sugar_conformation == "C2'-endo":
            if res_name in ("A", "G"):
                return self.angles_definition["pucker==A_G_C2p_endo"]
            if res_name in ("U", "T", "C"):
                return self.angles_definition["pucker==U_T_C_C2p_endo"]
            if res_name in ("DA", "DG"):
                return self.angles_definition["pucker==DA_DG_C2p_endo"]
            if res_name in ("DU", "DT", "DC"):
                return self.angles_definition["pucker==DU_DT_DC_C2p_endo"]
            raise Exception("Non-standard residue")

        if sugar_conformation == "C3'-endo":
            if res_name in ("A", "G"):
                return self.angles_definition["pucker==A_G_C3p_endo"]
            if res_name in ("U", "T", "C"):
                return self.angles_definition["pucker==U_T_C_C3p_endo"]
            if res_name in ("DA", "DG"):
                return self.angles_definition["pucker==DA_DG_C3p_endo"]
            if res_name in ("DU", "DT", "DC"):
                return self.angles_definition["pucker==DU_DT_DC_C3p_endo"]
            raise Exception("Non-standard residue")

        if sugar_conformation == "other":
            if res_name in ("A", "G"):
                return self.angles_definition["pucker==A_G_other"]
            if res_name in ("U", "T", "C"):
                # other is missing
                return self.basic_angles_definition["sugar_basic==U_T_C"]
            if res_name in ("DA", "DG"):
                return self.angles_definition["pucker==DA_DG_other"]
            if res_name in ("DU", "DT", "DC"):
                return self.angles_definition["pucker==DU_DT_DC_other"]
            raise Exception("Non-standard residue")

        if res_name in ("A", "G"):
            return self.basic_angles_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.basic_angles_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.basic_angles_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.basic_angles_definition["sugar_basic==DU_DT_DC"]

        raise Exception("Non-standard residue")
