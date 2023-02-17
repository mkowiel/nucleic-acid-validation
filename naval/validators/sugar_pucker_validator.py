from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.sugar_basic_validator import BASIC_SUGAR_ANGLES, BASIC_SUGAR_BONDS
from naval.validators.validator import NonStandardResidueException, Validator

# pylint: disable=too-many-lines
SUGAR_PUCER_BASED_SUGAR_BONDS = {
    "pucker==A_G_C2p_endo": [
        BondDefinition("pucker==A_G_C2p_endo", "C1'", "C2'", 0, 0, 1.528, 0.010, 211232, 1.524, 0.006, 1.494, 1.545, 1.470, 1.569),
        BondDefinition("pucker==A_G_C2p_endo", "C2'", "C3'", 0, 0, 1.528, 0.009, 211232, 1.524, 0.006, 1.497, 1.549, 1.431, 1.567),
        BondDefinition("pucker==A_G_C2p_endo", "C3'", "C4'", 0, 0, 1.526, 0.008, 211232, 1.526, 0.007, 1.495, 1.551, 1.465, 1.572),
        BondDefinition("pucker==A_G_C2p_endo", "C4'", "O4'", 0, 0, 1.453, 0.008, 211232, 1.453, 0.006, 1.427, 1.473, 1.405, 1.543),
        BondDefinition("pucker==A_G_C2p_endo", "C1'", "O4'", 0, 0, 1.413, 0.008, 211232, 1.412, 0.006, 1.385, 1.431, 1.360, 1.515),
        BondDefinition("pucker==A_G_C2p_endo", "C4'", "C5'", 0, 0, 1.510, 0.010, 211228, 1.507, 0.008, 1.474, 1.535, 1.449, 1.563),
        BondDefinition("pucker==A_G_C2p_endo", "C2'", "O2'", 0, 0, 1.410, 0.009, 211143, 1.411, 0.004, 1.388, 1.430, 1.360, 1.469),
        BondDefinition("pucker==A_G_C2p_endo", "C1'", "N9", 0, 0, 1.457, 0.011, 211220, 1.463, 0.007, 1.431, 1.489, 1.385, 1.514),
    ],
    "pucker==A_G_C3p_endo": [
        BondDefinition("pucker==A_G_C3p_endo", "C1'", "C2'", 0, 0, 1.533, 0.012, 1212816, 1.527, 0.005, 1.502, 1.545, 1.482, 1.561),
        BondDefinition("pucker==A_G_C3p_endo", "C2'", "C3'", 0, 0, 1.529, 0.010, 1212816, 1.521, 0.006, 1.494, 1.540, 1.471, 1.555),
        BondDefinition("pucker==A_G_C3p_endo", "C3'", "C4'", 0, 0, 1.522, 0.008, 1212816, 1.518, 0.006, 1.489, 1.538, 1.472, 1.552),
        BondDefinition("pucker==A_G_C3p_endo", "C4'", "O4'", 0, 0, 1.445, 0.009, 1212816, 1.450, 0.005, 1.425, 1.466, 1.409, 1.481),
        BondDefinition("pucker==A_G_C3p_endo", "C1'", "O4'", 0, 0, 1.413, 0.009, 1212816, 1.412, 0.005, 1.389, 1.433, 1.374, 1.449),
        BondDefinition("pucker==A_G_C3p_endo", "C4'", "C5'", 0, 0, 1.504, 0.009, 1212808, 1.506, 0.007, 1.478, 1.531, 1.461, 1.559),
        BondDefinition("pucker==A_G_C3p_endo", "C2'", "O2'", 0, 0, 1.416, 0.009, 1212731, 1.418, 0.004, 1.396, 1.434, 1.377, 1.457),
        BondDefinition("pucker==A_G_C3p_endo", "C1'", "N9", 0, 0, 1.470, 0.011, 1212763, 1.472, 0.007, 1.442, 1.497, 1.403, 1.519),
    ],
    "pucker==A_G_other": [
        BondDefinition("pucker==A_G_other", "C1'", "C2'", 0, 0, 1.541, 0.010, 69590, 1.527, 0.008, 1.492, 1.553, 1.463, 1.577),
        BondDefinition("pucker==A_G_other", "C2'", "C3'", 0, 0, 1.537, 0.007, 69590, 1.527, 0.010, 1.489, 1.561, 1.279, 1.588),
        BondDefinition("pucker==A_G_other", "C3'", "C4'", 0, 0, 1.529, 0.005, 69590, 1.522, 0.009, 1.481, 1.554, 1.413, 1.594),
        BondDefinition("pucker==A_G_other", "C4'", "O4'", 0, 0, 1.443, 0.009, 69590, 1.447, 0.009, 1.415, 1.473, 1.373, 1.539),
        BondDefinition("pucker==A_G_other", "C1'", "O4'", 0, 0, 1.421, 0.010, 69590, 1.407, 0.009, 1.375, 1.435, 1.341, 1.522),
        BondDefinition("pucker==A_G_other", "C4'", "C5'", 0, 0, 1.506, 0.007, 69575, 1.508, 0.009, 1.471, 1.541, 1.445, 1.566),
        BondDefinition("pucker==A_G_other", "C2'", "O2'", 0, 0, 1.413, 0.008, 69555, 1.416, 0.006, 1.388, 1.438, 1.360, 1.508),
        BondDefinition("pucker==A_G_other", "C1'", "N9", 0, 0, 1.454, 0.007, 69565, 1.470, 0.010, 1.432, 1.510, 1.374, 1.580),
    ],
    "pucker==U_T_C_C2p_endo": [
        BondDefinition("pucker==U_T_C_C2p_endo", "C1'", "C2'", 0, 0, 1.531, 0.010, 128831, 1.524, 0.006, 1.497, 1.546, 1.475, 1.570),
        BondDefinition("pucker==U_T_C_C2p_endo", "C2'", "C3'", 0, 0, 1.526, 0.008, 128831, 1.525, 0.006, 1.499, 1.550, 1.487, 1.570),
        BondDefinition("pucker==U_T_C_C2p_endo", "C3'", "C4'", 0, 0, 1.527, 0.011, 128831, 1.527, 0.006, 1.499, 1.550, 1.481, 1.571),
        BondDefinition("pucker==U_T_C_C2p_endo", "C4'", "O4'", 0, 0, 1.452, 0.006, 128831, 1.453, 0.005, 1.430, 1.473, 1.413, 1.490),
        BondDefinition("pucker==U_T_C_C2p_endo", "C1'", "O4'", 0, 0, 1.411, 0.008, 128831, 1.413, 0.005, 1.386, 1.431, 1.373, 1.456),
        BondDefinition("pucker==U_T_C_C2p_endo", "C4'", "C5'", 0, 0, 1.508, 0.008, 128829, 1.508, 0.007, 1.476, 1.535, 1.454, 1.551),
        BondDefinition("pucker==U_T_C_C2p_endo", "C2'", "O2'", 0, 0, 1.408, 0.008, 128821, 1.412, 0.004, 1.391, 1.429, 1.373, 1.452),
        BondDefinition("pucker==U_T_C_C2p_endo", "C1'", "N1", 0, 0, 1.471, 0.010, 128823, 1.471, 0.010, 1.436, 1.535, 1.416, 1.564),
    ],
    "pucker==U_T_C_C3p_endo": [
        BondDefinition("pucker==U_T_C_C3p_endo", "C1'", "C2'", 0, 0, 1.533, 0.009, 967220, 1.528, 0.005, 1.502, 1.544, 1.484, 1.560),
        BondDefinition("pucker==U_T_C_C3p_endo", "C2'", "C3'", 0, 0, 1.526, 0.009, 967220, 1.521, 0.006, 1.494, 1.540, 1.471, 1.558),
        BondDefinition("pucker==U_T_C_C3p_endo", "C3'", "C4'", 0, 0, 1.520, 0.009, 967220, 1.518, 0.006, 1.489, 1.538, 1.469, 1.557),
        BondDefinition("pucker==U_T_C_C3p_endo", "C4'", "O4'", 0, 0, 1.449, 0.009, 967220, 1.450, 0.005, 1.424, 1.466, 1.406, 1.479),
        BondDefinition("pucker==U_T_C_C3p_endo", "C1'", "O4'", 0, 0, 1.411, 0.007, 967220, 1.412, 0.005, 1.390, 1.434, 1.375, 1.449),
        BondDefinition("pucker==U_T_C_C3p_endo", "C4'", "C5'", 0, 0, 1.507, 0.007, 967203, 1.506, 0.006, 1.478, 1.530, 1.459, 1.549),
        BondDefinition("pucker==U_T_C_C3p_endo", "C2'", "O2'", 0, 0, 1.416, 0.008, 967196, 1.418, 0.004, 1.396, 1.434, 1.378, 1.455),
        BondDefinition("pucker==U_T_C_C3p_endo", "C1'", "N1", 0, 0, 1.488, 0.010, 967105, 1.480, 0.011, 1.440, 1.547, 1.417, 1.584),
    ],
    "pucker==DA_DG_C2p_endo": [
        BondDefinition("pucker==DA_DG_C2p_endo", "C1'", "C2'", 0, 0, 1.517, 0.008, 66089, 1.520, 0.006, 1.493, 1.550, 1.464, 1.591),
        BondDefinition("pucker==DA_DG_C2p_endo", "C2'", "C3'", 0, 0, 1.518, 0.009, 66089, 1.525, 0.007, 1.496, 1.553, 1.454, 1.606),
        BondDefinition("pucker==DA_DG_C2p_endo", "C3'", "C4'", 0, 0, 1.528, 0.008, 66089, 1.527, 0.006, 1.497, 1.553, 1.471, 1.575),
        BondDefinition("pucker==DA_DG_C2p_endo", "C4'", "O4'", 0, 0, 1.446, 0.009, 66089, 1.451, 0.005, 1.426, 1.477, 1.395, 1.528),
        BondDefinition("pucker==DA_DG_C2p_endo", "C1'", "O4'", 0, 0, 1.421, 0.011, 66089, 1.415, 0.008, 1.387, 1.458, 1.352, 1.490),
        BondDefinition("pucker==DA_DG_C2p_endo", "C4'", "C5'", 0, 0, 1.511, 0.008, 66057, 1.515, 0.007, 1.487, 1.543, 1.442, 1.574),
        BondDefinition("pucker==DA_DG_C2p_endo", "C1'", "N9", 0, 0, 1.456, 0.008, 66077, 1.461, 0.007, 1.428, 1.491, 1.385, 1.523),
    ],
    "pucker==DA_DG_C3p_endo": [
        BondDefinition("pucker==DA_DG_C3p_endo", "C1'", "C2'", 0, 0, 1.528, 0.011, 7861, 1.526, 0.009, 1.486, 1.561, 1.457, 1.577),
        BondDefinition("pucker==DA_DG_C3p_endo", "C2'", "C3'", 0, 0, 1.516, 0.010, 7861, 1.522, 0.007, 1.488, 1.559, 1.481, 1.584),
        BondDefinition("pucker==DA_DG_C3p_endo", "C3'", "C4'", 0, 0, 1.520, 0.012, 7861, 1.518, 0.008, 1.482, 1.553, 1.472, 1.571),
        BondDefinition("pucker==DA_DG_C3p_endo", "C4'", "O4'", 0, 0, 1.443, 0.008, 7861, 1.445, 0.007, 1.409, 1.483, 1.398, 1.510),
        BondDefinition("pucker==DA_DG_C3p_endo", "C1'", "O4'", 0, 0, 1.420, 0.013, 7861, 1.419, 0.008, 1.390, 1.457, 1.380, 1.484),
        BondDefinition("pucker==DA_DG_C3p_endo", "C4'", "C5'", 0, 0, 1.505, 0.008, 7858, 1.514, 0.008, 1.480, 1.549, 1.430, 1.565),
        BondDefinition("pucker==DA_DG_C3p_endo", "C1'", "N9", 0, 0, 1.467, 0.013, 7858, 1.461, 0.008, 1.425, 1.496, 1.397, 1.513),
    ],
    "pucker==DA_DG_other": [
        BondDefinition("pucker==DA_DG_other", "C1'", "C2'", 0, 0, 1.525, 0.011, 27049, 1.523, 0.008, 1.490, 1.553, 1.452, 1.583),
        BondDefinition("pucker==DA_DG_other", "C2'", "C3'", 0, 0, 1.529, 0.011, 27049, 1.529, 0.009, 1.497, 1.561, 1.433, 1.597),
        BondDefinition("pucker==DA_DG_other", "C3'", "C4'", 0, 0, 1.531, 0.010, 27049, 1.526, 0.009, 1.489, 1.557, 1.453, 1.579),
        BondDefinition("pucker==DA_DG_other", "C4'", "O4'", 0, 0, 1.438, 0.006, 27049, 1.445, 0.007, 1.412, 1.476, 1.380, 1.513),
        BondDefinition("pucker==DA_DG_other", "C1'", "O4'", 0, 0, 1.427, 0.011, 27049, 1.410, 0.009, 1.378, 1.448, 1.352, 1.481),
        BondDefinition("pucker==DA_DG_other", "C4'", "C5'", 0, 0, 1.504, 0.006, 27000, 1.514, 0.008, 1.478, 1.550, 1.442, 1.590),
        BondDefinition("pucker==DA_DG_other", "C1'", "N9", 0, 0, 1.454, 0.010, 27034, 1.461, 0.007, 1.422, 1.494, 1.391, 1.575),
    ],
    "pucker==DU_DT_DC_C2p_endo": [
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "C2'", 0, 0, 1.519, 0.010, 49389, 1.520, 0.006, 1.495, 1.550, 1.468, 1.625),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C2'", "C3'", 0, 0, 1.518, 0.010, 49389, 1.524, 0.008, 1.492, 1.553, 1.224, 1.604),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C3'", "C4'", 0, 0, 1.526, 0.010, 49389, 1.528, 0.006, 1.500, 1.557, 1.467, 1.620),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C4'", "O4'", 0, 0, 1.447, 0.010, 49389, 1.451, 0.006, 1.424, 1.481, 1.332, 1.559),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "O4'", 0, 0, 1.419, 0.010, 49389, 1.414, 0.008, 1.386, 1.452, 1.284, 1.499),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C4'", "C5'", 0, 0, 1.511, 0.011, 49335, 1.514, 0.007, 1.485, 1.545, 1.416, 1.588),
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "N1", 0, 0, 1.474, 0.013, 49340, 1.488, 0.016, 1.430, 1.554, 1.351, 1.621),
    ],
    "pucker==DU_DT_DC_C3p_endo": [
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "C2'", 0, 0, 1.518, 0.010, 8671, 1.527, 0.009, 1.478, 1.560, 1.358, 1.594),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", 0, 0, 1.518, 0.010, 8671, 1.523, 0.008, 1.491, 1.562, 1.283, 1.607),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", 0, 0, 1.518, 0.008, 8671, 1.518, 0.008, 1.483, 1.553, 1.467, 1.649),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "O4'", 0, 0, 1.446, 0.009, 8671, 1.446, 0.008, 1.411, 1.478, 1.288, 1.511),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "O4'", 0, 0, 1.413, 0.009, 8671, 1.420, 0.008, 1.387, 1.463, 1.372, 1.613),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "C5'", 0, 0, 1.509, 0.011, 8662, 1.515, 0.007, 1.486, 1.549, 1.477, 1.576),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "N1", 0, 0, 1.492, 0.010, 8654, 1.491, 0.017, 1.438, 1.560, 1.367, 1.607),
    ],
    "pucker==DU_DT_DC_other": [
        BondDefinition("pucker==DU_DT_DC_other", "C1'", "C2'", 0, 0, 1.526, 0.013, 42575, 1.523, 0.008, 1.492, 1.556, 1.462, 1.621),
        BondDefinition("pucker==DU_DT_DC_other", "C2'", "C3'", 0, 0, 1.516, 0.013, 42575, 1.531, 0.010, 1.498, 1.561, 1.236, 1.581),
        BondDefinition("pucker==DU_DT_DC_other", "C3'", "C4'", 0, 0, 1.533, 0.009, 42575, 1.527, 0.008, 1.491, 1.558, 1.466, 1.592),
        BondDefinition("pucker==DU_DT_DC_other", "C4'", "O4'", 0, 0, 1.438, 0.008, 42575, 1.444, 0.007, 1.412, 1.470, 1.372, 1.499),
        BondDefinition("pucker==DU_DT_DC_other", "C1'", "O4'", 0, 0, 1.413, 0.016, 42575, 1.409, 0.009, 1.378, 1.451, 1.296, 1.493),
        BondDefinition("pucker==DU_DT_DC_other", "C4'", "C5'", 0, 0, 1.509, 0.010, 42539, 1.514, 0.007, 1.481, 1.544, 1.418, 1.583),
        BondDefinition("pucker==DU_DT_DC_other", "C1'", "N1", 0, 0, 1.472, 0.008, 42520, 1.489, 0.017, 1.426, 1.561, 1.344, 1.617),
    ],
}


SUGAR_PUCER_BASED_SUGAR_ANGLES = {
    "pucker==A_G_C2p_endo": [
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.2, 1.0, 211232, 101.4, 1.0, 97.7, 105.3, 93.5, 106.7
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.7, 0.7, 211232, 102.7, 0.8, 99.9, 106.1, 96.9, 108.8
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 106.0, 0.7, 211232, 106.2, 0.5, 104.2, 107.8, 102.2, 109.7
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.7, 0.8, 211232, 109.6, 0.7, 106.8, 111.9, 101.4, 112.9
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 105.9, 0.8, 211232, 105.8, 0.7, 102.7, 108.8, 100.1, 111.9
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 112.1, 2.3, 211143, 111.6, 0.9, 106.5, 116.3, 99.7, 122.6
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 114.0, 1.9, 211143, 114.4, 1.0, 109.3, 118.5, 103.0, 123.3
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 109.9, 2.4, 211224, 109.5, 1.7, 103.7, 120.1, 90.5, 130.5
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.5, 1.9, 211224, 109.4, 1.6, 102.8, 118.1, 96.1, 124.1
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 115.6, 1.3, 211228, 115.2, 1.1, 109.6, 119.6, 103.7, 123.0
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 108.8, 1.2, 211228, 109.1, 1.0, 104.9, 113.2, 101.5, 116.6
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "N9", "C1'", "O4'", 0, 0, 0, 108.2, 1.0, 211220, 108.5, 1.5, 103.4, 115.1, 100.2, 119.2
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "N9", "C1'", "C2'", 0, 0, 0, 114.2, 1.1, 211220, 114.2, 1.5, 107.8, 121.1, 102.9, 125.8
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.9, 1.8, 211222, 111.5, 1.1, 106.8, 115.8, 103.7, 118.7
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C1'", "N9", "C4", 0, 0, 0, 127.3, 1.6, 211209, 126.7, 2.4, 118.7, 136.3, 114.8, 144.7
        ),
        AngleDefinition(
            "pucker==A_G_C2p_endo", "C1'", "N9", "C8", 0, 0, 0, 126.7, 1.7, 211210, 127.1, 2.2, 118.7, 134.7, 111.8, 140.9
        ),
    ],
    "pucker==A_G_C3p_endo": [
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.4, 1.0, 1212816, 101.3, 0.6, 98.9, 103.9, 96.5, 106.2
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.4, 0.9, 1212816, 102.5, 0.7, 98.9, 105.1, 95.9, 106.7
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.4, 0.8, 1212816, 104.0, 0.6, 101.3, 106.5, 98.3, 108.1
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.1, 0.6, 1212816, 109.8, 0.4, 108.1, 111.3, 106.2, 112.5
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 107.5, 0.6, 1212816, 107.5, 0.5, 105.3, 109.1, 103.5, 110.8
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 108.9, 2.4, 1212731, 108.5, 1.0, 104.5, 114.7, 97.8, 119.8
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 110.1, 1.8, 1212731, 110.9, 1.2, 106.4, 117.6, 101.7, 123.2
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 114.0, 1.5, 1212798, 113.5, 1.1, 107.7, 118.9, 101.0, 125.4
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 112.5, 2.0, 1212798, 112.7, 1.3, 106.4, 117.5, 100.5, 123.3
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 115.8, 1.4, 1212808, 115.9, 1.0, 110.8, 119.9, 106.0, 124.0
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.2, 1.1, 1212808, 109.9, 0.8, 106.5, 113.1, 103.5, 116.3
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "N9", "C1'", "O4'", 0, 0, 0, 108.3, 0.8, 1212763, 108.5, 0.9, 105.0, 112.3, 102.4, 115.1
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "N9", "C1'", "C2'", 0, 0, 0, 112.9, 1.4, 1212763, 112.1, 1.4, 107.2, 118.0, 103.5, 122.4
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 110.5, 2.0, 1212789, 111.4, 1.0, 107.2, 115.1, 104.4, 117.7
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C1'", "N9", "C4", 0, 0, 0, 125.8, 1.8, 1212760, 126.7, 1.9, 119.9, 134.2, 114.2, 139.6
        ),
        AngleDefinition(
            "pucker==A_G_C3p_endo", "C1'", "N9", "C8", 0, 0, 0, 128.1, 1.6, 1212738, 127.2, 1.8, 119.8, 133.4, 114.6, 140.5
        ),
    ],
    "pucker==A_G_other": [
        AngleDefinition("pucker==A_G_other", "C1'", "C2'", "C3'", 0, 0, 0, 103.3, 1.3, 69590, 102.1, 1.5, 97.3, 106.5, 93.5, 108.8),
        AngleDefinition("pucker==A_G_other", "C2'", "C3'", "C4'", 0, 0, 0, 103.0, 0.8, 69590, 103.5, 1.3, 98.8, 107.5, 91.8, 109.7),
        AngleDefinition("pucker==A_G_other", "C3'", "C4'", "O4'", 0, 0, 0, 104.7, 1.1, 69590, 105.1, 1.3, 99.1, 108.5, 94.2, 110.6),
        AngleDefinition("pucker==A_G_other", "C1'", "O4'", "C4'", 0, 0, 0, 108.2, 1.4, 69590, 109.7, 1.1, 103.0, 112.4, 86.4, 114.0),
        AngleDefinition("pucker==A_G_other", "C2'", "C1'", "O4'", 0, 0, 0, 107.0, 1.5, 69590, 107.3, 1.3, 101.5, 110.8, 95.7, 113.1),
        AngleDefinition("pucker==A_G_other", "C1'", "C2'", "O2'", 0, 0, 0, 110.7, 3.1, 69555, 109.3, 1.7, 101.4, 117.3, 90.2, 123.6),
        AngleDefinition("pucker==A_G_other", "C3'", "C2'", "O2'", 0, 0, 0, 112.7, 1.2, 69555, 111.9, 2.0, 105.1, 120.8, 97.8, 125.6),
        AngleDefinition("pucker==A_G_other", "C2'", "C3'", "O3'", 0, 0, 0, 114.1, 1.5, 69589, 112.6, 2.6, 101.3, 125.3, 95.2, 133.4),
        AngleDefinition("pucker==A_G_other", "C4'", "C3'", "O3'", 0, 0, 0, 111.7, 2.8, 69589, 111.4, 2.3, 101.0, 120.7, 94.0, 128.9),
        AngleDefinition("pucker==A_G_other", "C3'", "C4'", "C5'", 0, 0, 0, 115.1, 1.6, 69575, 115.4, 1.7, 108.1, 122.8, 100.2, 126.0),
        AngleDefinition("pucker==A_G_other", "C5'", "C4'", "O4'", 0, 0, 0, 108.7, 0.9, 69575, 109.3, 1.5, 103.3, 115.2, 99.7, 129.7),
        AngleDefinition("pucker==A_G_other", "N9", "C1'", "O4'", 0, 0, 0, 107.6, 0.8, 69565, 108.4, 2.0, 102.5, 116.5, 98.6, 130.0),
        AngleDefinition("pucker==A_G_other", "N9", "C1'", "C2'", 0, 0, 0, 115.2, 1.5, 69565, 113.5, 2.5, 105.2, 124.4, 99.0, 132.6),
        AngleDefinition("pucker==A_G_other", "C4'", "C5'", "O5'", 0, 0, 0, 111.7, 1.8, 69572, 111.2, 1.5, 105.2, 117.0, 102.0, 121.9),
        AngleDefinition("pucker==A_G_other", "C1'", "N9", "C4", 0, 0, 0, 126.3, 1.6, 69554, 127.2, 3.0, 115.8, 138.2, 108.1, 148.7),
        AngleDefinition("pucker==A_G_other", "C1'", "N9", "C8", 0, 0, 0, 127.5, 2.1, 69552, 126.9, 2.8, 116.9, 137.5, 108.3, 147.4),
    ],
    "pucker==U_T_C_C2p_endo": [
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.3, 0.8, 128831, 101.4, 0.9, 98.2, 105.2, 95.8, 106.1
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.7, 0.9, 128831, 102.7, 0.7, 100.1, 105.9, 98.0, 107.8
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 106.0, 1.0, 128831, 106.2, 0.4, 104.4, 107.8, 102.7, 108.8
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.8, 0.8, 128831, 109.6, 0.6, 107.2, 111.8, 104.9, 112.6
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 105.8, 0.7, 128831, 105.7, 0.7, 103.1, 108.7, 100.7, 110.5
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 112.4, 1.7, 128821, 111.7, 0.9, 106.3, 115.7, 98.8, 120.6
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 113.3, 2.9, 128821, 114.5, 1.0, 109.7, 119.1, 103.6, 123.6
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 108.7, 2.3, 128825, 109.5, 1.6, 104.0, 119.4, 99.2, 132.4
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 110.6, 2.3, 128825, 109.5, 1.5, 103.6, 117.9, 98.0, 123.3
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 115.4, 1.2, 128829, 115.3, 1.0, 110.3, 119.6, 105.1, 124.0
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.5, 1.1, 128829, 109.2, 0.9, 105.1, 112.8, 102.3, 117.1
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "N1", "C1'", "O4'", 0, 0, 0, 108.0, 0.7, 128823, 108.6, 1.3, 103.8, 114.2, 98.4, 119.5
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "N1", "C1'", "C2'", 0, 0, 0, 113.5, 1.1, 128823, 114.2, 1.5, 103.7, 120.0, 98.2, 124.9
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 110.0, 1.9, 128826, 111.6, 1.1, 107.0, 115.6, 104.1, 118.1
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C1'", "N1", "C2", 0, 0, 0, 118.6, 1.3, 128813, 118.4, 2.3, 110.4, 128.0, 98.8, 134.0
        ),
        AngleDefinition(
            "pucker==U_T_C_C2p_endo", "C1'", "N1", "C6", 0, 0, 0, 120.1, 1.2, 128813, 120.8, 2.0, 112.5, 128.1, 104.6, 136.9
        ),
    ],
    "pucker==U_T_C_C3p_endo": [
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.2, 0.9, 967220, 101.4, 0.6, 99.1, 103.7, 96.7, 106.3
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 101.7, 0.9, 967220, 102.4, 0.7, 99.1, 104.8, 96.0, 106.7
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.0, 0.7, 967220, 103.9, 0.6, 101.3, 106.2, 98.4, 107.9
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.7, 0.7, 967220, 109.8, 0.4, 108.1, 111.1, 106.4, 112.5
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 107.3, 0.6, 967220, 107.6, 0.5, 105.4, 109.0, 103.6, 110.3
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C1'", "C2'", "O2'", 0, 0, 0, 108.7, 2.4, 967196, 108.5, 0.9, 104.4, 114.3, 98.8, 119.4
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C3'", "C2'", "O2'", 0, 0, 0, 110.6, 2.3, 967196, 110.9, 1.2, 106.5, 117.8, 102.1, 122.7
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 113.5, 2.1, 967207, 113.5, 1.1, 107.8, 118.7, 102.1, 124.2
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 112.6, 1.9, 967207, 112.8, 1.3, 106.7, 117.5, 101.5, 122.6
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 116.4, 1.4, 967203, 115.9, 1.0, 110.7, 119.6, 105.0, 123.2
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.8, 0.8, 967203, 109.9, 0.8, 106.6, 112.9, 104.0, 115.9
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "N1", "C1'", "O4'", 0, 0, 0, 108.8, 0.6, 967105, 108.8, 0.8, 105.4, 112.8, 102.1, 116.3
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "N1", "C1'", "C2'", 0, 0, 0, 112.3, 1.1, 967105, 112.3, 1.3, 107.8, 118.1, 97.9, 122.6
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.2, 1.9, 967128, 111.4, 1.0, 107.4, 115.0, 104.6, 117.8
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C1'", "N1", "C2", 0, 0, 0, 116.6, 1.1, 967097, 118.5, 1.9, 111.8, 126.6, 107.0, 131.6
        ),
        AngleDefinition(
            "pucker==U_T_C_C3p_endo", "C1'", "N1", "C6", 0, 0, 0, 122.4, 0.8, 967094, 121.2, 1.6, 114.3, 127.2, 110.0, 132.1
        ),
    ],
    "pucker==DA_DG_C2p_endo": [
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 101.9, 1.1, 66089, 101.4, 1.4, 95.6, 105.7, 91.1, 108.2
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 103.0, 0.9, 66089, 103.1, 0.9, 98.8, 106.3, 95.2, 109.6
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 106.0, 0.9, 66089, 105.9, 0.7, 102.6, 108.4, 100.4, 113.3
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 109.6, 0.9, 66089, 109.7, 1.0, 105.2, 112.8, 101.5, 114.5
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 106.0, 0.9, 66089, 105.8, 1.0, 101.9, 109.0, 98.1, 110.9
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 109.6, 2.3, 66035, 111.8, 2.0, 101.6, 120.2, 96.5, 126.4
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.4, 1.9, 66035, 109.0, 1.9, 100.7, 117.4, 94.9, 126.1
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 114.4, 1.1, 66057, 114.9, 1.3, 108.6, 120.0, 99.9, 124.9
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 108.9, 1.1, 66057, 109.8, 1.0, 104.7, 115.3, 100.2, 121.2
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "N9 ", "C1'", "O4'", 0, 0, 0, 108.0, 1.0, 66077, 108.2, 1.4, 101.0, 114.9, 95.4, 122.1
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "N9 ", "C1'", "C2'", 0, 0, 0, 115.0, 1.0, 66077, 114.6, 1.6, 108.0, 121.0, 98.1, 128.0
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.0, 2.3, 65959, 110.6, 1.2, 104.4, 115.8, 98.5, 121.0
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C1'", "N9", "C4", 0, 0, 0, 126.8, 1.9, 66066, 126.5, 0.9, 122.1, 130.3, 118.5, 133.8
        ),
        AngleDefinition(
            "pucker==DA_DG_C2p_endo", "C1'", "N9", "C8", 0, 0, 0, 126.8, 1.9, 66066, 127.5, 0.9, 123.9, 131.7, 119.6, 135.7
        ),
    ],
    "pucker==DA_DG_C3p_endo": [
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 103.5, 1.1, 7861, 102.4, 1.4, 97.6, 106.7, 96.6, 108.2
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 101.8, 1.2, 7861, 102.0, 1.4, 97.2, 106.4, 92.4, 108.1
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.9, 1.0, 7861, 104.6, 1.1, 98.9, 108.3, 96.4, 110.0
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.0, 1.0, 7861, 109.6, 0.9, 106.1, 112.5, 105.0, 113.8
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 106.7, 0.5, 7861, 107.3, 0.8, 103.8, 110.5, 99.2, 112.5
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 113.5, 1.6, 7841, 111.7, 2.1, 100.7, 122.1, 93.2, 127.8
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.1, 0.7, 7841, 110.1, 1.8, 102.2, 120.3, 98.7, 125.5
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 116.3, 0.8, 7858, 115.0, 1.3, 107.5, 121.3, 104.7, 132.7
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.6, 1.0, 7858, 109.6, 1.2, 104.1, 116.1, 101.8, 120.5
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "N9 ", "C1'", "O4'", 0, 0, 0, 107.5, 0.7, 7858, 108.5, 1.1, 102.8, 113.4, 95.7, 118.7
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "N9 ", "C1'", "C2'", 0, 0, 0, 113.6, 0.8, 7858, 113.8, 1.6, 106.6, 120.5, 103.8, 122.2
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.9, 1.3, 7849, 110.5, 1.4, 104.1, 116.3, 101.2, 124.3
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C1'", "N9", "C4", 0, 0, 0, 125.7, 1.0, 7856, 126.5, 1.1, 121.0, 131.6, 119.9, 140.4
        ),
        AngleDefinition(
            "pucker==DA_DG_C3p_endo", "C1'", "N9", "C8", 0, 0, 0, 128.1, 1.3, 7856, 127.5, 1.0, 123.0, 132.7, 115.5, 133.7
        ),
    ],
    "pucker==DA_DG_other": [
        AngleDefinition("pucker==DA_DG_other", "C1'", "C2'", "C3'", 0, 0, 0, 103.8, 1.3, 27049, 103.0, 1.5, 97.0, 107.1, 91.3, 110.5),
        AngleDefinition("pucker==DA_DG_other", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 1.3, 27049, 103.9, 1.3, 98.0, 107.5, 93.9, 110.6),
        AngleDefinition("pucker==DA_DG_other", "C3'", "C4'", "O4'", 0, 0, 0, 104.8, 1.6, 27049, 105.6, 1.2, 99.8, 109.0, 96.6, 112.5),
        AngleDefinition("pucker==DA_DG_other", "C1'", "O4'", "C4'", 0, 0, 0, 107.3, 1.5, 27049, 108.7, 1.6, 102.1, 112.5, 95.4, 115.2),
        AngleDefinition("pucker==DA_DG_other", "C2'", "C1'", "O4'", 0, 0, 0, 105.7, 1.1, 27049, 106.0, 1.4, 100.4, 109.8, 95.5, 112.4),
        AngleDefinition("pucker==DA_DG_other", "C2'", "C3'", "O3'", 0, 0, 0, 112.3, 1.8, 26998, 112.1, 2.1, 101.5, 122.6, 92.5, 128.0),
        AngleDefinition("pucker==DA_DG_other", "C4'", "C3'", "O3'", 0, 0, 0, 110.2, 1.8, 26998, 109.4, 2.1, 100.1, 120.5, 91.6, 126.5),
        AngleDefinition("pucker==DA_DG_other", "C3'", "C4'", "C5'", 0, 0, 0, 114.9, 1.5, 27000, 114.8, 1.4, 106.8, 122.6, 97.9, 128.1),
        AngleDefinition("pucker==DA_DG_other", "C5'", "C4'", "O4'", 0, 0, 0, 108.7, 1.4, 27000, 109.6, 1.4, 103.1, 117.6, 99.4, 121.1),
        AngleDefinition("pucker==DA_DG_other", "N9 ", "C1'", "O4'", 0, 0, 0, 107.7, 1.3, 27034, 108.5, 1.5, 101.6, 116.6, 95.6, 126.4),
        AngleDefinition("pucker==DA_DG_other", "N9 ", "C1'", "C2'", 0, 0, 0, 115.4, 0.9, 27034, 114.2, 1.8, 106.1, 121.8, 97.9, 126.0),
        AngleDefinition("pucker==DA_DG_other", "C4'", "C5'", "O5'", 0, 0, 0, 110.9, 2.1, 26939, 110.5, 1.6, 103.4, 117.7, 98.5, 140.9),
        AngleDefinition("pucker==DA_DG_other", "C1'", "N9", "C4", 0, 0, 0, 127.6, 1.9, 27026, 126.6, 1.0, 121.6, 130.6, 114.5, 139.6),
        AngleDefinition("pucker==DA_DG_other", "C1'", "N9", "C8", 0, 0, 0, 126.2, 1.9, 27026, 127.4, 1.0, 123.5, 132.1, 114.2, 141.6),
    ],
    "pucker==DU_DT_DC_C2p_endo": [
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 102.4, 0.9, 49389, 101.3, 1.5, 95.6, 106.0, 91.8, 107.8
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.9, 0.9, 49389, 103.1, 0.9, 98.9, 106.4, 95.9, 108.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 105.9, 0.8, 49389, 105.9, 0.7, 102.6, 108.5, 98.7, 110.3
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.0, 0.6, 49389, 109.4, 1.0, 105.4, 112.4, 102.1, 115.0
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 106.2, 0.7, 49389, 105.6, 1.0, 101.4, 108.7, 97.8, 110.7
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 109.5, 2.3, 49347, 111.5, 2.0, 101.6, 120.1, 94.9, 124.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 109.8, 2.0, 49347, 109.2, 1.8, 101.0, 117.2, 96.4, 125.4
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 114.7, 1.4, 49335, 114.7, 1.4, 108.2, 120.7, 102.8, 126.1
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.7, 1.1, 49335, 109.7, 1.1, 104.1, 115.8, 98.8, 122.5
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "N1", "C1'", "O4'", 0, 0, 0, 107.8, 0.6, 49342, 108.4, 1.4, 101.7, 114.6, 91.8, 119.6
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "N1", "C1'", "C2'", 0, 0, 0, 114.0, 0.9, 49342, 114.9, 1.8, 108.0, 121.6, 103.7, 133.1
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 110.1, 1.7, 49228, 110.6, 1.2, 105.3, 116.5, 98.7, 122.7
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C1'", "N1", "C2", 0, 0, 0, 118.3, 1.3, 49331, 119.0, 0.9, 115.0, 122.6, 108.0, 127.8
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C2p_endo", "C1'", "N1", "C6", 0, 0, 0, 120.3, 1.1, 49331, 120.1, 0.9, 116.4, 123.8, 110.6, 132.2
        ),
    ],
    "pucker==DU_DT_DC_C3p_endo": [
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C1'", "C2'", "C3'", 0, 0, 0, 103.0, 1.1, 8671, 102.5, 1.3, 97.5, 106.5, 95.9, 108.4
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", "C4'", 0, 0, 0, 102.2, 0.8, 8671, 102.1, 1.4, 96.4, 106.4, 91.2, 107.4
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", "O4'", 0, 0, 0, 104.8, 0.8, 8671, 104.5, 1.1, 99.6, 107.8, 93.9, 109.3
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C1'", "O4'", "C4'", 0, 0, 0, 110.1, 0.7, 8671, 109.6, 0.9, 106.0, 112.5, 104.5, 115.6
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C2'", "C1'", "O4'", 0, 0, 0, 107.1, 0.5, 8671, 107.3, 0.7, 103.6, 110.6, 101.1, 112.0
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", "O3'", 0, 0, 0, 112.6, 2.7, 8626, 111.7, 2.0, 101.8, 122.3, 96.4, 127.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C4'", "C3'", "O3'", 0, 0, 0, 112.1, 2.0, 8626, 110.2, 1.8, 103.2, 121.0, 99.3, 131.5
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", "C5'", 0, 0, 0, 116.1, 1.1, 8662, 114.9, 1.2, 108.6, 120.5, 102.4, 122.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C5'", "C4'", "O4'", 0, 0, 0, 109.5, 0.7, 8662, 109.6, 1.1, 104.1, 115.5, 98.9, 121.2
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "N1", "C1'", "O4'", 0, 0, 0, 108.1, 0.8, 8654, 108.8, 1.2, 103.9, 114.6, 97.7, 118.2
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "N1", "C1'", "C2'", 0, 0, 0, 113.0, 1.3, 8654, 113.9, 1.5, 106.6, 121.4, 101.1, 127.1
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C4'", "C5'", "O5'", 0, 0, 0, 111.1, 2.5, 8651, 110.5, 1.3, 104.6, 116.8, 100.5, 118.2
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C1'", "N1", "C2", 0, 0, 0, 116.6, 1.1, 8650, 118.9, 1.1, 114.1, 122.9, 109.8, 125.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_C3p_endo", "C1'", "N1", "C6", 0, 0, 0, 122.0, 0.9, 8650, 120.3, 1.0, 116.7, 124.6, 112.9, 129.5
        ),
    ],
    "pucker==DU_DT_DC_other": [
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C1'", "C2'", "C3'", 0, 0, 0, 103.9, 1.5, 42575, 103.0, 1.4, 97.8, 107.1, 93.7, 108.5
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C2'", "C3'", "C4'", 0, 0, 0, 103.5, 1.3, 42575, 104.1, 1.1, 98.6, 107.6, 94.3, 109.5
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C3'", "C4'", "O4'", 0, 0, 0, 106.1, 1.5, 42575, 105.6, 1.1, 100.0, 108.6, 94.6, 110.5
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C1'", "O4'", "C4'", 0, 0, 0, 109.1, 0.9, 42575, 108.3, 1.5, 102.1, 112.2, 95.7, 114.1
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C2'", "C1'", "O4'", 0, 0, 0, 106.3, 1.1, 42575, 105.8, 1.3, 100.6, 109.7, 96.2, 111.3
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C2'", "C3'", "O3'", 0, 0, 0, 111.6, 2.6, 42529, 112.1, 1.9, 102.1, 121.5, 92.2, 127.7
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C4'", "C3'", "O3'", 0, 0, 0, 111.5, 1.4, 42529, 109.3, 2.0, 100.6, 119.2, 95.2, 128.0
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C3'", "C4'", "C5'", 0, 0, 0, 115.1, 1.2, 42539, 114.8, 1.3, 108.4, 121.1, 102.5, 127.2
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C5'", "C4'", "O4'", 0, 0, 0, 108.7, 1.0, 42539, 109.7, 1.2, 104.3, 115.7, 98.6, 122.2
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "N1", "C1'", "O4'", 0, 0, 0, 107.7, 0.8, 42520, 108.4, 1.5, 101.2, 114.8, 96.5, 119.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "N1", "C1'", "C2'", 0, 0, 0, 114.5, 1.7, 42520, 114.6, 1.8, 106.9, 122.3, 101.8, 129.1
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C4'", "C5'", "O5'", 0, 0, 0, 109.7, 2.3, 42490, 110.5, 1.4, 104.1, 117.3, 99.9, 123.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C1'", "N1", "C2", 0, 0, 0, 119.0, 1.7, 42502, 119.1, 0.9, 115.1, 123.9, 109.5, 130.9
        ),
        AngleDefinition(
            "pucker==DU_DT_DC_other", "C1'", "N1", "C6", 0, 0, 0, 119.4, 1.6, 42502, 120.0, 0.9, 115.5, 123.9, 108.3, 133.8
        ),
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
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

    def _atom_names_angles(self, res_name: str) -> List[AngleDefinition]:
        if res_name in ("A", "G"):
            return self.angles_definition["pucker==A_G_C2p_endo"]
        if res_name in ("U", "T", "C"):
            return self.angles_definition["pucker==U_T_C_C2p_endo"]
        if res_name in ("DA", "DG"):
            return self.angles_definition["pucker==DA_DG_C2p_endo"]
        if res_name in ("DU", "DT", "DC"):
            return self.angles_definition["pucker==DU_DT_DC_C2p_endo"]
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

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
            raise NonStandardResidueException(f"Non-standard residue: {res_name}")

        if sugar_conformation == "C3'-endo":
            if res_name in ("A", "G"):
                return self.bonds_definition["pucker==A_G_C3p_endo"]
            if res_name in ("U", "T", "C"):
                return self.bonds_definition["pucker==U_T_C_C3p_endo"]
            if res_name in ("DA", "DG"):
                return self.bonds_definition["pucker==DA_DG_C3p_endo"]
            if res_name in ("DU", "DT", "DC"):
                return self.bonds_definition["pucker==DU_DT_DC_C3p_endo"]
            raise NonStandardResidueException(f"Non-standard residue: {res_name}")

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
            raise NonStandardResidueException(f"Non-standard residue: {res_name}")

        if res_name in ("A", "G"):
            return self.basic_bonds_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.basic_bonds_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.basic_bonds_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.basic_bonds_definition["sugar_basic==DU_DT_DC"]

        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

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
            raise NonStandardResidueException(f"Non-standard residue: {res_name}")

        if sugar_conformation == "C3'-endo":
            if res_name in ("A", "G"):
                return self.angles_definition["pucker==A_G_C3p_endo"]
            if res_name in ("U", "T", "C"):
                return self.angles_definition["pucker==U_T_C_C3p_endo"]
            if res_name in ("DA", "DG"):
                return self.angles_definition["pucker==DA_DG_C3p_endo"]
            if res_name in ("DU", "DT", "DC"):
                return self.angles_definition["pucker==DU_DT_DC_C3p_endo"]
            raise NonStandardResidueException(f"Non-standard residue: {res_name}")

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
            raise NonStandardResidueException(f"Non-standard residue: {res_name}")

        if res_name in ("A", "G"):
            return self.basic_angles_definition["sugar_basic==A_G"]
        if res_name in ("U", "T", "C"):
            return self.basic_angles_definition["sugar_basic==U_T_C"]
        if res_name in ("DA", "DG"):
            return self.basic_angles_definition["sugar_basic==DA_DG"]
        if res_name in ("DU", "DT", "DC"):
            return self.basic_angles_definition["sugar_basic==DU_DT_DC"]

        raise NonStandardResidueException(f"Non-standard residue: {res_name}")
