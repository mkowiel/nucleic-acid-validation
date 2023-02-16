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
        BondDefinition("pucker==DU_DT_DC_C2p_endo", "C1'", "N1",  0, 0, 1.474, 0.013, 49340, 1.488, 0.016, 1.430, 1.554, 1.351, 1.621),
    ],
    "pucker==DU_DT_DC_C3p_endo": [
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "C2'", 0, 0, 1.518, 0.010, 8671, 1.527, 0.009, 1.478, 1.560, 1.358, 1.594),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C2'", "C3'", 0, 0, 1.518, 0.010, 8671, 1.523, 0.008, 1.491, 1.562, 1.283, 1.607),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C3'", "C4'", 0, 0, 1.518, 0.008, 8671, 1.518, 0.008, 1.483, 1.553, 1.467, 1.649),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "O4'", 0, 0, 1.446, 0.009, 8671, 1.446, 0.008, 1.411, 1.478, 1.288, 1.511),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "O4'", 0, 0, 1.413, 0.009, 8671, 1.420, 0.008, 1.387, 1.463, 1.372, 1.613),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C4'", "C5'", 0, 0, 1.509, 0.011, 8662, 1.515, 0.007, 1.486, 1.549, 1.477, 1.576),
        BondDefinition("pucker==DU_DT_DC_C3p_endo", "C1'", "N1",  0, 0, 1.492, 0.010, 8654, 1.491, 0.017, 1.438, 1.560, 1.367, 1.607),
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
