from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.validator import NonStandardResidueException, Validator

# TODO: maybe  "O3'", "C3'" relative positions should be -1, -1
PO4_BONDS = {
    "PO4==AA_0": [
        BondDefinition("PO4==AA_0", "OP1", "P", 0, 0, 1.487, 0.01, 272554, 1.485, 0.007, 1.448, 1.519, 1.381, 1.557),
        BondDefinition("PO4==AA_0", "OP2", "P", 0, 0, 1.483, 0.01, 272566, 1.486, 0.008, 1.444, 1.524, 1.359, 1.566),
        BondDefinition("PO4==AA_0", "O3'", "P", -1, 0, 1.580, 0.010, 272526, 1.605, 0.008, 1.561, 1.640, 1.491, 1.682),
        BondDefinition("PO4==AA_0", "O5'", "P", 0, 0, 1.603, 0.011, 272573, 1.591, 0.007, 1.558, 1.621, 1.533, 1.656),
        BondDefinition("PO4==AA_0", "O3'", "C3'", 0, 0, 1.422, 0.010, 272573, 1.418, 0.010, 1.382, 1.469, 1.346, 1.522),
        BondDefinition("PO4==AA_0", "O5'", "C5'", 0, 0, 1.428, 0.013, 272573, 1.418, 0.009, 1.370, 1.456, 1.284, 1.499),
    ],
    "PO4==AA_1": [
        BondDefinition("PO4==AA_1", "OP1", "P", 0, 0, 1.483, 0.010, 127494, 1.484, 0.008, 1.436, 1.527, 1.377, 1.568),
        BondDefinition("PO4==AA_1", "OP2", "P", 0, 0, 1.487, 0.010, 127501, 1.485, 0.008, 1.439, 1.529, 1.386, 1.586),
        BondDefinition("PO4==AA_1", "O3'", "P", -1, 0, 1.603, 0.011, 127498, 1.604, 0.010, 1.548, 1.641, 1.479, 1.682),
        BondDefinition("PO4==AA_1", "O5'", "P", 0, 0, 1.580, 0.010, 127502, 1.592, 0.008, 1.557, 1.624, 1.528, 1.675),
        BondDefinition("PO4==AA_1", "O3'", "C3'", 0, 0, 1.422, 0.010, 127502, 1.420, 0.011, 1.366, 1.465, 1.330, 1.522),
        BondDefinition("PO4==AA_1", "O5'", "C5'", 0, 0, 1.428, 0.013, 127502, 1.422, 0.009, 1.380, 1.462, 1.340, 1.502),
    ],
    "PO4==AA_2": [
        BondDefinition("PO4==AA_2", "OP1", "P", 0, 0, 1.487, 0.010, 81443, 1.485, 0.008, 1.440, 1.525, 1.376, 1.570),
        BondDefinition("PO4==AA_2", "OP2", "P", 0, 0, 1.483, 0.010, 81449, 1.485, 0.009, 1.435, 1.530, 1.372, 1.569),
        BondDefinition("PO4==AA_2", "O3'", "P", -1, 0, 1.603, 0.011, 81419, 1.605, 0.009, 1.557, 1.642, 1.505, 1.679),
        BondDefinition("PO4==AA_2", "O5'", "P", 0, 0, 1.580, 0.010, 81453, 1.591, 0.008, 1.552, 1.624, 1.521, 1.655),
        BondDefinition("PO4==AA_2", "O3'", "C3'", 0, 0, 1.422, 0.010, 81453, 1.420, 0.010, 1.377, 1.464, 1.315, 1.495),
        BondDefinition("PO4==AA_2", "O5'", "C5'", 0, 0, 1.428, 0.013, 81453, 1.420, 0.009, 1.354, 1.459, 1.259, 1.517),
    ],
    "PO4==AA_3": [
        BondDefinition("PO4==AA_3", "OP1", "P", 0, 0, 1.483, 0.010, 57134, 1.485, 0.009, 1.439, 1.526, 1.328, 1.604),
        BondDefinition("PO4==AA_3", "OP2", "P", 0, 0, 1.487, 0.010, 57138, 1.485, 0.008, 1.438, 1.527, 1.342, 1.629),
        BondDefinition("PO4==AA_3", "O3'", "P", -1, 0, 1.580, 0.010, 57125, 1.606, 0.010, 1.554, 1.648, 1.486, 1.725),
        BondDefinition("PO4==AA_3", "O5'", "P", 0, 0, 1.603, 0.011, 57138, 1.592, 0.008, 1.557, 1.623, 1.524, 1.655),
        BondDefinition("PO4==AA_3", "O3'", "C3'", 0, 0, 1.422, 0.010, 57138, 1.423, 0.012, 1.382, 1.489, 1.340, 1.519),
        BondDefinition("PO4==AA_3", "O5'", "C5'", 0, 0, 1.428, 0.013, 57138, 1.421, 0.008, 1.384, 1.460, 1.338, 1.488),
    ],
    "PO4==AS_0": [
        BondDefinition("PO4==AS_0", "OP1", "P", 0, 0, 1.484, 0.012, 1808496, 1.485, 0.007, 1.442, 1.520, 1.392, 1.597),
        BondDefinition("PO4==AS_0", "OP2", "P", 0, 0, 1.478, 0.010, 1808549, 1.485, 0.008, 1.441, 1.522, 1.382, 1.605),
        BondDefinition("PO4==AS_0", "O3'", "P", -1, 0, 1.599, 0.016, 1808621, 1.604, 0.008, 1.560, 1.635, 1.498, 1.676),
        BondDefinition("PO4==AS_0", "O5'", "P", 0, 0, 1.601, 0.016, 1808644, 1.592, 0.007, 1.561, 1.619, 1.537, 1.646),
        BondDefinition("PO4==AS_0", "O3'", "C3'", 0, 0, 1.438, 0.007, 1808644, 1.415, 0.008, 1.379, 1.454, 1.345, 1.495),
        BondDefinition("PO4==AS_0", "O5'", "C5'", 0, 0, 1.437, 0.017, 1808644, 1.419, 0.007, 1.386, 1.455, 1.345, 1.486),
    ],
    "PO4==AS_1": [
        BondDefinition("PO4==AS_1", "OP1", "P", 0, 0, 1.478, 0.010, 1808496, 1.485, 0.007, 1.442, 1.520, 1.392, 1.597),
        BondDefinition("PO4==AS_1", "OP2", "P", 0, 0, 1.484, 0.012, 1808549, 1.485, 0.008, 1.441, 1.522, 1.382, 1.605),
        BondDefinition("PO4==AS_1", "O3'", "P", -1, 0, 1.601, 0.016, 1808621, 1.604, 0.008, 1.560, 1.635, 1.498, 1.676),
        BondDefinition("PO4==AS_1", "O5'", "P", 0, 0, 1.599, 0.016, 1808644, 1.592, 0.007, 1.561, 1.619, 1.537, 1.646),
        BondDefinition("PO4==AS_1", "O3'", "C3'", 0, 0, 1.438, 0.007, 1808644, 1.415, 0.008, 1.379, 1.454, 1.345, 1.495),
        BondDefinition("PO4==AS_1", "O5'", "C5'", 0, 0, 1.437, 0.017, 1808644, 1.419, 0.007, 1.386, 1.455, 1.345, 1.486),
    ],
    "PO4==AS_2": [
        BondDefinition("PO4==AS_2", "OP1", "P", 0, 0, 1.484, 0.012, 70389, 1.485, 0.009, 1.433, 1.528, 1.350, 1.579),
        BondDefinition("PO4==AS_2", "OP2", "P", 0, 0, 1.478, 0.010, 70393, 1.484, 0.009, 1.436, 1.525, 1.348, 1.578),
        BondDefinition("PO4==AS_2", "O3'", "P", -1, 0, 1.601, 0.016, 70389, 1.606, 0.009, 1.560, 1.641, 1.471, 1.694),
        BondDefinition("PO4==AS_2", "O5'", "P", 0, 0, 1.599, 0.016, 70395, 1.592, 0.008, 1.556, 1.625, 1.530, 1.674),
        BondDefinition("PO4==AS_2", "O3'", "C3'", 0, 0, 1.438, 0.007, 70395, 1.424, 0.010, 1.383, 1.476, 1.354, 1.507),
        BondDefinition("PO4==AS_2", "O5'", "C5'", 0, 0, 1.437, 0.017, 70395, 1.420, 0.008, 1.379, 1.459, 1.331, 1.502),
    ],
    "PO4==AS_3": [
        BondDefinition("PO4==AS_3", "OP1", "P", 0, 0, 1.478, 0.010, 70389, 1.485, 0.009, 1.433, 1.528, 1.350, 1.579),
        BondDefinition("PO4==AS_3", "OP2", "P", 0, 0, 1.484, 0.012, 70393, 1.484, 0.009, 1.436, 1.525, 1.348, 1.578),
        BondDefinition("PO4==AS_3", "O3'", "P", -1, 0, 1.599, 0.016, 70389, 1.606, 0.009, 1.560, 1.641, 1.471, 1.694),
        BondDefinition("PO4==AS_3", "O5'", "P", 0, 0, 1.601, 0.016, 70395, 1.592, 0.008, 1.556, 1.625, 1.530, 1.674),
        BondDefinition("PO4==AS_3", "O3'", "C3'", 0, 0, 1.438, 0.007, 70395, 1.424, 0.010, 1.383, 1.476, 1.354, 1.507),
        BondDefinition("PO4==AS_3", "O5'", "C5'", 0, 0, 1.437, 0.017, 70395, 1.420, 0.008, 1.379, 1.459, 1.331, 1.502),
    ],
    "other==A_G": [
        BondDefinition("other==A_G", "OP1", "P", 0, 0, 1.485, 0.017, 1492318, 1.485, 0.008, 1.443, 1.522, 1.378, 1.599),
        BondDefinition("other==A_G", "OP2", "P", 0, 0, 1.485, 0.017, 1492363, 1.485, 0.008, 1.442, 1.524, 1.374, 1.606),
        BondDefinition("other==A_G", "O3'", "P", -1, 0, 1.607, 0.012, 1479008, 1.605, 0.008, 1.564, 1.638, 1.495, 1.686),
        BondDefinition("other==A_G", "O5'", "P", 0, 0, 1.593, 0.010, 1492383, 1.592, 0.007, 1.559, 1.620, 1.535, 1.654),
        BondDefinition("other==A_G", "O3'", "C3'", 0, 0, 1.420, 0.009, 1493987, 1.417, 0.009, 1.382, 1.460, 1.355, 1.499),
        BondDefinition("other==A_G", "O5'", "C5'", 0, 0, 1.426, 0.014, 1493976, 1.420, 0.007, 1.387, 1.455, 1.346, 1.484),
    ],
    "other==DA_DG": [
        BondDefinition("other==DA_DG", "OP1", "P", 0, 0, 1.485, 0.017, 97138, 1.482, 0.008, 1.426, 1.527, 1.385, 1.578),
        BondDefinition("other==DA_DG", "OP2", "P", 0, 0, 1.485, 0.017, 97157, 1.483, 0.009, 1.424, 1.532, 1.363, 1.586),
        BondDefinition("other==DA_DG", "O3'", "P", -1, 0, 1.607, 0.012, 94948, 1.603, 0.013, 1.520, 1.655, 1.384, 1.747),
        BondDefinition("other==DA_DG", "O5'", "P", 0, 0, 1.593, 0.010, 97146, 1.596, 0.009, 1.557, 1.634, 1.448, 1.719),
        BondDefinition("other==DA_DG", "O3'", "C3'", 0, 0, 1.427, 0.011, 100888, 1.423, 0.015, 1.356, 1.489, 1.308, 1.568),
        BondDefinition("other==DA_DG", "O5'", "C5'", 0, 0, 1.424, 0.011, 100790, 1.424, 0.015, 1.343, 1.476, 1.257, 1.538),
    ],
    "other==U_T_C": [
        BondDefinition("other==U_T_C", "OP1", "P", 0, 0, 1.485, 0.017, 1126572, 1.485, 0.007, 1.444, 1.521, 1.380, 1.569),
        BondDefinition("other==U_T_C", "OP2", "P", 0, 0, 1.485, 0.017, 1126601, 1.485, 0.008, 1.442, 1.521, 1.374, 1.565),
        BondDefinition("other==U_T_C", "O3'", "P", -1, 0, 1.607, 0.012, 1118723, 1.605, 0.008, 1.564, 1.638, 1.496, 1.683),
        BondDefinition("other==U_T_C", "O5'", "P", 0, 0, 1.593, 0.010, 1126582, 1.592, 0.007, 1.560, 1.619, 1.537, 1.648),
        BondDefinition("other==U_T_C", "O3'", "C3'", 0, 0, 1.420, 0.011, 1127583, 1.416, 0.009, 1.381, 1.457, 1.350, 1.497),
        BondDefinition("other==U_T_C", "O5'", "C5'", 0, 0, 1.429, 0.012, 1127502, 1.419, 0.007, 1.388, 1.454, 1.343, 1.481),
    ],
    "other==DU_DT_DC": [
        BondDefinition("other==DU_DT_DC", "OP1", "P", 0, 0, 1.485, 0.017, 96681, 1.482, 0.008, 1.425, 1.530, 1.380, 1.604),
        BondDefinition("other==DU_DT_DC", "OP2", "P", 0, 0, 1.485, 0.017, 96687, 1.483, 0.009, 1.423, 1.534, 1.364, 1.612),
        BondDefinition("other==DU_DT_DC", "O3'", "P", -1, 0, 1.607, 0.012, 94757, 1.604, 0.012, 1.524, 1.652, 1.431, 1.715),
        BondDefinition("other==DU_DT_DC", "O5'", "P", 0, 0, 1.593, 0.010, 96663, 1.596, 0.008, 1.558, 1.636, 1.502, 1.688),
        BondDefinition("other==DU_DT_DC", "O3'", "C3'", 0, 0, 1.430, 0.015, 100519, 1.423, 0.014, 1.358, 1.487, 1.305, 1.544),
        BondDefinition("other==DU_DT_DC", "O5'", "C5'", 0, 0, 1.425, 0.014, 100418, 1.423, 0.015, 1.345, 1.478, 1.256, 1.535),
    ],
}


PO4_ANGLES = {
    "PO4==AA_0": [
        AngleDefinition("PO4==AA_0", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 272547, 119.6, 1.5, 111.8, 127.8, 104.4, 133.1),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O3'", 0, 0, -1, 106.2, 1.1, 272554, 107.7, 2.0, 97.7, 116.6, 48.2, 155.2),
        AngleDefinition("PO4==AA_0", "OP1", "P", "O5'", 0, 0, 0, 110.2, 1.3, 272554, 108.3, 1.9, 99.0, 116.6, 91.7, 123.2),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O3'", 0, 0, -1, 112.2, 1.0, 272566, 108.5, 2.0, 100.1, 118.0, 24.4, 126.3),
        AngleDefinition("PO4==AA_0", "OP2", "P", "O5'", 0, 0, 0, 109.3, 0.9, 272566, 108.2, 2.0, 98.5, 117.0, 90.0, 123.8),
        AngleDefinition("PO4==AA_0", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 272573, 103.3, 1.3, 96.8, 108.5, 48.6, 123.6),
        AngleDefinition("PO4==AA_0", "P", "O3'", "C3'", 0, -1, -1, 120.8, 1.1, 272573, 120.2, 1.8, 113.3, 128.7, 34.7, 137.6),
        AngleDefinition("PO4==AA_0", "P", "O5'", "C5'", 0, 0, 0, 120.3, 1.9, 272573, 120.9, 1.0, 116.0, 125.3, 110.2, 130.3),
    ],
    "PO4==AA_1": [
        AngleDefinition("PO4==AA_1", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 127493, 119.8, 1.6, 111.7, 128.0, 96.9, 133.9),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O3'", 0, 0, -1, 109.3, 0.9, 127494, 108.2, 2.1, 98.9, 118.4, 82.8, 137.4),
        AngleDefinition("PO4==AA_1", "OP1", "P", "O5'", 0, 0, 0, 112.2, 1.0, 127494, 108.3, 2.2, 97.8, 118.1, 75.1, 129.6),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O3'", 0, 0, -1, 110.2, 1.3, 127501, 108.1, 2.0, 98.9, 117.9, 82.8, 135.2),
        AngleDefinition("PO4==AA_1", "OP2", "P", "O5'", 0, 0, 0, 106.2, 1.1, 127501, 107.9, 2.2, 96.6, 117.0, 85.6, 132.3),
        AngleDefinition("PO4==AA_1", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 127502, 103.2, 1.2, 96.5, 107.4, 79.4, 114.9),
        AngleDefinition("PO4==AA_1", "P", "O3'", "C3'", 0, -1, -1, 120.3, 1.9, 127502, 120.7, 1.6, 114.1, 130.4, 85.4, 138.5),
        AngleDefinition("PO4==AA_1", "P", "O5'", "C5'", 0, 0, 0, 120.8, 1.1, 127502, 120.6, 1.1, 114.0, 125.9, 108.4, 132.7),
    ],
    "PO4==AA_2": [
        AngleDefinition("PO4==AA_2", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 81439, 119.9, 1.7, 112.1, 129.1, 101.5, 136.8),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O3'", 0, 0, -1, 110.2, 1.3, 81443, 108.2, 2.2, 99.4, 118.5, 31.1, 149.0),
        AngleDefinition("PO4==AA_2", "OP1", "P", "O5'", 0, 0, 0, 106.2, 1.1, 81443, 107.9, 2.2, 95.7, 117.1, 88.5, 132.6),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O3'", 0, 0, -1, 109.3, 0.9, 81449, 108.1, 2.1, 98.7, 118.3, 71.4, 144.5),
        AngleDefinition("PO4==AA_2", "OP2", "P", "O5'", 0, 0, 0, 112.2, 1.0, 81449, 108.1, 2.2, 97.1, 117.9, 76.5, 123.5),
        AngleDefinition("PO4==AA_2", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 81453, 103.3, 1.4, 96.0, 107.9, 50.6, 122.0),
        AngleDefinition("PO4==AA_2", "P", "O3'", "C3'", 0, -1, -1, 120.3, 1.9, 81453, 120.5, 1.9, 112.5, 130.3, 46.1, 139.2),
        AngleDefinition("PO4==AA_2", "P", "O5'", "C5'", 0, 0, 0, 120.8, 1.1, 81453, 120.6, 1.1, 113.6, 125.4, 107.4, 130.5),
    ],
    "PO4==AA_3": [
        AngleDefinition("PO4==AA_3", "OP1", "P", "OP2", 0, 0, 0, 117.6, 1.2, 57134, 119.5, 1.7, 110.9, 127.9, 103.1, 134.2),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O3'", 0, 0, -1, 112.2, 1.0, 57134, 108.5, 2.2, 99.0, 118.1, 72.2, 134.0),
        AngleDefinition("PO4==AA_3", "OP1", "P", "O5'", 0, 0, 0, 109.3, 0.9, 57134, 108.6, 2.1, 98.0, 118.2, 88.8, 126.4),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O3'", 0, 0, -1, 106.2, 1.1, 57138, 107.4, 2.1, 97.2, 116.9, 67.4, 133.7),
        AngleDefinition("PO4==AA_3", "OP2", "P", "O5'", 0, 0, 0, 110.2, 1.3, 57138, 107.8, 2.2, 97.7, 117.1, 89.0, 123.6),
        AngleDefinition("PO4==AA_3", "O3'", "P", "O5'", -1, 0, 0, 99.9, 0.7, 57138, 103.9, 1.5, 97.5, 110.6, 65.5, 148.5),
        AngleDefinition("PO4==AA_3", "P", "O3'", "C3'", 0, -1, -1, 120.8, 1.1, 57138, 121.1, 2.0, 114.1, 131.6, 72.0, 142.2),
        AngleDefinition("PO4==AA_3", "P", "O5'", "C5'", 0, 0, 0, 120.3, 1.9, 57138, 120.7, 1.2, 113.7, 125.8, 107.0, 130.6),
    ],
    "PO4==AS_0": [
        AngleDefinition("PO4==AS_0", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 1808403, 119.7, 1.5, 111.9, 127.5, 103.3, 133.7),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O3'", 0, 0, -1, 104.5, 0.9, 1808496, 107.6, 2.0, 97.3, 116.8, 86.2, 124.1),
        AngleDefinition("PO4==AS_0", "OP1", "P", "O5'", 0, 0, 0, 110.3, 0.8, 1808496, 108.4, 2.0, 99.0, 117.2, 92.1, 124.2),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O3'", 0, 0, -1, 111.5, 1.1, 1808549, 108.4, 1.9, 100.0, 118.1, 90.7, 124.8),
        AngleDefinition("PO4==AS_0", "OP2", "P", "O5'", 0, 0, 0, 105.2, 0.8, 1808549, 108.0, 2.1, 98.0, 117.2, 90.3, 124.4),
        AngleDefinition("PO4==AS_0", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 1808644, 103.6, 1.0, 98.9, 107.3, 87.9, 115.0),
        AngleDefinition("PO4==AS_0", "P", "O3'", "C3'", 0, -1, -1, 121.5, 3.0, 1808644, 120.0, 1.2, 112.8, 126.1, 105.5, 133.0),
        AngleDefinition("PO4==AS_0", "P", "O5'", "C5'", 0, 0, 0, 121.6, 2.8, 1808644, 120.5, 1.0, 115.3, 124.5, 110.2, 129.5),
    ],
    "PO4==AS_1": [
        AngleDefinition("PO4==AS_1", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 1808403, 119.7, 1.5, 111.9, 127.5, 103.3, 133.7),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O3'", 0, 0, -1, 105.2, 0.8, 1808496, 107.6, 2.0, 97.3, 116.8, 86.2, 124.1),
        AngleDefinition("PO4==AS_1", "OP1", "P", "O5'", 0, 0, 0, 111.5, 1.1, 1808496, 108.4, 2.0, 99.0, 117.2, 92.1, 124.2),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O3'", 0, 0, -1, 110.3, 0.8, 1808549, 108.4, 1.9, 100.0, 118.1, 90.7, 124.8),
        AngleDefinition("PO4==AS_1", "OP2", "P", "O5'", 0, 0, 0, 104.5, 0.9, 1808549, 108.0, 2.1, 98.0, 117.2, 90.3, 124.4),
        AngleDefinition("PO4==AS_1", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 1808644, 103.6, 1.0, 98.9, 107.3, 87.9, 115.0),
        AngleDefinition("PO4==AS_1", "P", "O3'", "C3'", 0, -1, -1, 121.6, 2.8, 1808644, 120.0, 1.2, 112.8, 126.1, 105.5, 133.0),
        AngleDefinition("PO4==AS_1", "P", "O5'", "C5'", 0, 0, 0, 121.5, 3.0, 1808644, 120.5, 1.0, 115.3, 124.5, 110.2, 129.5),
    ],
    "PO4==AS_2": [
        AngleDefinition("PO4==AS_2", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 70387, 119.8, 1.8, 111.4, 128.6, 103.8, 133.9),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O3'", 0, 0, -1, 110.3, 0.8, 70389, 108.7, 2.3, 99.8, 119.4, 86.3, 132.3),
        AngleDefinition("PO4==AS_2", "OP1", "P", "O5'", 0, 0, 0, 104.5, 0.9, 70389, 108.0, 2.4, 95.7, 117.4, 87.5, 125.1),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O3'", 0, 0, -1, 105.2, 0.8, 70393, 107.1, 2.2, 96.1, 116.7, 74.3, 133.2),
        AngleDefinition("PO4==AS_2", "OP2", "P", "O5'", 0, 0, 0, 111.5, 1.1, 70393, 108.1, 2.3, 97.3, 117.8, 88.6, 125.6),
        AngleDefinition("PO4==AS_2", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 70395, 103.9, 1.1, 99.3, 109.2, 84.2, 118.6),
        AngleDefinition("PO4==AS_2", "P", "O3'", "C3'", 0, -1, -1, 121.6, 2.8, 70395, 120.8, 1.6, 115.0, 129.8, 101.4, 138.8),
        AngleDefinition("PO4==AS_2", "P", "O5'", "C5'", 0, 0, 0, 121.5, 3.0, 70395, 120.8, 1.0, 115.6, 125.3, 110.4, 130.3),
    ],
    "PO4==AS_3": [
        AngleDefinition("PO4==AS_3", "OP1", "P", "OP2", 0, 0, 0, 119.9, 1.6, 70387, 119.8, 1.8, 111.4, 128.6, 103.8, 133.9),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O3'", 0, 0, -1, 111.5, 1.1, 70389, 108.7, 2.3, 99.8, 119.4, 86.3, 132.3),
        AngleDefinition("PO4==AS_3", "OP1", "P", "O5'", 0, 0, 0, 105.2, 0.8, 70389, 108.0, 2.4, 95.7, 117.4, 87.5, 125.1),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O3'", 0, 0, -1, 104.5, 0.9, 70393, 107.1, 2.2, 96.1, 116.7, 74.3, 133.2),
        AngleDefinition("PO4==AS_3", "OP2", "P", "O5'", 0, 0, 0, 110.3, 0.8, 70393, 108.1, 2.3, 97.3, 117.8, 88.6, 125.6),
        AngleDefinition("PO4==AS_3", "O3'", "P", "O5'", -1, 0, 0, 104.2, 1.5, 70395, 103.9, 1.1, 99.3, 109.2, 84.2, 118.6),
        AngleDefinition("PO4==AS_3", "P", "O3'", "C3'", 0, -1, -1, 121.5, 3.0, 70395, 120.8, 1.6, 115.0, 129.8, 101.4, 138.8),
        AngleDefinition("PO4==AS_3", "P", "O5'", "C5'", 0, 0, 0, 121.6, 2.8, 70395, 120.8, 1.0, 115.6, 125.3, 110.4, 130.3),
    ],
    "other": [
        AngleDefinition("other", "OP1", "P", "OP2", 0, 0, 0, 119.6, 1.5, 2812583, 119.7, 1.6, 111.7, 127.7, 102.8, 134.2),
        AngleDefinition("other", "OP1", "P", "O3'", 0, 0, -1, 107.7, 3.2, 2787628, 107.7, 2.1, 97.4, 117.4, 61.4, 142.7),
        AngleDefinition("other", "OP1", "P", "O5'", 0, 0, 0, 108.1, 2.9, 2812541, 108.3, 2.0, 98.5, 117.3, 89.8, 125.4),
        AngleDefinition("other", "OP2", "P", "O3'", 0, 0, -1, 108.3, 3.2, 2787726, 108.3, 2.0, 99.4, 118.1, 62.4, 136.8),
        AngleDefinition("other", "OP2", "P", "O5'", 0, 0, 0, 108.3, 2.7, 2812640, 108.0, 2.1, 97.9, 117.2, 88.8, 124.9),
        AngleDefinition("other", "O3'", "P", "O5'", -1, 0, 0, 104.0, 1.9, 2787685, 103.5, 1.2, 98.0, 108.2, 63.4, 121.9),
        AngleDefinition("other", "P", "O3'", "C3'", 0, -1, -1, 119.7, 1.2, 2787486, 120.2, 1.5, 112.8, 128.1, 49.3, 138.0),
        AngleDefinition("other", "P", "O5'", "C5'", 0, 0, 0, 120.9, 1.6, 2812124, 120.6, 1.0, 115.1, 125.1, 109.2, 130.8),
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

    def _find_bond_definitions(self, res_name: str, altloc: str, atom1_name: str, atom2_name: str) -> List[BondDefinition]:
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
        raise NonStandardResidueException(f"Non-standard residue: {res_name}")

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
