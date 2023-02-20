from naval.nucleotide_geometry import NucleotideGeometry
from naval.restraint_definition import AngleDefinition, BondDefinition
from naval.validators.validator import Validator

BASES_BONDS = {
    "A": [
        BondDefinition("A/DA", "N1", "C2", 0, 0, 1.339, 0.007, 651941, 1.338, 0.006, 1.308, 1.367, 1.279, 1.406),
        BondDefinition("A/DA", "C2", "N3", 0, 0, 1.330, 0.007, 651916, 1.331, 0.006, 1.302, 1.355, 1.276, 1.405),
        BondDefinition("A/DA", "N3", "C4", 0, 0, 1.346, 0.006, 651915, 1.343, 0.007, 1.301, 1.371, 1.266, 1.406),
        BondDefinition("A/DA", "C4", "C5", 0, 0, 1.382, 0.008, 651938, 1.382, 0.007, 1.345, 1.413, 1.312, 1.457),
        BondDefinition("A/DA", "C5", "C6", 0, 0, 1.406, 0.008, 651941, 1.404, 0.010, 1.352, 1.437, 1.288, 1.476),
        BondDefinition("A/DA", "C6", "N1", 0, 0, 1.353, 0.007, 651940, 1.349, 0.007, 1.311, 1.377, 1.277, 1.409),
        BondDefinition("A/DA", "C5", "N7", 0, 0, 1.388, 0.007, 651941, 1.385, 0.007, 1.349, 1.411, 1.317, 1.501),
        BondDefinition("A/DA", "N7", "C8", 0, 0, 1.311, 0.007, 651916, 1.310, 0.005, 1.286, 1.332, 1.267, 1.378),
        BondDefinition("A/DA", "C8", "N9", 0, 0, 1.370, 0.008, 651915, 1.371, 0.007, 1.335, 1.403, 1.297, 1.451),
        BondDefinition("A/DA", "N9", "C4", 0, 0, 1.374, 0.007, 651939, 1.373, 0.009, 1.323, 1.406, 1.276, 1.448),
        BondDefinition("A/DA", "C6", "N6", 0, 0, 1.334, 0.007, 651939, 1.334, 0.006, 1.301, 1.360, 1.272, 1.483),
    ],
    "G": [
        BondDefinition("G/DG", "N1", "C2", 0, 0, 1.372, 0.006, 943090, 1.372, 0.007, 1.335, 1.399, 1.304, 1.433),
        BondDefinition("G/DG", "C2", "N3", 0, 0, 1.327, 0.005, 943089, 1.322, 0.007, 1.288, 1.352, 1.260, 1.389),
        BondDefinition("G/DG", "N3", "C4", 0, 0, 1.352, 0.006, 943089, 1.349, 0.007, 1.312, 1.376, 1.280, 1.402),
        BondDefinition("G/DG", "C4", "C5", 0, 0, 1.379, 0.006, 943092, 1.377, 0.007, 1.343, 1.403, 1.316, 1.444),
        BondDefinition("G/DG", "C5", "C6", 0, 0, 1.418, 0.008, 943088, 1.417, 0.009, 1.374, 1.450, 1.334, 1.486),
        BondDefinition("G/DG", "C6", "N1", 0, 0, 1.392, 0.006, 943088, 1.390, 0.007, 1.354, 1.418, 1.316, 1.449),
        BondDefinition("G/DG", "C5", "N7", 0, 0, 1.388, 0.006, 943093, 1.387, 0.006, 1.355, 1.411, 1.325, 1.445),
        BondDefinition("G/DG", "N7", "C8", 0, 0, 1.308, 0.006, 943093, 1.304, 0.005, 1.280, 1.326, 1.258, 1.385),
        BondDefinition("G/DG", "C8", "N9", 0, 0, 1.375, 0.006, 943091, 1.372, 0.007, 1.338, 1.398, 1.306, 1.443),
        BondDefinition("G/DG", "N9", "C4", 0, 0, 1.374, 0.006, 943090, 1.374, 0.008, 1.335, 1.407, 1.288, 1.437),
        BondDefinition("G/DG", "C6", "O6", 0, 0, 1.238, 0.007, 943053, 1.237, 0.007, 1.205, 1.270, 1.178, 1.322),
        BondDefinition("G/DG", "C2", "N2", 0, 0, 1.338, 0.007, 943088, 1.339, 0.006, 1.309, 1.364, 1.281, 1.395),
    ],
    "U": [
        BondDefinition("U/DU", "N1", "C2", 0, 0, 1.381, 0.009, 459255, 1.381, 0.010, 1.334, 1.424, 1.293, 1.474),
        BondDefinition("U/DU", "C2", "N3", 0, 0, 1.373, 0.008, 459255, 1.372, 0.008, 1.333, 1.405, 1.300, 1.455),
        BondDefinition("U/DU", "N3", "C4", 0, 0, 1.381, 0.008, 459255, 1.378, 0.008, 1.336, 1.411, 1.302, 1.454),
        BondDefinition("U/DU", "C4", "C5", 0, 0, 1.432, 0.008, 459255, 1.431, 0.008, 1.394, 1.461, 1.365, 1.506),
        BondDefinition("U/DU", "C5", "C6", 0, 0, 1.337, 0.008, 459253, 1.336, 0.006, 1.307, 1.359, 1.282, 1.490),
        BondDefinition("U/DU", "C6", "N1", 0, 0, 1.374, 0.008, 459253, 1.374, 0.008, 1.332, 1.404, 1.292, 1.444),
        BondDefinition("U/DU", "C2", "O2", 0, 0, 1.219, 0.008, 459256, 1.218, 0.007, 1.179, 1.249, 1.147, 1.282),
        BondDefinition("U/DU", "C4", "O4", 0, 0, 1.231, 0.008, 459252, 1.231, 0.007, 1.196, 1.266, 1.172, 1.337),
    ],
    "T": [
        BondDefinition("T/DT", "N1", "C2", 0, 0, 1.376, 0.008, 50206, 1.378, 0.007, 1.351, 1.420, 1.333, 1.485),
        BondDefinition("T/DT", "C2", "N3", 0, 0, 1.372, 0.007, 50206, 1.371, 0.005, 1.340, 1.393, 1.307, 1.463),
        BondDefinition("T/DT", "N3", "C4", 0, 0, 1.382, 0.008, 50206, 1.382, 0.005, 1.348, 1.405, 1.314, 1.434),
        BondDefinition("T/DT", "C4", "C5", 0, 0, 1.446, 0.008, 50206, 1.445, 0.006, 1.417, 1.482, 1.394, 1.525),
        BondDefinition("T/DT", "C5", "C6", 0, 0, 1.340, 0.007, 50206, 1.342, 0.005, 1.319, 1.374, 1.278, 1.462),
        BondDefinition("T/DT", "C6", "N1", 0, 0, 1.381, 0.007, 50206, 1.381, 0.005, 1.353, 1.407, 1.289, 1.449),
        BondDefinition("T/DT", "C2", "O2", 0, 0, 1.222, 0.008, 50202, 1.219, 0.005, 1.192, 1.248, 1.158, 1.274),
        BondDefinition("T/DT", "C4", "O4", 0, 0, 1.229, 0.008, 50204, 1.228, 0.005, 1.203, 1.254, 1.170, 1.298),
        BondDefinition("T/DT", "C7", "C5", 0, 0, 1.498, 0.006, 49950, 1.499, 0.005, 1.472, 1.525, 1.447, 1.586),
    ],
    "C": [
        BondDefinition("C/DC", "N1", "C2", 0, 0, 1.395, 0.009, 718608, 1.398, 0.009, 1.355, 1.435, 1.316, 1.483),
        BondDefinition("C/DC", "C2", "N3", 0, 0, 1.353, 0.007, 718608, 1.353, 0.007, 1.321, 1.381, 1.290, 1.413),
        BondDefinition("C/DC", "N3", "C4", 0, 0, 1.337, 0.008, 718607, 1.333, 0.007, 1.296, 1.359, 1.263, 1.392),
        BondDefinition("C/DC", "C4", "C5", 0, 0, 1.424, 0.010, 718606, 1.423, 0.007, 1.389, 1.449, 1.364, 1.480),
        BondDefinition("C/DC", "C5", "C6", 0, 0, 1.338, 0.008, 718605, 1.338, 0.005, 1.311, 1.357, 1.291, 1.386),
        BondDefinition("C/DC", "C6", "N1", 0, 0, 1.365, 0.007, 718607, 1.366, 0.007, 1.328, 1.391, 1.294, 1.418),
        BondDefinition("C/DC", "C2", "O2", 0, 0, 1.240, 0.008, 718606, 1.239, 0.006, 1.206, 1.269, 1.179, 1.297),
        BondDefinition("C/DC", "C4", "N4", 0, 0, 1.330, 0.008, 718603, 1.334, 0.006, 1.301, 1.357, 1.270, 1.384),
    ],
}
BASES_BONDS["DA"] = BASES_BONDS["A"]
BASES_BONDS["DG"] = BASES_BONDS["G"]
BASES_BONDS["DU"] = BASES_BONDS["U"]
BASES_BONDS["DT"] = BASES_BONDS["T"]
BASES_BONDS["DC"] = BASES_BONDS["C"]


BASES_ANGLES = {
    "A": [
        AngleDefinition("A/DA", "C6", "N1", "C2", 0, 0, 0, 118.6, 0.6, 651941, 118.5, 0.7, 114.4, 121.7, 109.5, 125.8),
        AngleDefinition("A/DA", "N1", "C2", "N3", 0, 0, 0, 129.4, 0.7, 651916, 129.3, 0.7, 126.5, 133.6, 118.2, 137.4),
        AngleDefinition("A/DA", "C2", "N3", "C4", 0, 0, 0, 110.5, 0.6, 651915, 110.6, 0.8, 104.7, 113.9, 99.1, 121.0),
        AngleDefinition("A/DA", "N3", "C4", "C5", 0, 0, 0, 126.9, 0.6, 651915, 126.8, 0.8, 123.2, 131.4, 117.8, 138.2),
        AngleDefinition("A/DA", "C4", "C5", "C6", 0, 0, 0, 117.1, 0.5, 651940, 117.0, 0.6, 114.5, 120.3, 112.1, 124.2),
        AngleDefinition("A/DA", "C5", "C6", "N1", 0, 0, 0, 117.5, 0.5, 651941, 117.7, 0.7, 113.9, 121.0, 109.5, 125.6),
        AngleDefinition("A/DA", "N3", "C4", "N9", 0, 0, 0, 127.2, 0.7, 651914, 127.4, 0.8, 123.0, 131.0, 115.2, 136.2),
        AngleDefinition("A/DA", "C6", "C5", "N7", 0, 0, 0, 132.2, 0.6, 651941, 132.2, 0.9, 126.7, 135.3, 119.4, 139.0),
        AngleDefinition("A/DA", "C5", "C4", "N9", 0, 0, 0, 105.9, 0.4, 651939, 105.8, 0.6, 103.1, 108.5, 100.6, 111.3),
        AngleDefinition("A/DA", "C4", "N9", "C8", 0, 0, 0, 105.7, 0.4, 651913, 105.8, 0.7, 102.1, 109.1, 98.9, 115.2),
        AngleDefinition("A/DA", "N9", "C8", "N7", 0, 0, 0, 113.9, 0.5, 651915, 113.9, 0.8, 110.6, 118.4, 105.9, 122.4),
        AngleDefinition("A/DA", "C8", "N7", "C5", 0, 0, 0, 103.8, 0.4, 651915, 103.8, 0.7, 98.8, 106.6, 93.8, 109.5),
        AngleDefinition("A/DA", "N7", "C5", "C4", 0, 0, 0, 110.6, 0.5, 651940, 110.8, 0.7, 108.0, 114.9, 104.2, 119.2),
        AngleDefinition("A/DA", "N6", "C6", "N1", 0, 0, 0, 118.6, 0.7, 651939, 118.7, 1.1, 113.6, 124.8, 108.5, 131.8),
        AngleDefinition("A/DA", "N6", "C6", "C5", 0, 0, 0, 123.9, 0.7, 651939, 123.6, 1.1, 117.8, 128.4, 112.0, 133.1),
    ],
    "G": [
        AngleDefinition("G/DG", "C6", "N1", "C2", 0, 0, 0, 125.4, 0.5, 943087, 125.0, 0.6, 121.6, 127.7, 117.4, 131.3),
        AngleDefinition("G/DG", "N1", "C2", "N3", 0, 0, 0, 123.6, 0.5, 943089, 124.1, 0.6, 121.4, 127.7, 118.9, 131.9),
        AngleDefinition("G/DG", "C2", "N3", "C4", 0, 0, 0, 112.0, 0.4, 943089, 111.9, 0.7, 108.2, 114.8, 102.5, 118.9),
        AngleDefinition("G/DG", "N3", "C4", "C5", 0, 0, 0, 128.6, 0.5, 943089, 128.5, 0.8, 125.1, 132.0, 120.9, 138.9),
        AngleDefinition("G/DG", "C4", "C5", "C6", 0, 0, 0, 118.9, 0.4, 943087, 119.0, 0.6, 116.4, 122.2, 113.4, 125.7),
        AngleDefinition("G/DG", "C5", "C6", "N1", 0, 0, 0, 111.5, 0.5, 943088, 111.5, 0.7, 108.0, 114.8, 102.9, 118.8),
        AngleDefinition("G/DG", "N3", "C4", "N9", 0, 0, 0, 125.8, 0.7, 943087, 126.1, 1.0, 122.0, 130.0, 112.6, 133.5),
        AngleDefinition("G/DG", "C6", "C5", "N7", 0, 0, 0, 130.3, 0.5, 943088, 130.2, 0.9, 125.8, 133.5, 120.8, 136.2),
        AngleDefinition("G/DG", "C5", "C4", "N9", 0, 0, 0, 105.6, 0.5, 943090, 105.4, 0.6, 102.7, 108.1, 100.6, 110.7),
        AngleDefinition("G/DG", "C4", "N9", "C8", 0, 0, 0, 106.2, 0.4, 943090, 106.3, 0.7, 103.1, 109.6, 99.8, 112.4),
        AngleDefinition("G/DG", "N9", "C8", "N7", 0, 0, 0, 113.2, 0.4, 943091, 113.2, 0.7, 110.0, 116.6, 107.5, 121.9),
        AngleDefinition("G/DG", "C8", "N7", "C5", 0, 0, 0, 104.2, 0.4, 943093, 104.3, 0.6, 100.8, 107.0, 93.5, 110.1),
        AngleDefinition("G/DG", "N7", "C5", "C4", 0, 0, 0, 110.8, 0.4, 943092, 110.8, 0.6, 108.2, 113.9, 105.2, 118.2),
        AngleDefinition("G/DG", "O6", "C6", "N1", 0, 0, 0, 120.1, 0.5, 943053, 120.0, 1.1, 115.1, 125.6, 110.3, 132.1),
        AngleDefinition("G/DG", "O6", "C6", "C5", 0, 0, 0, 128.4, 0.6, 943053, 128.5, 1.0, 123.4, 132.9, 118.0, 138.4),
        AngleDefinition("G/DG", "N2", "C2", "N1", 0, 0, 0, 116.5, 0.6, 943088, 116.2, 1.1, 110.7, 120.8, 102.4, 125.0),
        AngleDefinition("G/DG", "N2", "C2", "N3", 0, 0, 0, 119.9, 0.6, 943088, 119.8, 1.0, 114.9, 124.2, 109.1, 132.2),
    ],
    "U": [
        AngleDefinition("U/DU", "C6", "N1", "C2", 0, 0, 0, 121.1, 0.5, 459254, 120.8, 1.0, 116.5, 125.0, 112.0, 128.2),
        AngleDefinition("U/DU", "N1", "C2", "N3", 0, 0, 0, 114.9, 0.6, 459256, 115.2, 0.8, 111.8, 119.8, 108.1, 124.9),
        AngleDefinition("U/DU", "C2", "N3", "C4", 0, 0, 0, 127.0, 0.5, 459256, 126.9, 0.8, 122.6, 130.0, 118.5, 134.3),
        AngleDefinition("U/DU", "N3", "C4", "C5", 0, 0, 0, 114.5, 0.6, 459256, 114.6, 0.8, 110.3, 118.4, 104.0, 121.9),
        AngleDefinition("U/DU", "C4", "C5", "C6", 0, 0, 0, 119.7, 0.6, 459254, 119.8, 0.7, 116.8, 124.1, 113.2, 128.5),
        AngleDefinition("U/DU", "C5", "C6", "N1", 0, 0, 0, 122.7, 0.5, 459254, 122.7, 0.9, 118.2, 126.4, 113.5, 130.2),
        AngleDefinition("U/DU", "O2", "C2", "N1", 0, 0, 0, 122.8, 0.7, 459256, 123.0, 1.2, 117.4, 128.4, 112.1, 132.1),
        AngleDefinition("U/DU", "O2", "C2", "N3", 0, 0, 0, 122.3, 0.6, 459256, 121.8, 1.2, 115.6, 126.9, 110.9, 130.8),
        AngleDefinition("U/DU", "O4", "C4", "C5", 0, 0, 0, 126.0, 0.7, 459251, 126.0, 1.0, 121.4, 131.0, 117.5, 135.4),
        AngleDefinition("U/DU", "O4", "C4", "N3", 0, 0, 0, 119.5, 0.7, 459251, 119.5, 1.0, 114.2, 124.6, 109.5, 130.0),
    ],
    "T": [
        AngleDefinition("T/DT", "C6", "N1", "C2", 0, 0, 0, 121.2, 0.5, 50206, 121.2, 0.4, 119.0, 123.0, 103.8, 124.8),
        AngleDefinition("T/DT", "N1", "C2", "N3", 0, 0, 0, 114.7, 0.6, 50206, 114.7, 0.4, 113.0, 119.2, 111.6, 123.9),
        AngleDefinition("T/DT", "C2", "N3", "C4", 0, 0, 0, 127.1, 0.5, 50206, 127.2, 0.4, 122.8, 129.0, 119.0, 130.8),
        AngleDefinition("T/DT", "N3", "C4", "C5", 0, 0, 0, 115.2, 0.5, 50206, 115.3, 0.3, 113.6, 118.4, 111.3, 122.0),
        AngleDefinition("T/DT", "C4", "C5", "C6", 0, 0, 0, 118.1, 0.5, 50206, 118.1, 0.3, 116.4, 120.3, 105.4, 123.9),
        AngleDefinition("T/DT", "C5", "C6", "N1", 0, 0, 0, 123.6, 0.5, 50206, 123.5, 0.4, 120.5, 125.2, 117.6, 146.8),
        AngleDefinition("T/DT", "O2", "C2", "N1", 0, 0, 0, 123.0, 0.7, 50202, 123.3, 0.6, 120.2, 126.8, 116.3, 131.9),
        AngleDefinition("T/DT", "O2", "C2", "N3", 0, 0, 0, 122.3, 0.6, 50202, 122.0, 0.6, 117.5, 124.2, 112.2, 126.8),
        AngleDefinition("T/DT", "O4", "C4", "C5", 0, 0, 0, 125.0, 0.7, 50204, 123.3, 1.3, 119.9, 127.2, 117.3, 130.8),
        AngleDefinition("T/DT", "O4", "C4", "N3", 0, 0, 0, 119.8, 0.6, 50204, 121.4, 1.2, 117.5, 124.3, 114.4, 126.3),
        AngleDefinition("T/DT", "C7", "C5", "C4", 0, 0, 0, 118.7, 0.6, 49950, 119.9, 0.7, 116.9, 123.2, 114.2, 125.5),
        AngleDefinition("T/DT", "C7", "C5", "C6", 0, 0, 0, 123.2, 0.6, 49950, 122.1, 0.7, 118.6, 124.6, 115.1, 127.7),
    ],
    "C": [
        AngleDefinition("C/DC", "C6", "N1", "C2", 0, 0, 0, 120.3, 0.5, 718607, 120.1, 0.9, 116.6, 124.2, 112.7, 128.2),
        AngleDefinition("C/DC", "N1", "C2", "N3", 0, 0, 0, 119.1, 0.6, 718608, 119.3, 0.7, 116.1, 123.1, 112.5, 127.2),
        AngleDefinition("C/DC", "C2", "N3", "C4", 0, 0, 0, 120.1, 0.5, 718607, 119.9, 0.7, 116.2, 122.7, 112.8, 127.4),
        AngleDefinition("C/DC", "N3", "C4", "C5", 0, 0, 0, 121.6, 0.6, 718606, 121.9, 0.6, 119.0, 125.0, 114.6, 127.9),
        AngleDefinition("C/DC", "C4", "C5", "C6", 0, 0, 0, 117.5, 0.5, 718605, 117.5, 0.6, 114.9, 120.7, 112.1, 123.8),
        AngleDefinition("C/DC", "C5", "C6", "N1", 0, 0, 0, 121.2, 0.6, 718605, 121.1, 0.9, 117.1, 124.7, 113.7, 127.5),
        AngleDefinition("C/DC", "O2", "C2", "N1", 0, 0, 0, 118.8, 0.8, 718606, 119.1, 1.1, 113.9, 124.1, 109.1, 130.9),
        AngleDefinition("C/DC", "O2", "C2", "N3", 0, 0, 0, 122.0, 0.6, 718606, 121.6, 1.1, 116.4, 126.4, 110.0, 130.5),
        AngleDefinition("C/DC", "N4", "C4", "C5", 0, 0, 0, 120.3, 0.7, 718602, 120.1, 0.9, 115.5, 124.3, 111.7, 129.2),
        AngleDefinition("C/DC", "N4", "C4", "N3", 0, 0, 0, 118.1, 0.6, 718603, 118.0, 0.9, 113.5, 122.5, 109.4, 126.5),
    ],
}
BASES_ANGLES["DA"] = BASES_ANGLES["A"]
BASES_ANGLES["DG"] = BASES_ANGLES["G"]
BASES_ANGLES["DU"] = BASES_ANGLES["U"]
BASES_ANGLES["DT"] = BASES_ANGLES["T"]
BASES_ANGLES["DC"] = BASES_ANGLES["C"]


class BasesValidator(Validator):
    """
    Validator for nucleotde basees
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        super().__init__(geometry, csd_sig)

        self.bonds_definition = BASES_BONDS
        self.angles_definition = BASES_ANGLES
