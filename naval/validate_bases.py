import os
import sys
import numpy as np

from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_angle

BONDS = {
    "A": [
        ["N1", "C2", 1.339, 0.007, 544647, 1.338, 0.006, 1.307, 1.367, 1.278, 1.410],
        ["C2", "N3", 1.330, 0.007, 544637, 1.331, 0.006, 1.301, 1.355, 1.274, 1.411],
        ["N3", "C4", 1.346, 0.006, 544636, 1.343, 0.008, 1.299, 1.372, 1.264, 1.410],
        ["C4", "C5", 1.382, 0.008, 544646, 1.382, 0.008, 1.344, 1.414, 1.310, 1.466],
        ["C5", "C6", 1.406, 0.008, 544647, 1.404, 0.010, 1.351, 1.438, 1.302, 1.483],
        ["C6", "N1", 1.353, 0.007, 544647, 1.349, 0.007, 1.310, 1.378, 1.279, 1.414],
        ["C5", "N7", 1.388, 0.007, 544647, 1.385, 0.007, 1.348, 1.412, 1.316, 1.503],
        ["N7", "C8", 1.311, 0.007, 544631, 1.310, 0.005, 1.285, 1.333, 1.264, 1.380],
        ["C8", "N9", 1.370, 0.008, 544630, 1.371, 0.007, 1.334, 1.405, 1.285, 1.454],
        ["N9", "C4", 1.374, 0.007, 544645, 1.373, 0.009, 1.321, 1.407, 1.273, 1.454],
        ["C6", "N6", 1.334, 0.007, 544645, 1.334, 0.006, 1.301, 1.360, 1.270, 1.483],
    ],
    "G": [
        ["N1", "C2", 1.372, 0.006, 782271, 1.372, 0.007, 1.334, 1.400, 1.302, 1.433],
        ["C2", "N3", 1.327, 0.005, 782270, 1.322, 0.007, 1.288, 1.353, 1.260, 1.393],
        ["N3", "C4", 1.352, 0.006, 782270, 1.349, 0.007, 1.311, 1.377, 1.278, 1.402],
        ["C4", "C5", 1.379, 0.006, 782273, 1.377, 0.007, 1.342, 1.404, 1.312, 1.447],
        ["C5", "C6", 1.418, 0.008, 782271, 1.417, 0.009, 1.373, 1.451, 1.333, 1.488],
        ["C6", "N1", 1.392, 0.006, 782271, 1.390, 0.007, 1.353, 1.418, 1.314, 1.450],
        ["C5", "N7", 1.388, 0.006, 782274, 1.387, 0.006, 1.353, 1.411, 1.323, 1.447],
        ["N7", "C8", 1.308, 0.006, 782274, 1.304, 0.005, 1.280, 1.326, 1.257, 1.388],
        ["C8", "N9", 1.375, 0.006, 782272, 1.372, 0.007, 1.337, 1.399, 1.305, 1.444],
        ["N9", "C4", 1.374, 0.006, 782271, 1.374, 0.009, 1.333, 1.407, 1.287, 1.438],
        ["C6", "O6", 1.238, 0.007, 782242, 1.237, 0.007, 1.205, 1.271, 1.177, 1.321],
        ["C2", "N2", 1.338, 0.007, 782269, 1.339, 0.006, 1.308, 1.365, 1.280, 1.398],
    ],
    "U": [
        ["N1", "C2", 1.381, 0.009, 387309, 1.381, 0.011, 1.332, 1.426, 1.292, 1.478],
        ["C2", "N3", 1.373, 0.008, 387309, 1.372, 0.008, 1.332, 1.406, 1.299, 1.454],
        ["N3", "C4", 1.381, 0.008, 387310, 1.378, 0.009, 1.335, 1.412, 1.301, 1.455],
        ["C4", "C5", 1.432, 0.008, 387310, 1.430, 0.008, 1.394, 1.462, 1.362, 1.510],
        ["C5", "C6", 1.337, 0.008, 387310, 1.336, 0.006, 1.306, 1.360, 1.283, 1.473],
        ["C6", "N1", 1.374, 0.008, 387310, 1.374, 0.009, 1.330, 1.405, 1.290, 1.457],
        ["C2", "O2", 1.219, 0.008, 387309, 1.218, 0.007, 1.178, 1.250, 1.146, 1.287],
        ["C4", "O4", 1.231, 0.008, 387306, 1.231, 0.007, 1.195, 1.267, 1.172, 1.339],
    ],
    "T": [
        ["N1", "C2", 1.376, 0.008, 41088, 1.378, 0.007, 1.353, 1.421, 1.336, 1.488],
        ["C2", "N3", 1.372, 0.007, 41088, 1.371, 0.005, 1.339, 1.394, 1.305, 1.465],
        ["N3", "C4", 1.382, 0.008, 41088, 1.382, 0.005, 1.347, 1.405, 1.312, 1.437],
        ["C4", "C5", 1.446, 0.008, 41088, 1.445, 0.006, 1.418, 1.484, 1.402, 1.525],
        ["C5", "C6", 1.340, 0.007, 41088, 1.342, 0.005, 1.319, 1.379, 1.279, 1.463],
        ["C6", "N1", 1.381, 0.007, 41088, 1.381, 0.005, 1.353, 1.407, 1.315, 1.450],
        ["C2", "O2", 1.222, 0.008, 41086, 1.219, 0.005, 1.192, 1.248, 1.152, 1.274],
        ["C4", "O4", 1.229, 0.008, 41086, 1.228, 0.005, 1.203, 1.255, 1.163, 1.300],
        ["C7", "C5", 1.498, 0.006, 40843, 1.499, 0.004, 1.472, 1.524, 1.446, 1.559],
    ],
    "C": [
        ["N1", "C2", 1.395, 0.009, 596005, 1.398, 0.010, 1.354, 1.436, 1.313, 1.485],
        ["C2", "N3", 1.353, 0.007, 596005, 1.353, 0.007, 1.320, 1.381, 1.289, 1.417],
        ["N3", "C4", 1.337, 0.008, 596005, 1.332, 0.007, 1.295, 1.360, 1.261, 1.391],
        ["C4", "C5", 1.424, 0.010, 596004, 1.423, 0.007, 1.388, 1.449, 1.364, 1.477],
        ["C5", "C6", 1.338, 0.008, 596003, 1.338, 0.005, 1.311, 1.358, 1.291, 1.387],
        ["C6", "N1", 1.365, 0.007, 596005, 1.366, 0.007, 1.327, 1.392, 1.292, 1.418],
        ["C2", "O2", 1.240, 0.008, 596004, 1.240, 0.007, 1.205, 1.270, 1.177, 1.302],
        ["C4", "N4", 1.330, 0.008, 596001, 1.334, 0.006, 1.300, 1.357, 1.273, 1.388],
    ],
}
BONDS["DA"] = BONDS["A"]
BONDS["DG"] = BONDS["G"]
BONDS["DU"] = BONDS["U"]
BONDS["DT"] = BONDS["T"]
BONDS["DC"] = BONDS["C"]


ANGLES = {
    "A": [
        ["C6", "N1", "C2", 118.6, 0.6, 544647, 118.5, 0.7, 114.2, 121.8, 109.1, 126.0],
        ["N1", "C2", "N3", 129.4, 0.7, 544637, 129.3, 0.7, 126.4, 133.8, 118.2, 137.6],
        ["C2", "N3", "C4", 110.5, 0.6, 544636, 110.6, 0.9, 104.4, 114.1, 98.9, 121.0],
        ["N3", "C4", "C5", 126.9, 0.6, 544636, 126.8, 0.9, 123.1, 131.6, 117.6, 138.4],
        ["C4", "C5", "C6", 117.1, 0.5, 544646, 117.0, 0.7, 114.5, 120.4, 111.9, 124.9],
        ["C5", "C6", "N1", 117.5, 0.5, 544647, 117.7, 0.7, 113.8, 121.2, 109.4, 125.9],
        ["N3", "C4", "N9", 127.2, 0.7, 544635, 127.4, 0.9, 122.8, 131.1, 115.0, 136.4],
        ["C6", "C5", "N7", 132.2, 0.6, 544647, 132.1, 0.9, 126.6, 135.4, 119.3, 139.2],
        ["C5", "C4", "N9", 105.9, 0.4, 544645, 105.8, 0.6, 103.0, 108.6, 100.6, 111.5],
        ["C4", "N9", "C8", 105.7, 0.4, 544628, 105.8, 0.8, 102.0, 109.2, 98.9, 117.9],
        ["N9", "C8", "N7", 113.9, 0.5, 544630, 113.9, 0.8, 110.5, 118.5, 105.8, 122.5],
        ["C8", "N7", "C5", 103.8, 0.4, 544630, 103.7, 0.7, 98.6, 106.7, 93.5, 109.6],
        ["N7", "C5", "C4", 110.6, 0.5, 544646, 110.8, 0.7, 107.9, 115.0, 104.1, 119.3],
        ["N6", "C6", "N1", 118.6, 0.7, 544645, 118.7, 1.2, 113.5, 125.0, 108.2, 131.9],
        ["N6", "C6", "C5", 123.9, 0.7, 544645, 123.6, 1.1, 117.6, 128.5, 111.8, 133.5],
    ],
    "G": [
        ["C6", "N1", "C2", 125.4, 0.5, 782270, 125.0, 0.6, 121.5, 127.7, 117.2, 131.6],
        ["N1", "C2", "N3", 123.6, 0.5, 782270, 124.1, 0.7, 121.3, 127.8, 118.7, 132.1],
        ["C2", "N3", "C4", 112.0, 0.4, 782270, 111.9, 0.7, 108.0, 114.9, 102.2, 118.9],
        ["N3", "C4", "C5", 128.6, 0.5, 782270, 128.5, 0.8, 125.1, 132.1, 120.7, 139.1],
        ["C4", "C5", "C6", 118.9, 0.4, 782270, 119.0, 0.7, 116.3, 122.3, 112.9, 125.7],
        ["C5", "C6", "N1", 111.5, 0.5, 782271, 111.5, 0.7, 107.9, 114.9, 102.7, 118.7],
        ["N3", "C4", "N9", 125.8, 0.7, 782268, 126.1, 1.0, 121.9, 130.1, 112.5, 133.6],
        ["C6", "C5", "N7", 130.3, 0.5, 782271, 130.2, 0.9, 125.7, 133.5, 120.8, 136.3],
        ["C5", "C4", "N9", 105.6, 0.5, 782271, 105.4, 0.6, 102.6, 108.1, 100.5, 110.8],
        ["C4", "N9", "C8", 106.2, 0.4, 782271, 106.3, 0.7, 103.0, 109.7, 99.7, 112.5],
        ["N9", "C8", "N7", 113.2, 0.4, 782272, 113.2, 0.7, 110.0, 116.7, 107.5, 122.4],
        ["C8", "N7", "C5", 104.2, 0.4, 782274, 104.3, 0.6, 100.6, 107.1, 93.1, 110.3],
        ["N7", "C5", "C4", 110.8, 0.4, 782273, 110.8, 0.6, 108.1, 114.0, 104.9, 118.4],
        ["O6", "C6", "N1", 120.1, 0.5, 782242, 120.0, 1.1, 114.9, 125.8, 110.2, 132.4],
        ["O6", "C6", "C5", 128.4, 0.6, 782242, 128.5, 1.0, 123.2, 133.1, 118.0, 138.6],
        ["N2", "C2", "N1", 116.5, 0.6, 782269, 116.2, 1.1, 110.6, 120.9, 102.4, 125.2],
        ["N2", "C2", "N3", 119.9, 0.6, 782269, 119.8, 1.0, 114.8, 124.4, 108.9, 132.2],
    ],
    "U": [
        ["C6", "N1", "C2", 121.1, 0.5, 387309, 120.8, 1.0, 116.4, 125.1, 111.8, 128.5],
        ["N1", "C2", "N3", 114.9, 0.6, 387309, 115.2, 0.9, 111.7, 120.0, 108.0, 125.0],
        ["C2", "N3", "C4", 127.0, 0.5, 387309, 126.9, 0.8, 122.5, 130.1, 118.3, 134.4],
        ["N3", "C4", "C5", 114.5, 0.6, 387310, 114.6, 0.8, 110.2, 118.5, 103.9, 122.1],
        ["C4", "C5", "C6", 119.7, 0.6, 387310, 119.8, 0.8, 116.7, 124.3, 113.1, 128.7],
        ["C5", "C6", "N1", 122.7, 0.5, 387310, 122.7, 0.9, 118.1, 126.4, 113.3, 130.4],
        ["O2", "C2", "N1", 122.8, 0.7, 387309, 123.0, 1.3, 117.3, 128.5, 111.8, 132.3],
        ["O2", "C2", "N3", 122.3, 0.6, 387309, 121.8, 1.3, 115.4, 127.0, 110.7, 131.0],
        ["O4", "C4", "C5", 126.0, 0.7, 387305, 126.0, 1.0, 121.3, 131.2, 117.3, 136.3],
        ["O4", "C4", "N3", 119.5, 0.7, 387305, 119.5, 1.1, 114.1, 124.7, 109.3, 130.1],
    ],
    "T": [
        ["C6", "N1", "C2", 121.2, 0.5, 41088, 121.2, 0.3, 119.0, 123.0, 116.8, 124.7],
        ["N1", "C2", "N3", 114.7, 0.6, 41088, 114.7, 0.4, 113.0, 119.9, 111.7, 123.5],
        ["C2", "N3", "C4", 127.1, 0.5, 41088, 127.2, 0.4, 121.9, 129.0, 119.0, 130.7],
        ["N3", "C4", "C5", 115.2, 0.5, 41088, 115.3, 0.3, 113.7, 118.7, 111.1, 122.1],
        ["C4", "C5", "C6", 118.1, 0.5, 41088, 118.1, 0.3, 116.5, 120.3, 113.7, 124.0],
        ["C5", "C6", "N1", 123.6, 0.5, 41088, 123.5, 0.4, 120.4, 125.1, 117.3, 127.8],
        ["O2", "C2", "N1", 123.0, 0.7, 41086, 123.3, 0.5, 120.1, 126.7, 116.3, 132.0],
        ["O2", "C2", "N3", 122.3, 0.6, 41086, 122.0, 0.6, 117.5, 124.2, 113.2, 126.7],
        ["O4", "C4", "C5", 125.0, 0.7, 41086, 123.3, 1.3, 119.9, 127.3, 117.4, 130.0],
        ["O4", "C4", "N3", 119.8, 0.6, 41086, 121.4, 1.2, 117.2, 124.3, 115.0, 126.9],
        ["C7", "C5", "C4", 118.7, 0.6, 40843, 119.9, 0.7, 116.9, 123.1, 114.1, 125.4],
        ["C7", "C5", "C6", 123.2, 0.6, 40843, 122.1, 0.7, 118.5, 124.5, 114.9, 126.9],
    ],
    "C": [
        ["C6", "N1", "C2", 120.3, 0.5, 596004, 120.1, 0.9, 116.5, 124.3, 112.4, 128.4],
        ["N1", "C2", "N3", 119.1, 0.6, 596005, 119.3, 0.8, 116.0, 123.2, 112.4, 127.4],
        ["C2", "N3", "C4", 120.1, 0.5, 596005, 119.9, 0.7, 116.1, 122.8, 112.4, 127.8],
        ["N3", "C4", "C5", 121.6, 0.6, 596004, 121.9, 0.6, 118.9, 125.1, 114.3, 128.0],
        ["C4", "C5", "C6", 117.5, 0.5, 596003, 117.5, 0.6, 114.8, 120.8, 111.9, 123.9],
        ["C5", "C6", "N1", 121.2, 0.6, 596003, 121.2, 0.9, 117.0, 124.7, 113.5, 127.4],
        ["O2", "C2", "N1", 118.8, 0.8, 596004, 119.1, 1.2, 113.8, 124.2, 108.9, 131.2],
        ["O2", "C2", "N3", 122.0, 0.6, 596004, 121.6, 1.1, 116.3, 126.5, 109.3, 130.8],
        ["N4", "C4", "C5", 120.3, 0.7, 596000, 120.1, 0.9, 115.4, 124.4, 111.4, 129.4],
        ["N4", "C4", "N3", 118.1, 0.6, 596001, 118.0, 1.0, 113.4, 122.6, 109.2, 126.6],
    ],
}
ANGLES["DA"] = ANGLES["A"]
ANGLES["DG"] = ANGLES["G"]
ANGLES["DU"] = ANGLES["U"]
ANGLES["DT"] = ANGLES["T"]
ANGLES["DC"] = ANGLES["C"]


def read_structure(pdb_file_path):
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdbcode = os.path.basename(pdb_file_path)[0:4]
    return parser.get_structure(pdbcode, pdb_file_path)


def validate_bases(structure, csd_sig=3):
    pdbcode = structure.id

    res_names = list(BONDS.keys())

    bond_preferred = 0
    bond_allowed = 0
    bond_suspicious = 0
    bond_outlier = 0

    angle_preferred = 0
    angle_allowed = 0
    angle_suspicious = 0
    angle_outlier = 0

    for model in structure:
        model_id = model.get_id()
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_name = residue.get_resname()
                if res_name in res_names:
                    res_full_id = residue.get_id()

                    resseq = res_full_id[1]
                    inscode = res_full_id[2]

                    for definitions in BONDS[res_name]:
                        atom_group1 = chain[resseq][definitions[0]]
                        atom_group2 = chain[resseq][definitions[1]]

                        csd_3low = definitions[2] - csd_sig * definitions[3]
                        csd_3high = definitions[2] + csd_sig * definitions[3]

                        pdb_3low = definitions[7]
                        pdb_3high = definitions[8]

                        pdb_4low = definitions[9]
                        pdb_4high = definitions[10]

                        atoms1 = atom_group1.disordered_get_list() if atom_group1.is_disordered() else [atom_group1]
                        atoms2 = atom_group2.disordered_get_list() if atom_group2.is_disordered() else [atom_group2]

                        for atom1 in atoms1:
                            for atom2 in atoms2:
                                if (
                                    atom1.get_altloc() == atom2.get_altloc()
                                    or atom1.get_altloc() == " "
                                    or atom2.get_altloc() == " "
                                ):
                                    dist = round(atom2 - atom1, 4)

                                    in_csd = csd_3low <= dist and dist <= csd_3high
                                    in_3pdb = False if in_csd else (pdb_3low <= dist and dist <= pdb_3high)
                                    in_4pdb = False if (in_csd or in_3pdb) else (pdb_4low <= dist and dist <= pdb_4high)
                                    out = False if in_4pdb or in_3pdb or in_csd else True

                                    line = ",".join(
                                        str(_)
                                        for _ in (
                                            "bond",
                                            pdbcode,
                                            model_id,
                                            chain_id,
                                            res_name,
                                            resseq,
                                            inscode,
                                            atom1.get_name(),
                                            atom1.get_altloc(),
                                            atom2.get_name(),
                                            atom2.get_altloc(),
                                            dist,
                                            in_csd,
                                            in_3pdb,
                                            in_4pdb,
                                            out,
                                        )
                                    )
                                    print(line)
                                    bond_preferred += int(in_csd)
                                    bond_allowed += int(in_3pdb)
                                    bond_suspicious += int(in_4pdb)
                                    bond_outlier += int(out)

                    for definitions in ANGLES[res_name]:
                        atom_group1 = chain[resseq][definitions[0]]
                        atom_group2 = chain[resseq][definitions[1]]
                        atom_group3 = chain[resseq][definitions[2]]

                        csd_3low = definitions[3] - csd_sig * definitions[4]
                        csd_3high = definitions[3] + csd_sig * definitions[4]

                        pdb_3low = definitions[8]
                        pdb_3high = definitions[9]

                        pdb_4low = definitions[10]
                        pdb_4high = definitions[11]

                        atoms1 = atom_group1.disordered_get_list() if atom_group1.is_disordered() else [atom_group1]
                        atoms2 = atom_group2.disordered_get_list() if atom_group2.is_disordered() else [atom_group2]
                        atoms3 = atom_group3.disordered_get_list() if atom_group3.is_disordered() else [atom_group3]

                        for atom1 in atoms1:
                            for atom2 in atoms2:
                                for atom3 in atoms3:
                                    if (
                                        (
                                            atom1.get_altloc() == atom2.get_altloc()
                                            and atom2.get_altloc() == atom3.get_altloc()
                                        )
                                        or (atom1.get_altloc() == " " and atom2.get_altloc() == atom3.get_altloc())
                                        or (atom2.get_altloc() == " " and atom1.get_altloc() == atom3.get_altloc())
                                        or (atom3.get_altloc() == " " and atom1.get_altloc() == atom2.get_altloc())
                                        or (
                                            atom1.get_altloc() != " "
                                            and atom2.get_altloc() == " "
                                            and atom3.get_altloc() == " "
                                        )
                                        or (
                                            atom1.get_altloc() == " "
                                            and atom2.get_altloc() != " "
                                            and atom3.get_altloc() == " "
                                        )
                                        or (
                                            atom1.get_altloc() == " "
                                            and atom2.get_altloc() == " "
                                            and atom3.get_altloc() != " "
                                        )
                                    ):
                                        angle_value = calc_angle(
                                            atom1.get_vector(),
                                            atom2.get_vector(),
                                            atom3.get_vector(),
                                        )
                                        angle_value = round(np.rad2deg(angle_value), 2)

                                        in_csd = csd_3low <= angle_value and angle_value <= csd_3high
                                        in_3pdb = (
                                            False if in_csd else (pdb_3low <= angle_value and angle_value <= pdb_3high)
                                        )
                                        in_4pdb = (
                                            False
                                            if (in_csd or in_3pdb)
                                            else (pdb_4low <= angle_value and angle_value <= pdb_4high)
                                        )
                                        out = False if in_4pdb or in_3pdb or in_csd else True

                                        line = ",".join(
                                            str(_)
                                            for _ in (
                                                "angle",
                                                pdbcode,
                                                model_id,
                                                chain_id,
                                                res_name,
                                                resseq,
                                                inscode,
                                                atom1.get_name(),
                                                atom1.get_altloc(),
                                                atom2.get_name(),
                                                atom2.get_altloc(),
                                                atom3.get_name(),
                                                atom3.get_altloc(),
                                                angle_value,
                                                in_csd,
                                                in_3pdb,
                                                in_4pdb,
                                                out,
                                            )
                                        )
                                        print(line)
                                        angle_preferred += int(in_csd)
                                        angle_allowed += int(in_3pdb)
                                        angle_suspicious += int(in_4pdb)
                                        angle_outlier += int(out)
    bond_total = bond_preferred + bond_allowed + bond_suspicious + bond_outlier
    print(f"{pdbcode} Bases bond validation summary")
    print(f"\tPreferred: {100*bond_preferred/bond_total:.2f}%")
    print(f"\tAllowed: {100*bond_allowed/bond_total:.2f}%")
    print(f"\tSuspicious: {100*bond_suspicious/bond_total:.2f}%")
    print(f"\tOutlier: {100*bond_outlier/bond_total:.2f}%")

    angle_total = angle_preferred + angle_allowed + angle_suspicious + angle_outlier
    print(f"{pdbcode} Bases angle validation summary")
    print(f"\tPreferred: {100*angle_preferred/angle_total:.2f}%")
    print(f"\tAllowed: {100*angle_allowed/angle_total:.2f}%")
    print(f"\tSuspicious: {100*angle_suspicious/angle_total:.2f}%")
    print(f"\tOutlier: {100*angle_outlier/angle_total:.2f}%")


def main(pdb_filepath):
    structure = read_structure(pdb_filepath)
    validate_bases(structure)
    return 0


if __name__ == "__main__":
    pdb_filepath = sys.argv[1]
    main(pdb_filepath)
