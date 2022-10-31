import math
from typing import Dict, List, Optional

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_dihedral

from naval.nucleotide_definitions import PURINES_RES_NAMES
from naval.residue_cache_entry import ResidueCacheEntry


class NucleotideGeometry:
    """
    Class to keep cache torsion angles for given residue
    """

    # pylint: disable=too-many-public-methods
    # pylint: disable=too-many-instance-attributes
    def __init__(self, residue_entry: ResidueCacheEntry) -> None:
        """Simple container class to keep torsion angles of the residue.
        Calculates torison anles for all alternative conformations.
        Maps torsion angles to string classes.
        """
        self.residue_entry = residue_entry

        self.alpha: Dict[str, Optional[float]] = {}  # O3'(i-1)-P-O5'-C5'
        self.beta: Dict[str, Optional[float]] = {}  # P-O5'-C5'-C4'
        self.gamma: Dict[str, Optional[float]] = {}  # O5'-C5'-C4'-C3'
        self.delta: Dict[str, Optional[float]] = {}  # C5'-C4'-C3'-O3'
        self.epsilon: Dict[str, Optional[float]] = {}  # C4'-C3'-O3'-P(i+1)
        self.zeta: Dict[str, Optional[float]] = {}  # C3'-O3'-P(i+1)-O5'(i+1)
        self.chi: Dict[str, Optional[float]] = {}  # O4'-C1'-N1-C2 or O4'-C1'-N9-C4

        self.theta0: Dict[str, Optional[float]] = {}  # C4'-O4'-C1'-C2'
        self.theta1: Dict[str, Optional[float]] = {}  # O4'-C1'-C2'-C3'
        self.theta2: Dict[str, Optional[float]] = {}  # C1'-C2'-C3'-C4'
        self.theta3: Dict[str, Optional[float]] = {}  # C2'-C3'-C4'-O4'
        self.theta4: Dict[str, Optional[float]] = {}  # C3'-C4'-O4'-C1'

        self.alpha_conformation: Dict[str, Optional[str]] = {}  # ap/sc/-sc
        self.gamma_conformation: Dict[str, Optional[str]] = {}  # trans/gauche+/gauche-
        self.zeta_conformation: Dict[str, Optional[str]] = {}  # ap/sc/-sc
        self.chi_conformation: Dict[str, Optional[str]] = {}  # syn/anti

        self.tau_max: Dict[str, Optional[float]] = {}  # sugar pucker amplitude
        self.pseudorotation: Dict[str, Optional[float]] = {}  # the phase angle of pseudorotation
        self.sugar_conformation: Dict[str, Optional[str]] = {}  # sugar conformation (C2' endo or C3' endo)

    def pick_atoms(self, atom_name: str, relative_position: int):
        if relative_position == 0:
            relative_residue: Residue = self.residue_entry.residue
        elif relative_position == 1:
            relative_residue = self.residue_entry.get_next().residue
        elif relative_position == -1:
            relative_residue = self.residue_entry.get_prev().residue
        else:
            raise ValueError("Invalid relative_position")

        atom_group: Atom = relative_residue[atom_name]
        return atom_group.disordered_get_list() if atom_group.is_disordered() else [atom_group]

    @staticmethod
    def _round_torsion(atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom):
        torsion = calc_dihedral(atom1.get_vector(), atom2.get_vector(), atom3.get_vector(), atom4.get_vector())
        return round(np.rad2deg(torsion), 1)

    def _calculate_ordered_torsions(self, atom_names: List[str], atom_relative_positions: List[int]):
        try:
            atom1 = self.pick_atoms(atom_names[0], atom_relative_positions[0])[0]
            atom2 = self.pick_atoms(atom_names[1], atom_relative_positions[1])[0]
            atom3 = self.pick_atoms(atom_names[2], atom_relative_positions[2])[0]
            atom4 = self.pick_atoms(atom_names[3], atom_relative_positions[3])[0]
        except KeyError:
            return {"": None}
        return {"": self._round_torsion(atom1, atom2, atom3, atom4)}

    def _calculate_disordered_torsions(
        self, atom_names: List[str], atom_relative_positions: List[int]
    ) -> Dict[str, Optional[float]]:
        torsions = {}

        try:
            atoms1 = self.pick_atoms(atom_names[0], atom_relative_positions[0])
            atoms2 = self.pick_atoms(atom_names[1], atom_relative_positions[1])
            atoms3 = self.pick_atoms(atom_names[2], atom_relative_positions[2])
            atoms4 = self.pick_atoms(atom_names[3], atom_relative_positions[3])
        except KeyError:
            return {"": None}

        # pylint: disable=too-many-nested-blocks
        for atom1 in atoms1:
            for atom2 in atoms2:
                for atom3 in atoms3:
                    for atom4 in atoms4:
                        alt_locs = set(
                            [
                                atom1.get_altloc().strip(),
                                atom2.get_altloc().strip(),
                                atom3.get_altloc().strip(),
                                atom4.get_altloc().strip(),
                            ]
                        )
                        alt_locs.discard("")
                        # if not mixed (for example only "", or only one alternative fonformation "" and "A")
                        if len(alt_locs) <= 1:
                            torsion = self._round_torsion(atom1, atom2, atom3, atom4)
                            alt_loc = ""
                            if len(alt_locs) == 1:
                                alt_loc = alt_locs.pop()
                            torsions[alt_loc] = torsion
        return torsions

    def calculate_torsions(self, atom_names: List[str], atom_relative_positions: List[int]) -> Dict[str, Optional[float]]:
        if self.residue_entry.residue.is_disordered() == 0:
            return self._calculate_ordered_torsions(atom_names, atom_relative_positions)
        return self._calculate_disordered_torsions(atom_names, atom_relative_positions)

    @classmethod
    def _pseudorotation_with_sd(cls, theta0, theta1, theta2, theta3, theta4):
        """
        Calculate pseudorotation angle
        :param theta: list of theta torsion values, for example theta[0] = torsion(C4', O4' C1', C2')
        :return: P in deg, standard deviation of P, Tm, standard deviation of Tm
        """
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-locals

        # the initial definition  is Theta(1) = C1-C2-C3-C4, Theta(2) = C2-C3-C4-O4, etc.
        _theta = [theta2, theta3, theta4, theta0, theta1]

        sum_sin = 0.0
        sum_cos = 0.0

        for i_t, _t in enumerate(_theta):
            _x = 0.8 * math.pi * i_t
            sum_sin += _t * math.sin(_x)
            sum_cos += _t * math.cos(_x)

        pseudo_deg = math.degrees(math.atan2(-sum_sin, sum_cos))

        if pseudo_deg < 0.0:
            pseudo_deg += 360.0

        pseudo_rad = math.radians(pseudo_deg)
        _tm = 0.4 * (math.cos(pseudo_rad) * sum_cos - math.sin(pseudo_rad) * sum_sin)

        _st = 0.0
        thc = [0.0, 0.0, 0.0, 0.0, 0.0]

        for i_t, _t in enumerate(_theta):
            thc[i_t] = _tm * math.cos(pseudo_rad + (0.8 * math.pi * i_t))
            _d = _t - thc[i_t]
            _st += _d * _d

        sd_tm = math.sqrt(0.4 * _st / 3.0)
        sd_p = sd_tm / math.radians(_tm)
        return round(pseudo_deg, 1), sd_p, round(_tm, 1), sd_tm

    def calculate_alpha(self):
        self.alpha = self.calculate_torsions(("O3'", "P", "O5'", "C5'"), (-1, 0, 0, 0))

    def calculate_alpha_conformation(self):
        for alt_loc, angle in self.alpha.items():
            if angle:
                if 30 <= angle <= 110:
                    self.alpha_conformation[alt_loc] = "sc+"
                elif -110 <= angle <= -30:
                    self.alpha_conformation[alt_loc] = "sc-"
                elif angle <= -130 or angle > 110:
                    self.alpha_conformation[alt_loc] = "ap"
                else:
                    self.alpha_conformation[alt_loc] = "other"
            else:
                self.alpha_conformation[alt_loc] = "undefined"

    def calculate_beta(self):
        self.beta = self.calculate_torsions(("P", "O5'", "C5'", "C4'"), (0, 0, 0, 0))

    def calculate_gamma(self):
        self.gamma = self.calculate_torsions(("O5'", "C5'", "C4'", "C3'"), (0, 0, 0, 0))

    def calculate_gamma_conformation(self):
        for alt_loc, angle in self.gamma.items():
            if angle:
                if 30 <= angle <= 90:
                    self.gamma_conformation[alt_loc] = "gauche+"
                elif -90 <= angle <= -30:
                    self.gamma_conformation[alt_loc] = "gauche-"
                elif 150 <= angle or angle <= -150:
                    self.gamma_conformation[alt_loc] = "trans"
                else:
                    self.gamma_conformation[alt_loc] = "other"
            else:
                self.gamma_conformation[alt_loc] = "undefined"

    def calculate_delta(self):
        self.delta = self.calculate_torsions(("C5'", "C4'", "C3'", "O3'"), (0, 0, 0, 0))

    def calculate_epsilon(self):
        self.epsilon = self.calculate_torsions(("C4'", "C3'", "O3'", "P"), (0, 0, 0, 1))

    def calculate_zeta(self):
        self.zeta = self.calculate_torsions(("C3'", "O3'", "P", "O5'"), (0, 0, 1, 1))

    def calculate_zeta_conformation(self):
        for alt_loc, angle in self.zeta.items():
            if angle:
                if 30 <= angle <= 110:
                    self.zeta_conformation[alt_loc] = "sc+"
                elif -110 <= angle <= -30:
                    self.zeta_conformation[alt_loc] = "sc-"
                elif angle <= -130 or angle > 110:
                    self.zeta_conformation[alt_loc] = "ap"
                else:
                    self.zeta_conformation[alt_loc] = "other"
            else:
                self.zeta_conformation[alt_loc] = "undefined"

    def calculate_theta_and_pseudorotation(self):
        self.theta0 = self.calculate_torsions(("C4'", "O4'", "C1'", "C2'"), (0, 0, 0, 0))
        self.theta1 = self.calculate_torsions(("O4'", "C1'", "C2'", "C3'"), (0, 0, 0, 0))
        self.theta2 = self.calculate_torsions(("C1'", "C2'", "C3'", "C4'"), (0, 0, 0, 0))
        self.theta3 = self.calculate_torsions(("C2'", "C3'", "C4'", "O4'"), (0, 0, 0, 0))
        self.theta4 = self.calculate_torsions(("C3'", "C4'", "O4'", "C1'"), (0, 0, 0, 0))
        self.calculate_pseudorotation()

    def calculate_pseudorotation(self):
        alt_locs = set(self.theta0.keys())
        alt_locs.update(self.theta1.keys())
        alt_locs.update(self.theta2.keys())
        alt_locs.update(self.theta3.keys())
        alt_locs.update(self.theta4.keys())

        if "" in alt_locs and len(alt_locs) == 1:
            if (
                self.theta0[""] is not None
                and self.theta1[""] is not None
                and self.theta2[""] is not None
                and self.theta3[""] is not None
                and self.theta4[""] is not None
            ):
                pseudorotation, _, tau_max, _ = self._pseudorotation_with_sd(
                    self.theta0[""], self.theta1[""], self.theta2[""], self.theta3[""], self.theta4[""]
                )
                self.pseudorotation[""] = pseudorotation
                self.tau_max[""] = tau_max
            else:
                self.pseudorotation[""] = None
                self.tau_max[""] = None

        alt_locs.discard("")
        for alt_loc in alt_locs:
            if (
                self.theta0.get(alt_loc, self.theta0.get("", None)) is not None
                and self.theta1.get(alt_loc, self.theta1.get("", None)) is not None
                and self.theta2.get(alt_loc, self.theta2.get("", None)) is not None
                and self.theta3.get(alt_loc, self.theta3.get("", None)) is not None
                and self.theta4.get(alt_loc, self.theta4.get("", None)) is not None
            ):
                pseudorotation, _, tau_max, _ = self._pseudorotation_with_sd(
                    self.theta0.get(alt_loc, self.theta0.get("", None)),
                    self.theta1.get(alt_loc, self.theta1.get("", None)),
                    self.theta2.get(alt_loc, self.theta2.get("", None)),
                    self.theta3.get(alt_loc, self.theta3.get("", None)),
                    self.theta4.get(alt_loc, self.theta4.get("", None)),
                )
                self.pseudorotation[alt_loc] = pseudorotation
                self.tau_max[alt_loc] = tau_max
            else:
                self.pseudorotation[alt_loc] = None
                self.tau_max[alt_loc] = None

    def calulate_sugar_conformation(self):
        for alt_loc, angle in self.pseudorotation.items():
            if angle:
                if 140 <= angle <= 190:
                    self.sugar_conformation[alt_loc] = "C2'-endo"
                elif 0 <= angle <= 36:
                    self.sugar_conformation[alt_loc] = "C3'-endo"
                else:
                    self.sugar_conformation[alt_loc] = "other"
            else:
                self.sugar_conformation[alt_loc] = "undefined"

    def calculate_chi(self):
        atoms_names = ["O4'", "C1'", "N1", "C2"]
        if self.residue_entry.res_name in PURINES_RES_NAMES:
            atoms_names = ["O4'", "C1'", "N9", "C4"]
        self.chi = self.calculate_torsions(atoms_names, (0, 0, 0, 0))

    def calculate_chi_conformation(self):
        for alt_loc, angle in self.chi.items():
            if angle:
                if -90 <= angle <= 90:
                    self.chi_conformation[alt_loc] = "syn"
                else:
                    self.chi_conformation[alt_loc] = "anti"
            else:
                self.chi_conformation[alt_loc] = "undefined"

    def calculate_conformation(self):
        self.calculate_alpha()
        self.calculate_alpha_conformation()
        self.calculate_beta()
        self.calculate_gamma()
        self.calculate_gamma_conformation()
        self.calculate_delta()
        self.calculate_epsilon()
        self.calculate_zeta()
        self.calculate_zeta_conformation()
        self.calculate_theta_and_pseudorotation()
        self.calulate_sugar_conformation()
        self.calculate_chi()
        self.calculate_chi_conformation()

    @staticmethod
    def _print_torsion(name, torsion, conformation=None):
        for alt_loc in sorted(torsion.keys()):
            _alt_loc = alt_loc if alt_loc else " "
            _conformation = conformation[alt_loc] if conformation else ""
            print(f"#\t {name} {_alt_loc} {torsion[alt_loc]} {_conformation}")

    def prepare_report_torsion(self):
        print(
            f"# Chain: {self.residue_entry.chain.get_id()}, "
            + f"Residue id: {self.residue_entry.resseq}, "
            + f"Residue name: {self.residue_entry.res_name}"
        )
        self._print_torsion("Alpha   (O3'(i-1)-P-O5'-C5'):         ", self.alpha, self.alpha_conformation)
        self._print_torsion("Beta    (P-O5'-C5'-C4'):              ", self.beta)
        self._print_torsion("Gamma   (O5'-C5'-C4'-C3'):            ", self.gamma, self.gamma_conformation)
        self._print_torsion("Delta   (C5'-C4'-C3'-O3'):            ", self.delta)
        self._print_torsion("Epsilon (C4'-C3'-O3'-P(i+1)):         ", self.epsilon)
        self._print_torsion("Zeta    (C3'-O3'-P(i+1)-O5'(i+1)):    ", self.zeta, self.zeta_conformation)
        self._print_torsion("Chi     (O4'-C1'-N1-C2/O4'-C1'-N9-C4):", self.chi, self.chi_conformation)
        self._print_torsion("Theta0  (C4'-O4'-C1'-C2'):            ", self.theta0)
        self._print_torsion("Theta1  (O4'-C1'-C2'-C3'):            ", self.theta1)
        self._print_torsion("Theta2  (C1'-C2'-C3'-C4'):            ", self.theta2)
        self._print_torsion("Theta3  (C2'-C3'-C4'-O4'):            ", self.theta3)
        self._print_torsion("Theta4  (C3'-C4'-O4'-C1'):            ", self.theta4)
        self._print_torsion("Tau_max:                              ", self.tau_max)
        self._print_torsion("Pseudorotation:                       ", self.pseudorotation, self.sugar_conformation)
