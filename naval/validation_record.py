from typing import Optional

from Bio.PDB.Atom import Atom

from naval.nucleotide_geometry import NucleotideGeometry


class ValidationRecord:
    """
    Container class to keep the results of the validation
    """

    # pylint: disable=too-many-instance-attributes

    __slots__ = (
        "validation_type",
        "name",
        "geometry",
        "atom1",
        "atom2",
        "atom3",
        "calculated_value",
        "target_value",
        "target_sigma",
        "csd_left",
        "csd_right",
        "pdb_lv3_left",
        "pdb_lv3_right",
        "pdb_lv4_left",
        "pdb_lv4_right",
        "in_csd",
        "in_pdb_lv3",
        "in_pdb_lv4",
        "outlier",
    )

    # pylint: disable=too-many-arguments
    def __init__(
        self,
        validation_type: str,
        name: str,
        geometry: NucleotideGeometry,
        atom1: Atom,
        atom2: Atom,
        atom3: Atom,
        calculated_value: float,
        target_value: float,
        target_sigma: float,
        # percentiles equvalent of 3 sigma 0.9973% of population (1 per 370)
        pdb_lv3_left: Optional[float],
        pdb_lv3_right: Optional[float],
        # percentiles equvalent of 4 sigma 0.9999% of population (1 per 15787)
        pdb_lv4_left: Optional[float],
        pdb_lv4_right: Optional[float],
    ) -> None:
        if validation_type not in ("angle", "bond"):
            raise ValueError("Validation type nees to one of ['angle', 'bond']")
        self.validation_type = validation_type
        self.name = name
        self.geometry = geometry
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.calculated_value = calculated_value
        self.target_value = target_value
        self.target_sigma = target_sigma

        self.csd_left = self.target_value - 3 * self.target_sigma
        self.csd_right = self.target_value + 3 * self.target_sigma

        self.pdb_lv3_left = pdb_lv3_left
        self.pdb_lv3_right = pdb_lv3_right

        self.pdb_lv4_left = pdb_lv4_left
        self.pdb_lv4_right = pdb_lv4_right

        self.in_csd = self.is_in_csd_range()
        self.in_pdb_lv3 = self.is_in_pdb_lv3()
        self.in_pdb_lv4 = self.is_in_pdb_lv4()
        self.outlier = self.is_outlier()

    def is_in_csd_range(self):
        return self.csd_left <= self.calculated_value <= self.csd_right

    def is_in_pdb_lv3(self):
        return False if self.in_csd else (self.pdb_lv3_left <= self.calculated_value <= self.pdb_lv3_right)

    def is_in_pdb_lv4(self):
        return (
            False
            if (self.in_csd or self.in_pdb_lv3)
            else (self.pdb_lv4_left <= self.calculated_value <= self.pdb_lv4_right)
        )

    def is_outlier(self):
        return not (self.in_pdb_lv4 or self.in_pdb_lv3 or self.in_csd)
