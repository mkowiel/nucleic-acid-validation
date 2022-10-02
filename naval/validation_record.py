from typing import Optional

from Bio.PDB.Atom import Atom

from naval.nucleotide_geometry import NucleotideGeometry


class ValidationRecord:
    """
    Container class to keep the results of the validation
    """

    # TODO: we could link to restraint definition (BondDefinition and AngleDefinition)
    #       instead of coping the raw data

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
        "csd_preferred_left",
        "csd_preferred_right",
        "pdb_allowed_left",
        "pdb_allowed_right",
        "pdb_suspicious_left",
        "pdb_suspicious_right",
        "preferred",
        "allowed",
        "suspicious",
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
        pdb_allowed_left: Optional[float],
        pdb_allowed_right: Optional[float],
        # percentiles equvalent of 4 sigma 0.9999% of population (1 per 15787)
        pdb_suspicious_left: Optional[float],
        pdb_suspicious_right: Optional[float],
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

        self.csd_preferred_left = self.target_value - 3 * self.target_sigma
        self.csd_preferred_right = self.target_value + 3 * self.target_sigma

        self.pdb_allowed_left = pdb_allowed_left
        self.pdb_allowed_right = pdb_allowed_right

        self.pdb_suspicious_left = pdb_suspicious_left
        self.pdb_suspicious_right = pdb_suspicious_right

        self.preferred = self.is_preferred()
        self.allowed = self.is_allowed()
        self.suspicious = self.is_suspicious()
        self.outlier = self.is_outlier()

    def __str__(self):
        return (
            f"{self.validation_type} {self.name} {self.atom1} {self.atom2}"
            f" {self.atom3} {self.calculated_value:.3f} {self.target_value}"
        )

    def is_preferred(self):
        return self.csd_preferred_left <= self.calculated_value <= self.csd_preferred_right

    def is_allowed(self):
        return False if self.preferred else (self.pdb_allowed_left <= self.calculated_value <= self.pdb_allowed_right)

    def is_suspicious(self):
        return (
            False
            if (self.preferred or self.allowed)
            else (self.pdb_suspicious_left <= self.calculated_value <= self.pdb_suspicious_right)
        )

    def is_outlier(self):
        return not (self.suspicious or self.allowed or self.preferred)
