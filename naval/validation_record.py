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
        pdb_allowed_left: float,
        pdb_allowed_right: float,
        # percentiles equvalent of 4 sigma 0.9999% of population (1 per 15787)
        pdb_suspicious_left: float,
        pdb_suspicious_right: float,
    ) -> None:
        if validation_type not in ("angle", "bond"):
            raise ValueError("Validation type nees to one of ['angle', 'bond']")
        self.validation_type: str = validation_type
        self.name: str = name
        self.geometry: NucleotideGeometry = geometry
        self.atom1: Atom = atom1
        self.atom2: Atom = atom2
        self.atom3: Optional[Atom] = atom3
        self.calculated_value: float = calculated_value
        self.target_value: float = target_value
        self.target_sigma: float = target_sigma

        self.csd_preferred_left: float = self.target_value - 3 * self.target_sigma
        self.csd_preferred_right: float = self.target_value + 3 * self.target_sigma

        self.pdb_allowed_left: float = pdb_allowed_left
        self.pdb_allowed_right: float = pdb_allowed_right

        self.pdb_suspicious_left: float = pdb_suspicious_left
        self.pdb_suspicious_right: float = pdb_suspicious_right

    def __str__(self) -> str:
        return (
            f"{self.validation_type} {self.name} {self.atom1} {self.atom2}"
            f" {self.atom3} {self.calculated_value:.3f} {self.target_value}"
        )

    def is_preferred(self) -> bool:
        return self.csd_preferred_left <= self.calculated_value <= self.csd_preferred_right

    def is_allowed(self) -> bool:
        return False if self.is_preferred() else (self.pdb_allowed_left <= self.calculated_value <= self.pdb_allowed_right)

    def is_suspicious(self) -> bool:
        return (
            False
            if (self.is_preferred() or self.is_allowed())
            else (self.pdb_suspicious_left <= self.calculated_value <= self.pdb_suspicious_right)
        )

    def is_outlier(self) -> bool:
        return not (self.is_suspicious() or self.is_allowed() or self.is_preferred())

    @property
    def label(self) -> str:
        if self.is_preferred():
            return "CSD-preferred"  # CSD-acceptable ?
        if self.is_allowed():
            return "PDB-acceptable"
        if self.is_suspicious():
            return "PDB-suspicious"
        return "PDB-outlier"


class TorsionRecord:
    """
    Container class to keep the results of the torsion angles
    """

    # pylint: disable=too-few-public-methods

    __slots__ = (
        "validation_type",
        "name",
        "geometry",
        "alt_loc",
        "calculated_value",
        "calculated_value_label",
    )

    def __init__(
        self,
        validation_type: str,
        name: str,
        geometry: NucleotideGeometry,
        alt_loc: str,
        calculated_value: float,
        calculated_value_label: str,
    ) -> None:
        # pylint: disable=too-many-arguments
        if validation_type not in ("torsion", "pseudorotation"):
            raise ValueError("Validation type nees to one of ['torsion', 'pseudorotation']")
        self.validation_type = validation_type
        self.name = name
        self.geometry = geometry
        self.alt_loc = alt_loc
        self.calculated_value = calculated_value
        self.calculated_value_label = calculated_value_label
