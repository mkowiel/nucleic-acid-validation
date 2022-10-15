from typing import List

from naval.nucleotide_geometry import NucleotideGeometry
from naval.validation_record import TorsionRecord


class GeometryValidator:
    """
    Base validator class
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, geometry: NucleotideGeometry) -> None:
        self.geometry = geometry

    def _add_torsions(self, torsion_type, name, torsion, conformation):
        for alt_loc in sorted(torsion.keys()):
            _alt_loc = alt_loc if alt_loc else ""
            _conformation = conformation[alt_loc] if conformation else ""
            yield TorsionRecord(torsion_type, name, self.geometry, _alt_loc, torsion[alt_loc], _conformation)

    def _validate_torsion(self) -> List[TorsionRecord]:
        records = []

        records.extend(self._add_torsions("torsion", "alpha", self.geometry.alpha, self.geometry.alpha_conformation))
        records.extend(self._add_torsions("torsion", "beta", self.geometry.beta, None))
        records.extend(self._add_torsions("torsion", "gamma", self.geometry.gamma, self.geometry.gamma_conformation))
        records.extend(self._add_torsions("torsion", "delta", self.geometry.delta, None))
        records.extend(self._add_torsions("torsion", "epsilon", self.geometry.epsilon, None))
        records.extend(self._add_torsions("torsion", "zeta", self.geometry.zeta, self.geometry.zeta_conformation))
        records.extend(self._add_torsions("torsion", "chi", self.geometry.chi, self.geometry.chi_conformation))
        records.extend(self._add_torsions("torsion", "theta0", self.geometry.theta0, None))
        records.extend(self._add_torsions("torsion", "theta1", self.geometry.theta1, None))
        records.extend(self._add_torsions("torsion", "theta2", self.geometry.theta2, None))
        records.extend(self._add_torsions("torsion", "theta3", self.geometry.theta3, None))
        records.extend(self._add_torsions("torsion", "theta4", self.geometry.theta4, None))
        records.extend(self._add_torsions("pseudorotation", "tau_max", self.geometry.tau_max, None))
        records.extend(
            self._add_torsions(
                "pseudorotation", "pseudorotation", self.geometry.pseudorotation, self.geometry.sugar_conformation
            )
        )
        return records

    def validate(self) -> List[TorsionRecord]:
        records = self._validate_torsion()
        return records
