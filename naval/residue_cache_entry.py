from __future__ import annotations

from typing import Optional

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from naval.nucleotide_definitions import NUCLEOTIDE_RES_NAMES

if False:  # pylint: disable=using-constant-test
    # trick for mypy to avoid cyclic imports
    # pylint: disable=cyclic-import
    from naval.nucleotide_geometry import NucleotideGeometry


class ResidueCacheEntry:
    """
    Class to keep cache torsion angles for given residue
    """

    # pylint: disable=too-many-public-methods
    # pylint: disable=too-many-instance-attributes
    def __init__(self, pdbcode: str, model, chain: Chain, residue: Residue) -> None:
        """Simple container class to keep torsion angles of the residue.
        Calculates torison anles for all alternative conformations.
        Maps torsion angles to string classes.
        """
        self.pdbcode = pdbcode
        self.model = model
        self.chain = chain
        self.residue = residue
        self.res_name = self.residue.get_resname()
        self.res_full_id = self.residue.get_id()
        self.resseq = self.res_full_id[1]
        self.inscode = self.res_full_id[2]

        # TODO, move to Geometry/Residue base class
        self.next_res: Optional[ResidueCacheEntry] = None
        self.prev_res: Optional[ResidueCacheEntry] = None
        self.geometry: "Optional[NucleotideGeometry]" = None

    def is_terminal(self):
        return not self.has_next() or not self.has_prev()

    def has_next(self):
        return self.next_res is not None

    def has_prev(self):
        return self.prev_res is not None

    def get_next(self) -> ResidueCacheEntry:
        if self.next_res is not None:
            return self.next_res
        raise KeyError("Does not have next residue cache entry")

    def get_prev(self) -> ResidueCacheEntry:
        if self.prev_res is not None:
            return self.prev_res
        raise KeyError("Does not have prev residue cache entry")

    def is_nucleotide(self):
        return self.res_name in NUCLEOTIDE_RES_NAMES
