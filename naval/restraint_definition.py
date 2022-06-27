from typing import Optional


class BondDefinition:
    """
    Simple container class for bond definitions
    """

    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-instance-attributes

    __slots__ = (
        "name",
        "atom1",
        "atom2",
        "csd_target",
        "csd_std",
        "pdb_count",
        "pdb_mean",
        "pdb_std",
        "pdb_3low",
        "pdb_3high",
        "pdb_4low",
        "pdb_4high",
    )

    def __init__(
        self,
        name: str,
        atom1: str,
        atom2: str,
        csd_target: float,
        csd_std: float,
        pdb_count: Optional[int],
        pdb_mean: Optional[float],
        pdb_std: Optional[float],
        pdb_3low: Optional[float],
        pdb_3high: Optional[float],
        pdb_4low: Optional[float],
        pdb_4high: Optional[float],
    ):
        # pylint: disable=too-many-arguments
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2
        self.csd_target = csd_target
        self.csd_std = csd_std
        self.pdb_count = pdb_count
        self.pdb_mean = pdb_mean
        self.pdb_std = pdb_std
        self.pdb_3low = pdb_3low
        self.pdb_3high = pdb_3high
        self.pdb_4low = pdb_4low
        self.pdb_4high = pdb_4high


class AngleDefinition:
    """
    Simple container class for angle definitions
    """

    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-instance-attributes

    __slots__ = (
        "name",
        "atom1",
        "atom2",
        "atom3",
        "csd_target",
        "csd_std",
        "pdb_count",
        "pdb_mean",
        "pdb_std",
        "pdb_3low",
        "pdb_3high",
        "pdb_4low",
        "pdb_4high",
    )

    def __init__(
        self,
        name: str,
        atom1: str,
        atom2: str,
        atom3: str,
        csd_target: float,
        csd_std: float,
        pdb_count: Optional[int],
        pdb_mean: Optional[float],
        pdb_std: Optional[float],
        pdb_3low: Optional[float],
        pdb_3high: Optional[float],
        pdb_4low: Optional[float],
        pdb_4high: Optional[float],
    ):
        # pylint: disable=too-many-arguments
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.csd_target = csd_target
        self.csd_std = csd_std
        self.pdb_count = pdb_count
        self.pdb_mean = pdb_mean
        self.pdb_std = pdb_std
        self.pdb_3low = pdb_3low
        self.pdb_3high = pdb_3high
        self.pdb_4low = pdb_4low
        self.pdb_4high = pdb_4high
