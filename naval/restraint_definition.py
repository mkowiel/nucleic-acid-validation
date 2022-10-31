class BondDefinition:
    """
    Simple container class for bond definitions
    """

    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-instance-attributes

    __slots__ = (
        "name",
        "atom1_name",
        "atom2_name",
        "atom1_relative_res_position",
        "atom2_relative_res_position",
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
        atom1_name: str,
        atom2_name: str,
        atom1_relative_res_position: int,
        atom2_relative_res_position: int,
        csd_target: float,
        csd_std: float,
        pdb_count: int,
        pdb_mean: float,
        pdb_std: float,
        # TODO: rename
        pdb_3low: float,
        pdb_3high: float,
        pdb_4low: float,
        pdb_4high: float,
    ):
        # pylint: disable=too-many-arguments
        self.name = name
        self.atom1_name = atom1_name
        self.atom2_name = atom2_name
        self.atom1_relative_res_position = atom1_relative_res_position
        self.atom2_relative_res_position = atom2_relative_res_position
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
        "atom1_name",
        "atom2_name",
        "atom3_name",
        "atom1_relative_res_position",
        "atom2_relative_res_position",
        "atom3_relative_res_position",
        "csd_target",
        "csd_std",
        "pdb_count",
        "pdb_mean",
        "pdb_std",
        # TODO rename
        "pdb_3low",
        "pdb_3high",
        "pdb_4low",
        "pdb_4high",
    )

    def __init__(
        self,
        name: str,
        atom1_name: str,
        atom2_name: str,
        atom3_name: str,
        atom1_relative_res_position: int,
        atom2_relative_res_position: int,
        atom3_relative_res_position: int,
        csd_target: float,
        csd_std: float,
        pdb_count: int,
        pdb_mean: float,
        pdb_std: float,
        pdb_3low: float,
        pdb_3high: float,
        pdb_4low: float,
        pdb_4high: float,
    ):
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-locals
        self.name = name
        self.atom1_name = atom1_name
        self.atom2_name = atom2_name
        self.atom3_name = atom3_name
        self.atom1_relative_res_position = atom1_relative_res_position
        self.atom2_relative_res_position = atom2_relative_res_position
        self.atom3_relative_res_position = atom3_relative_res_position
        self.csd_target = csd_target
        self.csd_std = csd_std
        self.pdb_count = pdb_count
        self.pdb_mean = pdb_mean
        self.pdb_std = pdb_std
        self.pdb_3low = pdb_3low
        self.pdb_3high = pdb_3high
        self.pdb_4low = pdb_4low
        self.pdb_4high = pdb_4high
