from naval.nucleotide_geometry import NucleotideGeometry


class Validator:
    """
    Base validator class
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, geometry: NucleotideGeometry, csd_sig: float = 3) -> None:
        self.geometry = geometry
        self.csd_sig = csd_sig

    def validate(self):
        pass
