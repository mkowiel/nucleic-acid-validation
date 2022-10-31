from unittest.mock import Mock

from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import Vector

from naval.nucleotide_geometry import NucleotideGeometry
from naval.residue_cache_entry import ResidueCacheEntry
from naval.validators.bases_validator import BasesValidator
from naval.validators.po4_validator import Po4Validator
from naval.validators.sugar_basic_validator import BasicSugarValidator
from naval.validators.sugar_pucker_validator import SugarPuckerBasedSugarValidator


def prepare_geometry(resname, atom1, atom2, atom3, atom1_element, atom2_element, atom3_element):
    # pylint: disable=too-many-arguments
    res = Residue((" ", 1, " "), resname, 1)
    atom1 = Atom(atom1, Vector(1.45, 0, 0.00), 10, 1.0, " ", atom1, 1, atom1_element)
    atom2 = Atom(atom2, Vector(0.00, 0, 0.00), 10, 1.0, " ", atom2, 2, atom2_element)
    atom3 = Atom(atom3, Vector(0.00, 0, 1.45), 10, 1.0, " ", atom3, 3, atom3_element)

    res.add(atom1)
    res.add(atom2)
    res.add(atom3)

    residue_cache_entry = ResidueCacheEntry("1aa1", Mock(), Mock(), res)
    return NucleotideGeometry(residue_cache_entry)


def test_bases_validator():
    geometry = prepare_geometry("A", "C4", "C5", "C6", "C", "C", "C")

    validator = BasesValidator(geometry)
    validation_records = validator.validate()
    assert len(validation_records) == 3


def test_po4_validator():
    geometry = prepare_geometry("A", "OP1", "P", "OP2", "O", "P", "O")

    validator = Po4Validator(geometry)
    validation_records = validator.validate()
    assert len(validation_records) == 3


def test_basic_sugar_validator():
    for resname in ["A", "DA", "U", "DU"]:
        geometry = prepare_geometry(resname, "C1'", "C2'", "C3'", "C", "C", "C")

        validator = BasicSugarValidator(geometry)
        validation_records = validator.validate()
        assert len(validation_records) == 3


def test_sugar_pucker_based_sugar_validator():
    for resname in ["A", "DA", "U", "DU"]:
        for conforamtion in ["C2'-endo", "C3'-endo", "other", "undefined", None]:
            geometry = prepare_geometry(resname, "C1'", "C2'", "C3'", "C", "C", "C")
            geometry.sugar_conformation = {"": conforamtion}

            validator = SugarPuckerBasedSugarValidator(geometry)
            validation_records = validator.validate()
            assert len(validation_records) == 3
