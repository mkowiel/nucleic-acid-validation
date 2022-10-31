from unittest.mock import Mock

from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

from naval.nucleotide_geometry import NucleotideGeometry
from naval.printer import AnglesCsvPrinter, BondsCsvPrinter, GeometryCsvPrinter
from naval.residue_cache_entry import ResidueCacheEntry
from naval.validation_record import TorsionRecord, ValidationRecord


def prepare_validation_records():
    res = Residue((" ", 1, " "), "A", 1)
    atom1 = Atom("OP1", (1.6, 0, 0), 10, 1.0, " ", "OP1", 1, "O")
    atom2 = Atom("P", (0, 0, 0), 10, 1.0, " ", "P", 2, "P")
    atom3 = Atom("OP2", (0, 0, 1.6), 10, 1.0, " ", "OP2", 3, "O")

    res.add(atom1)
    res.add(atom2)
    res.add(atom3)

    records = [
        ValidationRecord("bond", "test", Mock(), atom1, atom2, None, 1.43, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55),
        ValidationRecord("angle", "test", Mock(), atom1, atom2, atom3, 1.55, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55),
    ]
    return records


def test_bonds_csv_printer():
    records = prepare_validation_records()

    printer = BondsCsvPrinter()
    lines = printer.print(records)
    assert len(lines) == 1 + 1
    assert "bond" in lines[1]


def test_angle_csv_printer():
    records = prepare_validation_records()

    printer = AnglesCsvPrinter()
    lines = printer.print(records)
    assert len(lines) == 1 + 1
    assert "angle" in lines[1]


def prepare_torsion_records():
    res = Residue((" ", 1, " "), "A", 1)
    residue_cache_entry = ResidueCacheEntry("1aa1", Mock(), Mock(), res)
    geometry = NucleotideGeometry(residue_cache_entry)

    records = [
        TorsionRecord("torsion", "test", geometry, "A", 180, "OK"),
    ]
    return records


def test_geometry_csv_printer():
    records = prepare_torsion_records()

    printer = GeometryCsvPrinter()
    lines = printer.print(records)
    assert len(lines) == 1 + 1
    assert "torsion" in lines[1]
