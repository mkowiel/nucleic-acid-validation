import os

from naval.printer import CsvPrinter
from naval.validate import read_structure, validate_structure


def test_read_structure():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    assert struct.id == "1d8g"


def filter_records_resseq(records, record_type, chain, resseq):
    return [
        record
        for record in records
        if (
            record_type == record.validation_type
            and chain == record.geometry.residue_entry.chain.get_id()
            and (resseq == record.geometry.residue_entry.resseq)
        )
    ]


def filter_records_atoms(records, record_type, chain, atom1, atom2, atom3):
    # pylint: disable=too-many-arguments
    return [
        record
        for record in records
        if (
            record_type == record.validation_type
            and chain == record.geometry.residue_entry.chain.get_id()
            and (atom1 == record.atom1.name)
            and (atom2 == record.atom2.name)
            and (atom3 == record.atom3.name)
        )
    ]


def filter_records_bond(records, chain, resseq1, atom1, resseq2, atom2):
    # pylint: disable=too-many-arguments
    return [
        record
        for record in records
        if (
            "bond" == record.validation_type
            and chain == record.geometry.residue_entry.chain.get_id()
            and (atom1 == record.atom1.name)
            and (atom2 == record.atom2.name)
            and (resseq1 == record.atom1.get_parent().get_id())
            and (resseq2 == record.atom2.get_parent().get_id())
        )
    ]


def test_iterate_struct_1d8g():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    records = validate_structure(struct)
    assert len(records) > 0

    print("\n".join(CsvPrinter().print(records)[0:50]))

    # 8 (base) + (C3'-O3')*2(disorder) + (C5'-O5')*2(disorder) + 7 (sugar) * 2(disorder)
    assert len(filter_records_resseq(records, "bond", "A", 1)) == (8 + 2 * 2 + (7 * 2))
    # 10 (base) + 0*2 (PO4) + 14 (sugar) * 2(disorder)
    assert len(filter_records_resseq(records, "angle", "A", 1)) == 10 + 0 + 14 * 2

    assert len(filter_records_atoms(records, "angle", "A", "C6", "N1", "C2")) == 11
    assert len(filter_records_atoms(records, "angle", "A", "OP1", "P", "OP2")) == 15


def test_iterate_struct_5hr7():
    struct = read_structure(os.path.dirname(__file__) + "/examples/5hr7.pdb")
    records = validate_structure(struct)
    assert len(records) > 0

    assert len(filter_records_atoms(records, "angle", "D", "OP1", "P", "OP2")) == 71 - 2 + 1
    assert len(filter_records_atoms(records, "angle", "D", "OP1", "P", "O3'")) == 71 - 2 + 1 - 1

    assert len(filter_records_bond(records, "D", (" ", 19, " "), "OP1", (" ", 19, " "), "P")) == 1

    # if atoms are taken from the same res then it will be greater than 0
    assert len(filter_records_bond(records, "D", (" ", 19, " "), "O3'", (" ", 19, " "), "P")) == 0
    assert len(filter_records_bond(records, "D", (" ", 19, " "), "O3'", (" ", 20, " "), "P")) == 1

    # we between residue 20 and 21 there is 20A inserted, we should have one bond between 20 and 20A and 20A and 21
    assert len(filter_records_bond(records, "D", (" ", 20, " "), "O3'", (" ", 21, " "), "P")) == 0
    assert len(filter_records_bond(records, "D", (" ", 20, " "), "O3'", (" ", 20, "A"), "P")) == 1
    assert len(filter_records_bond(records, "D", (" ", 20, "A"), "O3'", (" ", 21, " "), "P")) == 1
