import os
from naval.validate import validate_structure
from naval.validate import read_structure
from naval.printer import CsvPrinter


def test_read_structure():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    assert struct.id == "1d8g"


def filter_records_resseq(records, record_type, chain, resseq):
    return [
        record
        for record in records
        if (
            record_type == record.validation_type
            and chain == record.geometry.chain.get_id()
            and (resseq == record.geometry.resseq)
        )
    ]


def filter_records_atoms(records, record_type, chain, atom1, atom2, atom3):
    # pylint: disable=too-many-arguments
    return [
        record
        for record in records
        if (
            record_type == record.validation_type
            and chain == record.geometry.chain.get_id()
            and (atom1 == record.atom1.name)
            and (atom2 == record.atom2.name)
            and (atom3 == record.atom3.name)
        )
    ]


def test_iterate_struct():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    records = validate_structure(struct)
    assert len(records) > 0

    print("\n".join(CsvPrinter().print(records)))

    # 8 (base) + (C3'-O3')*2(disorder) + (C5'-O5')*2(disorder) + sugar * 2(disorder)
    assert len(filter_records_resseq(records, "bond", "A", 1)) == (8 + 2 * 2)  # sugar * 2(disorder)
    # 10 (base) + 0 (PO4) + sugar * 2(disorder)
    assert len(filter_records_resseq(records, "angle", "A", 1)) == 10  # (base == 10)

    assert len(filter_records_atoms(records, "angle", "A", "C6", "N1", "C2")) == 11
    assert len(filter_records_atoms(records, "angle", "A", "OP1", "P", "OP2")) == 15
