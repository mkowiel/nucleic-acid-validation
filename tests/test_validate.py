import os
from naval.validate import validate_structure
from naval.validate import read_structure
from naval.printer import CsvPrinter


def test_read_structure():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    assert struct.id == "1d8g"


def filter_records(records, record_type, chain, resseq):
    return [
        record
        for record in records
        if (
            record_type == record.validation_type
            and chain == record.geometry.chain.get_id()
            and (resseq == record.geometry.resseq)
        )
    ]


def test_iterate_struct():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    records = validate_structure(struct)
    assert len(records) > 0

    assert len(filter_records(records, "bond", "A", 1)) == 8
    assert len(filter_records(records, "angle", "A", 1)) == 10

    print("\n".join(CsvPrinter().print(records)))
