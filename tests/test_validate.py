import os
from naval.validate import iterate_struct
from naval.validate import read_structure
from naval.printer import CsvPrinter


def test_read_structure():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    assert struct.id == "1d8g"


def test_iterate_struct():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    records = iterate_struct(struct)
    assert len(records) > 0
    print("\n".join(CsvPrinter().print(records)))
