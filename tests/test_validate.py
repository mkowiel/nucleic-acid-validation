import os
from naval.validate import iterate_struct
from naval.validate import read_structure


def test_read_structure():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    assert struct.id == "1d8g"


def test_iterate_struct():
    struct = read_structure(os.path.dirname(__file__) + "/examples/1d8g.pdb")
    assert len(iterate_struct(struct)) > 0
