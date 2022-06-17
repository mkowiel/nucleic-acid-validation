from unittest.mock import Mock

from naval.validation_record import ValidationRecord


def test_is_csd():
    record = ValidationRecord("bond", Mock(), Mock(), Mock(), Mock(), 1.43, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.in_csd is True
    assert record.in_pdb_lv3 is False
    assert record.in_pdb_lv4 is False
    assert record.outlier is False


def test_is_in_pdb_lv3():
    record = ValidationRecord("bond", Mock(), Mock(), Mock(), Mock(), 1.49, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.in_csd is False
    assert record.in_pdb_lv3 is True
    assert record.in_pdb_lv4 is False
    assert record.outlier is False


def test_is_in_pdb_lv4():
    record = ValidationRecord("bond", Mock(), Mock(), Mock(), Mock(), 1.54, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.in_csd is False
    assert record.in_pdb_lv3 is False
    assert record.in_pdb_lv4 is True
    assert record.outlier is False


def test_is_outlier():
    record = ValidationRecord("bond", Mock(), Mock(), Mock(), Mock(), 1.60, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.in_csd is False
    assert record.in_pdb_lv3 is False
    assert record.in_pdb_lv4 is False
    assert record.outlier is True
