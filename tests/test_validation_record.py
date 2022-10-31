from unittest.mock import Mock

from naval.validation_record import ValidationRecord


def test_is_preferred():
    record = ValidationRecord("bond", "test", Mock(), Mock(), Mock(), Mock(), 1.43, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.is_preferred() is True
    assert record.is_allowed() is False
    assert record.is_suspicious() is False
    assert record.is_outlier() is False


def test_is_allowed():
    record = ValidationRecord("bond", "test", Mock(), Mock(), Mock(), Mock(), 1.49, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.is_preferred() is False
    assert record.is_allowed() is True
    assert record.is_suspicious() is False
    assert record.is_outlier() is False


def test_is_suspicious():
    record = ValidationRecord("bond", "test", Mock(), Mock(), Mock(), Mock(), 1.54, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.is_preferred() is False
    assert record.is_allowed() is False
    assert record.is_suspicious() is True
    assert record.is_outlier() is False


def test_is_outlier():
    record = ValidationRecord("bond", "test", Mock(), Mock(), Mock(), Mock(), 1.60, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55)
    assert record.is_preferred() is False
    assert record.is_allowed() is False
    assert record.is_suspicious() is False
    assert record.is_outlier() is True
