from unittest.mock import Mock

from naval.validation_record import ValidationRecord
from naval.printer import CsvPrinter


def test_csv_printer():
    records = [
        ValidationRecord("bond", Mock(), Mock(), Mock(), Mock(), 1.43, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55),
        ValidationRecord("angle", Mock(), Mock(), Mock(), Mock(), 1.55, 1.40, 0.02, 1.30, 1.50, 1.25, 1.55),
    ]

    printer = CsvPrinter()
    lines = printer.print(records)
    assert len(lines) == 3
    assert "bond" in lines[1]
    assert "angle" in lines[2]
