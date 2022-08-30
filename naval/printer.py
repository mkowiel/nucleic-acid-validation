from typing import List

from naval.validation_record import ValidationRecord


class Printer:
    """
    Base class for all Printers.
    """

    @classmethod
    def format_header(cls):
        return ""

    @classmethod
    def format_record(cls, record: ValidationRecord):
        return record.validation_type

    @classmethod
    def print(cls, records: List[ValidationRecord]):
        lines = []
        header = cls.format_header()
        if header:
            lines.append(header)
        for record in records:
            lines.append(cls.format_record(record))
        return lines


class CsvPrinter(Printer):
    """
    CSV printer converts Validation Records to lines of text.
    """

    @classmethod
    def format_header(cls):
        return (
            "type,pdbcode,modeid,chain,res_name,resseq,inscode,atom1_name,atom1_altloc,"
            "atom2_name,atom3_altloc,atom1_name,atom3_altloc,calculated,target,preferred,allowed,suspicious,outlier"
        )

    @classmethod
    def format_record(cls, record: ValidationRecord):
        line = ",".join(
            str(_)
            for _ in (
                record.validation_type,
                record.geometry.pdbcode,
                record.geometry.model.get_id(),
                record.geometry.chain.get_id(),
                record.geometry.res_name,
                record.geometry.resseq,
                record.geometry.inscode,
                record.atom1.get_name(),
                record.atom1.get_altloc(),
                record.atom2.get_name(),
                record.atom2.get_altloc(),
                record.atom3.get_name() if record.atom3 else "",
                record.atom3.get_altloc() if record.atom3 else "",
                record.calculated_value,
                record.target_value,
                record.preferred,
                record.allowed,
                record.suspicious,
                record.outlier,
            )
        )
        return line
