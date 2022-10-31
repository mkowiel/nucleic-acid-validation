from typing import List, Tuple

from naval.validation_record import TorsionRecord, ValidationRecord


class Printer:
    """
    Base class for all Printers.
    """

    supported_record_types: Tuple[str] = tuple()  # type: ignore

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
            if record.validation_type in cls.supported_record_types:
                lines.append(cls.format_record(record))
        return lines


class BondsCsvPrinter(Printer):
    """
    CSV printer converts Validation Records to lines of text.
    """

    supported_record_types: Tuple[str] = ("bond",)

    @classmethod
    def format_header(cls):
        return (
            "type,pdbcode,model_id,chain,"
            "atom1_res_name,atom1_resid,atom1_name,atom1_altloc,"
            "atom2_res_name,atom2_resid,atom2_name,atom2_altloc,"
            "calculated,target,validation_label"
        )

    @classmethod
    def format_record(cls, record: ValidationRecord):
        line = ",".join(
            str(_)
            for _ in (
                record.validation_type,
                record.geometry.residue_entry.pdbcode,
                record.geometry.residue_entry.model.get_id(),
                record.geometry.residue_entry.chain.get_id(),
                record.atom1.get_parent().get_resname(),
                str(record.atom1.get_parent().get_id()[1]) + record.atom1.get_parent().get_id()[2].strip(),
                record.atom1.get_name(),
                record.atom1.get_altloc().strip(),
                record.atom2.get_parent().get_resname(),
                str(record.atom2.get_parent().get_id()[1]) + record.atom2.get_parent().get_id()[2].strip(),
                record.atom2.get_name(),
                record.atom2.get_altloc().strip(),
                round(record.calculated_value, 3),
                round(record.target_value, 3),
                record.label,
            )
        )
        return line


class AnglesCsvPrinter(Printer):
    """
    CSV printer converts Validation Records to lines of text.
    """

    supported_record_types: Tuple[str] = ("angle",)

    @classmethod
    def format_header(cls):
        return (
            "type,pdbcode,model_id,chain,"
            "atom1_res_name,atom1_resid,atom1_name,atom1_altloc,"
            "atom2_res_name,atom2_resid,atom2_name,atom2_altloc,"
            "atom3_res_name,atom3_resid,atom3_name,atom3_altloc,"
            "calculated,target,validation_label"
        )

    @classmethod
    def format_record(cls, record: ValidationRecord):
        line = ",".join(
            str(_)
            for _ in (
                record.validation_type,
                record.geometry.residue_entry.pdbcode,
                record.geometry.residue_entry.model.get_id(),
                record.geometry.residue_entry.chain.get_id(),
                record.atom1.get_parent().get_resname(),
                str(record.atom1.get_parent().get_id()[1]) + record.atom1.get_parent().get_id()[2].strip(),
                record.atom1.get_name(),
                record.atom1.get_altloc().strip(),
                record.atom2.get_parent().get_resname(),
                str(record.atom2.get_parent().get_id()[1]) + record.atom2.get_parent().get_id()[2].strip(),
                record.atom2.get_name(),
                record.atom2.get_altloc().strip(),
                record.atom3.get_parent().get_resname() if record.atom3 else "",
                str(record.atom3.get_parent().get_id()[1]) + record.atom3.get_parent().get_id()[2].strip()
                if record.atom3
                else "",
                record.atom3.get_name() if record.atom3 else "",
                record.atom3.get_altloc().strip() if record.atom3 else "",
                round(record.calculated_value, 1),
                round(record.target_value, 1),
                record.label,
            )
        )
        return line


class GeometryCsvPrinter:
    """
    CSV printer converts Torsion Records to lines of text.
    """

    @classmethod
    def format_header(cls):
        return "type,pdbcode,model_id,chain,res_name,resid,altloc,name,calculated,validation_label"

    @classmethod
    def format_record(cls, record: TorsionRecord):
        line = ",".join(
            str(_)
            for _ in (
                record.validation_type,
                record.geometry.residue_entry.pdbcode,
                record.geometry.residue_entry.model.get_id(),
                record.geometry.residue_entry.chain.get_id(),
                record.geometry.residue_entry.res_name,
                str(record.geometry.residue_entry.resseq) + record.geometry.residue_entry.inscode.strip(),
                record.alt_loc.strip(),
                record.name,
                round(record.calculated_value, 1) if record.calculated_value else "",
                record.calculated_value_label,
            )
        )
        return line

    @classmethod
    def print(cls, records: List[TorsionRecord]):
        lines = []
        header = cls.format_header()
        if header:
            lines.append(header)
        for record in records:
            lines.append(cls.format_record(record))
        return lines
