# nucleic-acid-validation

Tool for validation of RNA/DNA bonds and angles geometry.

# Description

The library uses Biopython library to parse a pdb structure.

Calulate nucleotide gemoetry (torson angles and sugar pseudorotation) and validates bonds and angles agains
the CSD-derived and PDB-derived target values for given geometry. CSD-based target values based on publications:

M.Kowiel, D.Brzezinski, M.Jaskolski (2016).
*Conformation-dependent restraints for polynucleotides: I. clustering of the geometry of the phosphodiester group.*
Nucleic Acids Res. 44, 8479–8489. https://doi.org/10.1093/nar/gkw717 OpenAccess

M.Gilski, J.Zhao, M.Kowiel, D.Brzezinski, D.H.Turner, M.Jaskolski (2019).
*Accurate geometrical restraints for Watson–Crick base pairs.*
Acta Cryst. B75, 235-245. https://doi.org/10.1107/S2052520619002002 OpenAccess

M.Kowiel, D.Brzezinski, M. Gilski, M.Jaskolski (2020).
*Conformation-dependent restraints for polynucleotides: The sugar moiety.*
Nucleic Acids Res. 48, 962–973. https://doi.org/10.1093/nar/gkz1122 OpenAccess


# Tests

Tests are configured with `tox` library. For more detail check `Makefile`. The test pipeline should work
on linux environment.

To execute tests suit run:

    make test

# Changelog

## 0.0.4

- 4 tier naming change: CSD-preferred, PDB-acceptable, PDB-suspicious, PDB-outlier
- Separate output file for bonds, angles, and torsion validation records.
- Improvements in README.md

## 0.0.3

- Validate Base based on residue name, PO4 based on alhpa and zeta torsion angles, and sugar part based only
on residue name.

# TODO:

- add sugar conformation depenent values
- add options to select sugar validation level (residue, pseudorotation, conformation dependent)
- do not copy the target values to validation records, we can reuse target definition (may be problematic for functional dependencies)
- more tests

