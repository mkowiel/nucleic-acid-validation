# nucleic-acid-validation

Tool for validation of RNA/DNA bonds and angles geometry.

# Installation

At the moment there are 2 ways to install the library, install from github or from local version of the code


## Install from github

    # set up the virtual environment (python3 -m venv <envdir>)
    python3 -m venv navalenv
    # activate the virtual environment (source <envdir>/bin/activate)
    source navalenv/bin/activate
    # install library directly form github (the latest version from main branch)
    python -m pip install git+https://github.com/mkowiel/nucleic-acid-validation.git

or install the stable version (for example 0.0.6)

    python -m pip install git+https://github.com/mkowiel/nucleic-acid-validation.git@0.0.6

## Install from clone

    git clone https://github.com/mkowiel/nucleic-acid-validation.git
    cd nucleic-acid-validation
    python -m pip install .

# Execute

(optionally) activate the virtual environment (source <envdir>/bin/activate)

    source navalenv/bin/activate

Run the validation tool

    naval <struct.cif|struct.pdb> <bonds.csv> <angles.csv> <geometry.csv>

For example for 3p4j structure

    naval 3p4j.cif 3p4j_bonds.csv 3p4j_angles.csv 3p4j_geometry.csv

# Output format

The validation of nucleotides bonds and angles are stored in a `.csv` format.

- **type**: `bond`, `angle`, `torsion` or `pseudorotation`
- **pdbcode**: pdbcode extracted from the file name
- **model_id**: model id
- **chain**: chain name
- **atom1_res_name**: residue name
- **atom1_resid**: residue id
- **atom1_name**: atom label (for example `N1`, `C4'`,`P`, `OP1`)
- **atom1_altloc**: alternative conformation
- **atom2_res_name**: residue name
- **atom2_resid**: residue id
- **atom2_name**: atom label
- **atom2_altloc**: alternative conformation
- **atom3_res_name**: residue name  (for angles only)
- **atom3_resid**: residue id  (for angles only)
- **atom3_name**: atom label  (for angles only)
- **atom3_altloc**: alternative conformation (for angles only)
- **calculated:** calculate value for bond length, angle or torsion angle for given atoms
- **target**: expected bond length, angle or torsion angle for given atoms (based on the CSD)
- **validation_label**: a 4-tier validation category (`CSD-preferred`, `PDB-acceptable`, `PDB-suspicious` or `PDB-outlier`)
- **validator_name**: internal validator name, based on detected conformation (for debug purposes)

# Description

The library uses Biopython library to parse a structure in a mmCif or pdb format. Then calculates nucleotide geometry (torsion angles and sugar pucker pseudorotation) and validates bonds and angles based on the CSD-derived target values and PDB-derived distributions for conformation dependent group. CSD-based target values are based on the publications:

M.Kowiel, D.Brzezinski, M.Jaskolski (2016).
*Conformation-dependent restraints for polynucleotides: I. clustering of the geometry of the phosphodiester group.*
Nucleic Acids Res. 44, 8479–8489. https://doi.org/10.1093/nar/gkw717 OpenAccess

M.Gilski, J.Zhao, M.Kowiel, D.Brzezinski, D.H.Turner, M.Jaskolski (2019).
*Accurate geometrical restraints for Watson–Crick base pairs.*
Acta Cryst. B75, 235-245. https://doi.org/10.1107/S2052520619002002 OpenAccess

M.Kowiel, D.Brzezinski, M. Gilski, M.Jaskolski (2020).
*Conformation-dependent restraints for polynucleotides: The sugar moiety.*
Nucleic Acids Res. 48, 962–973. https://doi.org/10.1093/nar/gkz1122 OpenAccess

# Run unit tests

Tests are configured with `tox` library. For more detail check `Makefile`. The test pipeline should work
on linux environment.

To execute tests suit run:

    make test

If the make is missing you can execute the test manually

    # (activate the virtual environment)
    # make sure the python has the required libraries
    python -m pip install --upgrade tox
    # run tests via tox
    python -m tox

# Fix linting (for developers)

    make black isort

# Changelog

## 0.0.6

- script file
- improved readme

## 0.0.5

- Bonds and angle stats based on structures PDB/[NAKB](https://nakb.org/) (Release date >= 2010-01-01, Resolution <= 3.6A, Method = X-Ray, Without alternative conformation)
- Neighbors atom search bug fixes

## 0.0.4

- Sugar validation based on residue name and sugar pucker (C2'-endo, C3'-endo, other)
- 4 tier naming change: CSD-preferred, PDB-acceptable, PDB-suspicious, PDB-outlier
- Separate output file for bonds, angles, and torsion validation records.
- Improvements in README.md

## 0.0.3

- Validate Base based on residue name, PO4 based on alpha and zeta torsion angles, and sugar part based only
on residue name.

# TODO:

- add options to select sugar validation level (residue, pseudorotation, conformation dependent)
- do not copy the target values to validation records, we can reuse target definition (may be problematic for functional dependencies)
- more tests
- handle terminal sugars
- handle terminal PO4
