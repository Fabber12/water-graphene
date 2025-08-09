#!/usr/bin/env python3

import sys
import math

"""
Created by F. Tarulli, Politecnico di Torino, Italy - February 27, 2025

    Remove lines with atom types 1 or 2 (water molecules) and update the atom count
    based on the given oxidation percentage.

"""
def round_as_matlab(x):
    if x - math.floor(x) == 0.5:
        return math.ceil(x)
    else:
        return round(x)

if len(sys.argv) != 4:
    print("Usage: python remove_dump_lines.py input.dump output.dump oxidation_percent")
    sys.exit(1)

input_file, output_file, oxidation = sys.argv[1], sys.argv[2], sys.argv[3]

OH_groups = round_as_matlab(1972 * int(oxidation) / 100)
atoms = str(1972 + 2 * OH_groups)
replace_next = False

print(f"Total atoms: {atoms} (C from graphene + O atoms + H from OH groups)")

with open(input_file, "r") as fin, open(output_file, "w") as fout:
    for line in fin:
        if replace_next:
            fout.write(atoms + "\n")
            replace_next = False
            continue

        if line.strip() == "ITEM: NUMBER OF ATOMS":
            fout.write(line)
            replace_next = True
            continue

        if line and line[0].isdigit():
            tokens = line.split(maxsplit=2)
            if len(tokens) > 1 and tokens[1] in ("1", "2"):
                continue

        fout.write(line)

