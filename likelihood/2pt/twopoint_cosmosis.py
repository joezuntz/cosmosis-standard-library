"""
This file reads the type_table.txt file, which
tells cosmosis the mapping between the names of 
types of two-point function and the sections in cosmosis data blocks.

This file is used by two others, the save_2pt module and the 2pt_like module
"""

import os
import numpy as np


def load_type_table():
    dirname = os.path.split(__file__)[0]
    table_name = os.path.join(dirname, "type_table.txt")
    table = {}
    for line in open(table_name):
        line=line.strip()
        if line.startswith("#") or line=="":
            continue
        type1, type2, section, x, y = line.split()
        table[(type1,type2)] = (section, x, y)
    return table

type_table = load_type_table()

def theory_names(spectrum):
    section, x_name, y_name = type_table[(spectrum.type1.name,spectrum.type2.name)]
    return section, x_name, y_name
