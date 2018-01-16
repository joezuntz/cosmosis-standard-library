from astropy.table import Table
import os
import numpy as np


def load_type_table():
    dirname = os.path.split(__file__)[0]
    table_name = os.path.join(dirname, "type_table.txt")
    type_table = Table.read(table_name, format="ascii")
    table = {}
    for (type1, type2, section, x, y) in type_table:
        table[(type1, type2)] = (section, x, y)
    return table


type_table = load_type_table()


def theory_names(spectrum):
    section, x_name, y_name = type_table[(
        spectrum.type1.name, spectrum.type2.name)]
    return section, x_name, y_name
