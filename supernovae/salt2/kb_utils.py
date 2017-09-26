from builtins import range
from collections import OrderedDict
import numpy as np

INT_META = ['NFAKE', 'CHIP', 'NFAKE_FOUND', 'NFAKE_ACCEPT']
STR_META = ['BAND']

INT_COLS = ['NID', 'EXPNUM', 'CHIP', 'FAKEID', 'GALID', 'SAMEHOST', 'OBJID',
            'REJECT', 'FILTER', 'NPOS0', 'NPOS1', 'NPOS2', 'NNEG0', 'NNEG1',
            'NNEG2',  'REJECT', 'ERR', 'TILE', 'NSTAR']
STR_COLS = ['BAND', 'CID']


def dict_to_ndarray(d):
    """Convert a dictionary of lists (of equal length) to a structured
    numpy.ndarray"""

    # first convert all lists to 1-d arrays, in order to let numpy
    # figure out the necessary size of the string arrays.
    for key in d:
        d[key] = np.array(d[key])

    # Determine dtype of output array.
    dtypelist = []
    for key in d:
        dtypelist.append((key, d[key].dtype))

    # Initialize ndarray and then fill it.
    firstkey = list(d.keys())[0]
    col_len = len(d[firstkey])
    result = np.empty(col_len, dtype=dtypelist)
    for key in d:
        result[key] = d[key]

    return result


def read_datafile(filename, default_tablename=None, output_ndarray=False):
    """Read a data file with this (newline-independent) format:

    * Header with 'KEYWORD: value' pairs (optional) 
    * One or more tables. The start of a new table is indicated by a
      keyword that starts with 'NVAR'. If the format is
      'NVAR_TABLENAME', then 'TABLENAME' is taken to be the name of the table,
      and datalines in the table must be started with 'TABLENAME:'
    * A keyword starting with 'VARNAMES' must directly follow the
      'NVAR_TABLENAME' definition. 

    Parameters
    ----------
    filename : str
        Filename of object to read.
    default_tablename : str

    Returns
    -------
    meta : OrderedDict
        Metadata.
    tables : dictionary of multiple `numpy.ndarray`s or dictionaries.
        The data.
    """

    meta = OrderedDict()  # initialize structure to hold metadata.
    tables = {}  # initialize structure to hold data.

    infile = open(filename, 'r')
    words = infile.read().split()
    infile.close()
    i = 0
    nvar = None
    tablename = None
    while i < len(words):
        word = words[i]

        # If the word starts with 'NVAR', we are starting a new table.
        if word.startswith('NVAR'):
            nvar = int(words[i + 1])

            # Infer table name. The name will be used to designate a data row.
            if '_' in word:
                pos = word.find('_') + 1
                tablename = word[pos:].rstrip(':')
            elif default_tablename is None:
                raise ValueError('table name not given as part of NVAR keyword.'
                                 ' You must give a `default_tablename` so that'
                                 ' data rows can be identified in the file.')
            else:
                tablename = default_tablename

            table = OrderedDict()
            tables[tablename] = table

            i += 2

        # If the word starts with 'VARNAMES', the following `nvar` words
        # define the column names of the table.
        elif word.startswith('VARNAMES'):

            # Check that nvar is defined and that no column names are defined
            # for the current table.
            if nvar is None or len(table) > 0:
                raise Exception('NVAR must directly precede VARNAMES')

            # Read the column names
            for j in range(i + 1, i + 1 + nvar):
                table[words[j]] = []
            i += nvar + 1

        # If the word matches the current tablename, we are reading a data row.
        elif word.rstrip(':') == tablename:
            for j, colname in enumerate(table.keys()):
                table[colname].append(words[i + 1 + j])
            i += nvar + 1

        # Otherwise, we are reading metadata or some comment
        # If the word ends with ":", it is metadata.
        elif word[-1] == ':':
            name = word[:-1]  # strip off the ':'
            if len(words) >= i + 2:
                meta[name] = words[i + 1]
            else:
                meta[name] = None
            i += 2
        else:
            # It is some comment; continue onto next word.
            i += 1

    # Convert Values in metadata
    for key, val in list(meta.items()):
        if key in INT_META:
            meta[key] = int(val)
        elif key in STR_META:
            pass
        else:
            meta[key] = float(val)

    # Convert values in columns
    for tablename, table in list(tables.items()):
        for key in list(table.keys()):
            if key in INT_COLS:
                table[key] = [int(float(x)) for x in table[key]]  # HACK
            elif key in STR_COLS:
                pass
            else:
                table[key] = [float(x) for x in table[key]]

    # Convert tables to ndarrays
    if output_ndarray:
        for tablename in list(tables.keys()):
            tables[tablename] = dict_to_ndarray(tables[tablename])

    return meta, tables
