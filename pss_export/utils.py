'''
Utility functions for PSS exports.
'''

import csv
import logging

logger = logging.getLogger(__name__)


def sbo(term):
    ''' 'SBO:0000177' from 177 or '0000177'; None stays None '''
    return None if term is None else f"SBO:{int(term):07d}"


# lists are joined with '|' (as TAIR's files and GAF): names contain commas (AHK2,3,4) and gene symbols ';' (PIP1;3)
LIST_SEPARATOR = '|'


def clean(value, where=''):
    ''' A value for a tab-separated file: whitespace collapsed, lists joined with '|'. The database's list
    values have no '|' (checked when they are entered, and by a query in regular-updates.cypher); one that
    does anyway is written with '/' instead, with a warning, so that the list can still be split. '''
    if value is None:
        return ''
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (list, tuple)):
        values = [clean(v) for v in value]
        for v in values:
            if LIST_SEPARATOR in v:
                logger.warning(f"'{LIST_SEPARATOR}' in a list value, written with '/' ({where}): {v}")
        return LIST_SEPARATOR.join(v.replace(LIST_SEPARATOR, '/') for v in values)
    return ' '.join(str(value).split())


def write_tsv(filename, columns, rows):
    ''' Write rows (dicts) as a tab-separated file with a header '''
    with open(filename, 'w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\', lineterminator='\n')
        writer.writerow(columns)
        writer.writerows([clean(row.get(c), f"{row.get('id', row.get('reaction_id'))}, {c}") for c in columns]
                         for row in rows)


# the model-level note of a model "with model fixes" (SBML, TabularQual)
MODEL_FIXES_NOTE = ("{n} applied by pss-export (translation products in the active form, transport reactions "
                    "between compartments): this model has reactions and changes that are not in PSS (marked in "
                    "their notes)")
