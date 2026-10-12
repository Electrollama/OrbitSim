import csv
import os
from Satellite import Body

"""
time units: Days [Earth or Kerbin]
distance units: mega-meters [Mm]
mass units: Mega-grams aka metric ton
RL: mu [10^10 m^3/s^2] * 74.650 = mu [Mm^3/days^2]
KSP: mu [10^10 m^3/s^2] * 4.6656 = mu [Mm^3/days^2]
        or mass[10^20 kg] * 49.829
"""

Prefixes = {
    'T': 12, 'G': 9, 'M': 6, 'k': 3,
    'm': -3, 'u': -6, 'n': -9, 'p': -12
}

def read_entry(entry):
    if not isinstance(entry, str):
        return entry
    entry = entry.strip()
    if not entry:
        return entry
    if entry[0].isdigit() or (entry[0] == '-' and len(entry) > 1 and entry[1].isdigit()):
        if entry[-1] in Prefixes:
            try:
                value = float(entry[:-1])
                power = Prefixes[entry[-1]]
                return value * (10 ** power)
            except ValueError:
                pass
        else:
            try:
                return float(entry)
            except ValueError:
                pass
    return entry


def read_table(file_name=None):
    if file_name is None:
        file_name = input('File Name: ')
    
    candidates = [
        file_name,
        file_name + '.csv',
        os.path.join('System Tables', file_name),
        os.path.join('System Tables', file_name + '.csv'),
        os.path.join(os.path.dirname(__file__), file_name),
        os.path.join(os.path.dirname(__file__), file_name + '.csv'),
        os.path.join(os.path.dirname(__file__), 'System Tables', file_name),
        os.path.join(os.path.dirname(__file__), 'System Tables', file_name + '.csv'),
    ]

    resolved_path = None
    for cand in candidates:
        if os.path.isfile(cand):
            resolved_path = cand
            break

    if resolved_path is None:
        resolved_path = file_name if file_name.endswith('.csv') else file_name + '.csv'

    Result = {}
    errorcount = 0
    columns = []
    names = []

    with open(resolved_path, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if not row or not any(row):
                continue
            if i == 0 or not columns:
                columns = [c.strip() for c in row]
            else:
                Entries = {}
                try:
                    name = str(row[0]).strip()
                    names.append(name)
                    for j, c in enumerate(columns):
                        if j < len(row):
                            Entries[c] = read_entry(row[j])
                        else:
                            Entries[c] = None
                    Result[name] = Entries
                except (ValueError, IndexError):
                    errorcount += 1

    Result['file_name'] = resolved_path
    Result['columns'] = columns
    Result['names'] = names
    return Result


def load_system(file_name):
    """
    Read a table of Kepler elements. Compile by initializing objects
    representing the system.
    :param file_name: File path for the table
    :return: a Body object representing the parent body
    """
    Table = read_table(file_name)
    parent_name = Table['names'][0]
    parent_Entry = dict(Table[parent_name])
    parent_Entry['orbit'] = None
    System = Body(parent_Entry)
    Lookup = {parent_name: System}

    for body_name in Table['names'][1:]:
        Entry = dict(Table[body_name])
        Parent_Body = Lookup[Entry['parent']]
        orbit = dict(Entry)
        orbit['mu'] = Parent_Body.mu
        Entry['orbit'] = orbit
        Lookup[body_name] = Body(Entry, Parent_Body)

    return System
