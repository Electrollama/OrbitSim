import csv
from Satellite import Body

"""
time units: Days [Earth]
distance units: mega-meters [Mm]
mass units: Mega-grams aka metric ton
mu [m^3/s^2] * 7.4650e-9 = mu [Mm^3/days^2]
"""

def load_system(file_name):
    """
    Read a table of Kepler elements. Compile by initializing objects
    representing the system.
    :param file_name: File path for the table, a txt file
    :return: a Body object representing the parent body
    """
    Table = read_table(file_name)
    # Initialize the system parent
    parent = Table['names'][0] # first row is the parent
    parent_Entry = Table[parent] # dictionary for the parent
    parent_Entry['orbit'] = None
    System = Body(parent_Entry)
    Lookup = {parent: System}
    # Compile the rest of the table
    for body in Table['names'][1:]:
        Entry = Table[body]
        Parent_Body = Lookup[Entry['parent']]
        # replace mu value with that of its parent
        orbit = dict(Entry)
        orbit['mu'] = Parent_Body.mu
        Entry['orbit'] = orbit
        # create body
        Lookup[body] = Body(Entry, Parent_Body)
    # return the system parent body
    return System

### Low-Level ###

Prefixes = {'T': 12, 'G': 9, 'M': 6, 'k': 3,
            'm': -3, 'u': -6, 'n': -9, 'p': -12}

def read_entry(entry):
    result = entry
    if entry[0].isdigit() or (entry[0] == '-' and entry[1].isdigit()):
        if entry[-1] in Prefixes.keys():
            try:
                value = float(entry[:-1])
                power = Prefixes[entry[-1]]
                result = value * 10**power
            except ValueError:
                pass
        else:
            try:
                value = float(entry[0:])
                result = value
            except ValueError:
                pass
    return result

def read_table(file_name=None):
    print('Reading file: {}'.format(file_name))
    if file_name == None:
        file_name = input('File Name:')
    if not file_name[-4:] == '.csv':
        file_name += '.csv'
    Result = {}
    errorcount = 0
    columns = []
    names = []
    #Read celestial bodies file
    with open(file_name, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i == 0:
                columns = list(row)
            else:
                Entries = {}
                try:
                    name = str(row[0])
                    names.append(name)
                    for j, c in enumerate(columns):
                        Entries[c] = read_entry(row[j])
                    Result[name] = Entries
                except (ValueError, IndexError):
                    errorcount += 1
    print('  {} entries found. {} skipped.'.format(len(names), errorcount))
    Result['file_name'] = file_name
    Result['columns'] = columns
    Result['names'] = names
    return Result

def test_load():
    system = load_system('Celestial_Bodies_Demo')
    print(system)
    print(system.satellites['Mun'].orbit)