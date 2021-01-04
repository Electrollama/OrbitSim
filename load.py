import csv
from Satellite import Body
from tkinter import Tk, filedialog
from os import listdir, path, makedirs, getcwd
from numbers import Number

"""
time units: Days [Earth or Kerbin]
distance units: mega-meters [Mm]
mass units: Mega-grams aka metric ton
RL: mu [10^10 m^3/s^2] * 74.650 = mu [Mm^3/days^2]
KSP: mu [10^10 m^3/s^2] * 4.6656 = mu [Mm^3/days^2]
		or mass[10^20 kg] * 49.829
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
    parent_name = Table['names'][0] # first row is the parent
    parent_Entry = Table[parent_name] # dictionary for the parent
    parent_Entry['orbit'] = None
    System = Body(parent_Entry)
    Lookup = {parent_name: System}
    # Compile the rest of the table
    for body_name in Table['names'][1:]:
        Entry = Table[body_name]
        Parent_Body = Lookup[Entry['parent']]
        # replace mu value with that of its parent
        orbit = dict(Entry)
        orbit['mu'] = Parent_Body.mu
        Entry['orbit'] = orbit # copy of the dictionary
        # create body
        Lookup[body_name] = Body(Entry, Parent_Body)
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
    print(system.bodies['Mun'].orbit)


def write_file(Data, filepath=None, delim='\t'):
    """

    :param Data:
    :param folder:
    :param header:
    :param delim:
    :param save_clean: If false, save similar to the raw data file. If true, save a cleaner version
    :return:
    """
    if 'columns' not in Data.keys():
        raise KeyError("Data dictionary needs key: 'columns' to write a file.")
    if 'header' not in Data.keys():
        Data['header'] = []
    if filepath is None:
        filepath = file_dialog('save')
    print('Writing file:', filepath)
    with open(filepath, 'w', newline='') as f:
        data_writer = csv.writer(f, delimiter=delim)
        for h in Data['header']:
            data_writer.writerow([h, Data[h]])  # includes columns
        data_writer.writerow(Data['columns'])
        n_rows = max([len(Data[col]) for col in Data['columns']])
        for i in range(n_rows):  # rows
            row = []
            for col in Data['columns']:  # columns
                if i < len(Data[col]):
                    entry = Data[col][i]
                    # if isinstance(entry, Number):
                    #    entry = format(entry, '.6g')
                    row.append(entry)  # value
                else:
                    row.append('NaN')
            data_writer.writerow(row)
    print('Writing Done.')
    return


def file_dialog(operation='file', types=(('CSV', '.csv'),)):
    root = Tk()
    if operation == 'folder':
        folder_path = filedialog.askdirectory(parent=root,
                                              initialdir=getcwd())
        root.destroy()
        files = listdir(folder_path)
        files = [folder_path + '/' + f for f in files]
        folder_name = path.basename(folder_path)
        print("Opened {} files from {}".format(len(files), folder_name))
        print("   path: {}".format(folder_path))
        return files, folder_name
    elif operation == 'file':
        file_path = filedialog.askopenfilename(filetypes=types,
                                               initialdir=getcwd())
        root.destroy()
        print("Opened: {}".format(file_path))
        return file_path
    elif operation in ['save', 'write']:
        file_path = filedialog.asksaveasfile(filetypes=types,
                                             initialdir=getcwd()).name
        root.destroy()
        print("Created: {}".format(file_path))
        return file_path
    else:
        raise ValueError("'{}' is not a defined operation.".format(operation))