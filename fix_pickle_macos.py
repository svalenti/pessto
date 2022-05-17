import os
import importlib

init_file = importlib.util.find_spec("pyraf")
pyraf_path = os.path.dirname(init_file.origin)

# -------------------
# this fixes the pickle protocol issue
file2fix = os.path.join(pyraf_path, 'sqliteshelve.py')
# read the file and add the fix
with open(file2fix, "rt") as file:
    data = file.read()
    #data = data.replace('pickle.HIGHEST_PROTOCOL', 'protocol=4')
    data = data.replace('import pickle', 'import pickle5 as pickle')

# open the input file in write mode to actually fix it
with open(file2fix, "wt") as file:
    file.write(data)
