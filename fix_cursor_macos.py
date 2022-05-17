import os
import shutil
import requests
# importing pyraf might fail due to a pickle protocol issue
# this is solved with fix_pickle_macos.py
import pyraf
pyraf_path = pyraf.__path__[0]

# ---------------------------
# this fixes the cursor issue
fixed_file_url = 'https://raw.githubusercontent.com/svalenti/pessto/python3/Ptkplot.py'
data = requests.get(fixed_file_url)

with open('Ptkplot.py', 'wb')as file:
    file.write(data.content)

# update the file with the solution
ptkplot_path = os.path.join(pyraf_path, 'Ptkplot.py')
shutil.move('Ptkplot.py', ptkplot_path)
