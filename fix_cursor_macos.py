# This script fixes the issue the latest MacOS
# have with the iraf display: the cursor freezes
import os
import wget
import pyraf
import shutil

fixed_file_url = 'https://github.com/svalenti/pessto/blob/python3/Ptkplot.py'
wget.download(fixed_file_url)

pyraf_path = pyraf.__path__[0]
ptkplot_path = os.path.join(pyraf_path, 'Ptkplot.py')
shutil.move('Ptkplot.py', ptkplot_path)
