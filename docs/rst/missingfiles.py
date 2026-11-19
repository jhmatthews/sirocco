#%%

from astropy.table import Table
import numpy

t = Table.read('param_table.txt', format='ascii')
filename = t['File']
path_1 = t['Path_1']
path_2 = t['Path_2']

for i in range(len(filename)):
    # i
    if not isinstance(path_1[i], numpy.str_):
        print('Add to Docs: {}'.format(filename[i]))
    if not isinstance(path_2[i], numpy.str_):
        print('Remove from Docs: {}'.format(filename[i]))
# %%
