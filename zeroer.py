import numpy as np
import os
import matplotlib.pyplot as plt

import adminfunctions as ad

directory = os.getcwd() + '/spectra_not_zeroed'

fileno = ad.listfiles(directory)
fileselect = ad.inputfileselector(fileno)
f = ad.openfilefromlist(fileselect, directory)

xarr, yarr = ad.file_reader(directory + '/' + f, delete_0 = False)

xmax = max(xarr)

newx = np.arange(xmax) + 0.5
newy = []
xprev = None
i = 0

for x, y in zip(xarr, yarr):
    if (x == xarr[i - 1] + 1) | (x == 0.5):
        newy.append(int(y))
    else:
        diff = int(x - xarr[i - 1])
        for j in range(diff):
            newy.append(0)
    i = i + 1 

os.chdir(os.getcwd() + '/spectra')

with open(f + '_zeroed', 'w+') as f2:
    for x, y in zip(newx, newy):
        f2.write(str(x) + ' ' + str(y) + '\n')
    
    

   
