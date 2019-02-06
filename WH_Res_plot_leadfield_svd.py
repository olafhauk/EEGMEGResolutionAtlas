"""
=========================================================
Plot singular values for whitened leadfields
from WH_Res_LeadfieldSVD.py
E.g.: run WH_Res_plot_leadfield.py WH_MNE_Resolution_config
=========================================================

"""
# OH, September 2018

print __doc__

import sys
import importlib
import numpy as np

from matplotlib import pyplot as plt

module_name = sys.argv[1]

C = importlib.import_module(module_name)
reload(C)

plt.ion()

# maximum index for plotting
plot_to = 200

legend = []
lines = []

fig, ax = plt.subplots(1,1)

cols = ['r', 'b', 'g']

for [ii,modal] in enumerate(['EEGMEG', 'MEG', 'EEG']):

    fname_txt = C.fname_STC(C, 'SensitivityMaps', '', modal + '_SVD.txt')

    print('Loading singular values from: %s' % fname_txt)
    
    sing_vals = np.loadtxt(fname_txt)

    avg_vals = sing_vals[:,:plot_to].mean(axis=0)

    line, = ax.plot(avg_vals, cols[ii], linewidth=2)

    lines.append(line)

    legend.append(modal)

    ax.xaxis.set_tick_params(labelsize=24)
    ax.yaxis.set_tick_params(labelsize=24)

fig.legend(lines, legend, 'upper right')

plt.tight_layout()
