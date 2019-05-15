
#!/imaging/local/software/anaconda/latest/x86_64/bin/python

import sys
import os
import os.path as op
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.image import imread

import glob

import mne
from mne.report import Report

# Directory with subjects as sub-directories
data_dir = '/group/erp/data/olaf.hauk/MEG/WakemanHensonEMEG/Resolution/ResolutionMetrics/'

# sub-directory for subject
subdir = 'Figures'

# list of images files for HTML in data_dir/subject
# without extension
filelist = ['PSF_MNE-dSPM_011_locerr_peak_EEGMEG_loo0_dep0-00000.0-lat-lh']

# plot in interactive mode
plt.ion()

for ff in filelist:

        report = Report(subject=subdir, title=ff)

        # image file to be read
        filename = op.join(data_dir, subdir, ff+'.jpg')

        print('Reading image from %s.' % filename)

        # read image file
        img = imread(filename)

        # subplots creates matplotlib.figure.Figure() instance for report
        fig, ax = plt.subplots(2,2, num=2)

        # plot image
        plt.imshow(img)

        ax[1,1].set_title('Happy')

        # get rid of axis labels and ticks, not required for image
        ax[1,1].set_xticks([])
        ax[1,1].set_yticks([])

        # add image to HTML report
        report.add_figs_to_section(fig, 'Moonshine', section='Sunshine', scale=1)
        
        fname_html_out = op.join(data_dir, subdir, filename+'.html')

        print 'Saving HTML report to {0}'.format(fname_html_out)
        report.save(fname_html_out, overwrite=True)

