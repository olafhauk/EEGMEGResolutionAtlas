"""
=========================================================
Grand-average of morphed STCs for resolution metrics
for WH data set.
Doesn't run in parallel mode.
e.g.: run WH_Res_AvgSTC.py WH_MNE_Resolution_config SensitivityMaps SensMap RMS
or: run WH_Res_AvgSTC.py WH_MNE_Resolution_config ResolutionMetrics locerr_peak
=========================================================

"""
# OH, July 2018

print __doc__

import os
import os.path as op

import sys
sys.path = [
 '/home/olaf/MEG/WakemanHensonEMEG/ScriptsResolution', # following list created by trial and error
 '/imaging/local/software/mne_python/latest_v0.16',
 '/imaging/local/software/anaconda/2.4.1/2/bin',
 '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/',
 '/imaging/local/software/anaconda/2.4.1/2/envs/mayavi_env/lib/python2.7/site-packages',
 '/imaging/local/software/anaconda/2.4.1/2/envs/mayavi_env/lib/python2.7/site-packages/pysurfer-0.8.dev0-py2.7.egg',
 '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/h5io-0.1.dev0-py2.7.egg',
 '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/lib-dynload',
 '/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages'
 ]

import importlib
import glob

import numpy as np

import mne
print('MNE Version: %s\n\n' % mne.__version__) # just in case

## get analysis parameters from config file

module_name = sys.argv[1]

C = importlib.import_module(module_name)
reload(C)

# get functions for metrics etc.
R = importlib.import_module('WH_Resolution_Functions')
reload(R)

# Type of inverse method, 'norm' | 'lcmv'
inv_types = ['norm', 'lcmv']

# hack to have variables via qsub
stc_path, stc_type, metric = '', '', ''

# read variables specified via qsub
if len(sys.argv)>3: # if additional variable specified

    metric = sys.argv[3] # (e.g. 'RMS', 'SNR')

else:

    metric = 'RMS'


# Maxfilter parameters for filenames
st_duration = C.res_st_duration
origin = C.res_origin

# create dir for average if necessary
fname_avg = C.fname_STC(C, 'SensitivityMaps', C.stc_morph, '')
if not op.exists(fname_avg):
    os.mkdir(fname_avg)

for modality in ['EEGMEG', 'MEG', 'EEG', 'EEGMEG-MEG', 'EEGMEG-EEG']: # EEG/MEG/EEGMEG

    stcs = []

    for sbj in C.subjs:

        subject = 'Sub%02d' % sbj                

        fname_morph = C.fname_STC(C, 'SensitivityMaps', subject, 'SensMap_' + modality + '_' + metric + '_mph')

        # read existing source estimate
        print('Reading: %s.' % fname_morph)
        stc = mne.read_source_estimate(fname_morph, subject)

        stcs.append(stc)

    # average STCs across subjects
    print('Averaging %d STC files.' % len(stcs))

    avg = np.average([s.data for s in stcs], axis=0)

    # turn average into source estimate object
    avg_stc = mne.SourceEstimate(avg, stcs[0].vertices, stcs[0].tmin, stcs[0].tstep)

    fname_avg = C.fname_STC(C, 'SensitivityMaps', C.stc_morph, 'SensMap_' + modality + '_' + metric)

    print('###\nWriting grand-average STC file %s.\n###' % fname_avg)

    avg_stc.save(fname_avg)
# Done