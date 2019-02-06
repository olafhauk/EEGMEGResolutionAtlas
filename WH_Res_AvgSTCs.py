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
 '/imaging/local/software/mne_python/latest_v0.15',
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
if len(sys.argv)>2: # if additional variable specified

    stc_path = sys.argv[2] # pathname (e.g. 'SensitivityMap')

if len(sys.argv) > 3:

    stc_type = sys.argv[3] # beginning of filename (e.g. 'LocErrPeak')

if len(sys.argv) > 4:

    metric = sys.argv[4] # end of filename (e.g. 'RMS')

    print('###\nChosen variables: %s, %s, %s.\n###' % (stc_path, stc_type, metric))

# Maxfilter parameters for filenames
st_duration = C.res_st_duration
origin = C.res_origin

# create dir for average if necessary
fname_avg = C.fname_STC(C, stc_path, C.stc_morph, '')
if not op.exists(fname_avg):
    os.mkdir(fname_avg)

for modality in ['EEGMEG', 'MEG', 'EEG', 'EEGMEG-MEG', 'EEGMEG-EEG']: # EEG/MEG/EEGMEG

    # contrasts for inverse methods only computed for EEGMEG
    if modality == 'EEGMEG':

        res_inv_types = C.res_inv_types + ['MNE-dSPM', 'MNE-sLORETA', 'dSPM-sLORETA',
                                            'MNE-LCMV', 'dSPM-LCMV', 'sLORETA-LCMV', 'MNE-MNE40', 'MNE-MNE80']

    else:

        res_inv_types = C.res_inv_types

    for inv_type in res_inv_types: # 'MNE', 'LCMV' etc.

        # for CTFs and PSFs
        for functype in ['CTF', 'PSF']:

            # iterate over inverse operator types
                for loose in C.inv_loose: # orientation constraint

                    for depth in C.inv_depth: # depth weighting

                        stcs = [] # Will contain STCs per subject for averaging

                        if loose == None: loose = 0

                        loo_str = '_loo%s' % str(int(100*loose))

                        if depth == None: depth = 0

                        dep_str = '_dep%s' % str(int(100*depth))

                        if inv_type[:7] == 'MNE-MNE': # exception for depth-weighted MNE

                            mytext = functype + '_' + inv_type + '_' + stc_type + '_' + modality + loo_str

                        else:

                            mytext = functype + '_' + inv_type + '_' + stc_type + '_' + modality + loo_str + dep_str

                        if metric != '':

                            mytext = mytext + '_' + metric

                        for sbj in C.subjs:

                            subject = 'Sub%02d' % sbj                

                            fname_morph = C.fname_STC(C, stc_path, subject, mytext + '_mph')

                            # read existing source estimate
                            print('Reading: %s.' % fname_morph)
                            stc = mne.read_source_estimate(fname_morph, subject)

                            stcs.append(stc)

                        # average STCs across subjects
                        print('Averaging %d STC files.' % len(stcs))

                        avg = np.average([s.data for s in stcs], axis=0)

                        # turn average into source estimate object
                        avg_stc = mne.SourceEstimate(avg, stcs[0].vertices, stcs[0].tmin, stcs[0].tstep)

                        fname_avg = C.fname_STC(C, stc_path, C.stc_morph, mytext)

                        print('###\nWriting grand-average STC file %s.\n###' % fname_avg)

                        avg_stc.save(fname_avg)

                        if inv_type[:7] == 'MNE-MNE': # depth-weighted MNE only for one depth

                            break

                    if inv_type[:7] == 'MNE-MNE': # depth-weighted MNE only for one loose

                        break
# Done