"""
=========================================================
Compute singular values for whitened leadfields, for all subjects.
E.g.: run WH_Res_LeadfieldSVD.py WH_MNE_Resolution_config
=========================================================

"""
# OH, September 2018

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

# subjects to process
sbj_ids = C.subjs

# for filenames
st_duration = C.res_st_duration
origin = C.res_origin

# dict for modalities, will contain lists across subjects
sing_vals = {'EEGMEG': [], 'MEG': [], 'EEG': []}

for sbj in sbj_ids:

    subject = 'Sub%02d' % sbj

    print('###\nWorking hard on %s.\n###' % (subject))

    fwd_fname = C.fname_ForwardSolution(C, subject, 'EEGMEG')

    print('###\nReading EEGMEG forward solutions: %s .\n###' % (fwd_fname))

    fwd = mne.read_forward_solution(fwd_fname)

    fwd = mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True)

    # covariance matrix (filter with wildcard)
    fname_cov = C.fname_cov(C, subject, st_duration, origin, C.res_cov_latwin, C.inv_method, '*')

    # method may be underspecified, since it may be ranked differently for different subjects
    fname_cov = glob.glob(fname_cov)[0] # be careful if multiple options present

    print('###\nReading noise covariance matrix from: %s.\n###' % (fname_cov))

    noise_cov = mne.read_cov(fname_cov)

    # iterate over different combinations of sensors
    for (eeg,meg,modal) in [(True,True,'EEGMEG'), (False,True,'MEG'), (True,False,'EEG')]:

        fwd_use = mne.pick_types_forward(fwd, meg=meg, eeg=eeg)

        info = fwd_use['info']
        ch_names = info['ch_names']

        # restrict to channels in forward solution
        noise_cov_use = mne.cov.pick_channels_cov(noise_cov, ch_names)

        # fwd doesn't have projs, add from noise_cov
        info['projs'] = noise_cov_use['projs']
        info['comps'] = '' # dummy to avoid crash

        if C.inv_method == 'empirical': # if unregularised

            noise_cov_use = mne.cov.regularize(noise_cov_use, info, mag=C.res_lambda_empirical['mag'],
                                                    grad=C.res_lambda_empirical['grad'], eeg=C.res_lambda_empirical['eeg'])
        
        # SVD of (diagonally) whitened leadfield
        s = R.leadfield_svd(fwd_use, noise_cov_use, diag=True)        

        sing_vals[modal].append(s)

for modal in ['EEGMEG', 'MEG', 'EEG']:

    fname_txt = C.fname_STC(C, 'SensitivityMaps', '', modal + '_SVD.txt')

    print('Saving singular values to: %s' % fname_txt)
    np.savetxt(fname_txt, np.array(sing_vals[modal]))

# Done