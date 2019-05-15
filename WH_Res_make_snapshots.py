"""
=========================================================
Make snapshots from STC files using mne_make_movie
Set up MNE environment in linux shell first
=========================================================

"""
# OH, May 2019

print __doc__

import os
from os import path as op

outpath = '/group/erp/data/olaf.hauk/MEG/WakemanHensonEMEG/Resolution/ResolutionMetrics/Figures'  # output directory for images
inpath = '/group/erp/data/olaf.hauk/MEG/WakemanHensonEMEG/Resolution/ResolutionMetrics/fsaverage'

# List of STC files, with thresholds for visualisation in mne_make_movie (min, mid, max)
filelist = [
    # ### Peak localisation Error, raw distributions
    # ['PSF_MNE_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_MNE_011_locerr_peak_EEGMEG_loo0_dep80', 0, 2.5, 5.],
    # ['PSF_sLORETA_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_dSPM_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_LCMV_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_MNE_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_MNE_011_locerr_peak_EEGMEG_loo0_dep80', 0, 2.5, 5.],
    # ['CTF_sLORETA_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_dSPM_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_LCMV_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],

    ['PSF_MNE_011_locerr_peak_MEG_loo0_dep0', 0, 2.5, 5.],
    ['PSF_MNE_011_locerr_peak_EEGMEG-MEG_loo0_dep0', 0, 0.5, 1.],

    # ### Peak localisation Error, constrasts
    # ['PSF_MNE-dSPM_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_MNE-sLORETA_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_MNE-LCMV_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_dep0-80_011_locerr_peak_EEGMEG_loo0', 0, 2.5, 5.],
    # ['CTF_MNE-dSPM_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_MNE-sLORETA_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_MNE-LCMV_011_locerr_peak_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_dep0-80_011_locerr_peak_EEGMEG_loo0', 0, 2.5, 5.],

    
    # ### Spatial Deviation, raw distributions
    # ['PSF_MNE_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['PSF_MNE_011_width_sd_EEGMEG_loo0_dep80', 0, 1.5, 3.],
    # ['PSF_sLORETA_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['PSF_dSPM_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['PSF_LCMV_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_MNE_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_MNE_011_width_sd_EEGMEG_loo0_dep80', 0, 1.5, 3.],
    # ['CTF_sLORETA_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_dSPM_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_LCMV_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],

    ['PSF_MNE_011_width_sd_MEG_loo0_dep0', 0, 1.5, 3.],
    ['PSF_MNE_011_width_sd_EEGMEG-MEG_loo0_dep0', 0, 1.5, 3.],

    # ### Spatial Deviation, constrasts
    # ['PSF_MNE-dSPM_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['PSF_MNE-sLORETA_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['PSF_MNE-LCMV_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['PSF_dep0-80_011_width_sd_EEGMEG_loo0', 0, 1.5, 3.],
    # ['CTF_MNE-dSPM_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_MNE-sLORETA_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_MNE-LCMV_011_width_sd_EEGMEG_loo0_dep0', 0, 1.5, 3.],
    # ['CTF_dep0-80_011_width_sd_EEGMEG_loo0', 0, 1.5, 3.],


    # ### Maximum Radius, raw distributions
    # ['PSF_MNE_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['PSF_MNE_011_width_maxrad_EEGMEG_loo0_dep80', 0, 4., 8.],
    # ['PSF_sLORETA_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['PSF_dSPM_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['PSF_LCMV_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['CTF_MNE_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['CTF_MNE_011_width_maxrad_EEGMEG_loo0_dep80', 0, 4., 8.],
    # ['CTF_sLORETA_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['CTF_dSPM_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],
    # ['CTF_LCMV_011_width_maxrad_EEGMEG_loo0_dep0', 0, 4., 8.],

    # ### Maximum Radius, constrasts
    # ['PSF_MNE-dSPM_011_width_maxrad_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_MNE-sLORETA_011_width_maxrad_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_MNE-LCMV_011_width_maxrad_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['PSF_dep0-80_011_width_maxrad_EEGMEG_loo0', 0, 2.5, 5.],
    # ['CTF_MNE-dSPM_011_width_maxrad_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_MNE-sLORETA_011_width_maxrad_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_MNE-LCMV_011_width_maxrad_EEGMEG_loo0_dep0', 0, 2.5, 5.],
    # ['CTF_dep0-80_011_width_maxrad_EEGMEG_loo0', 0, 2.5, 5.],

    
    # ### Overall Amplitude, raw distributions
    # ['PSF_MNE_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['PSF_MNE_011_amplitude_sum_EEGMEG_loo0_dep80', 0, .5, 1.],
    # ['PSF_sLORETA_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['PSF_dSPM_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['PSF_LCMV_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['CTF_MNE_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['CTF_MNE_011_amplitude_sum_EEGMEG_loo0_dep80', 0, .5, 1.],
    # ['CTF_sLORETA_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['CTF_dSPM_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],
    # ['CTF_LCMV_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .5, 1.],

    ['PSF_MNE_011_amplitude_sum_MEG_loo0_dep0', 0, .5, 1.],
    ['PSF_MNE_011_amplitude_sum_EEGMEG-MEG_loo0_dep0', 0, .25, 0.5]
    

    # ### Overall Amplitude, constrasts
    # ['PSF_MNE-dSPM_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .25, .5],
    # ['PSF_MNE-sLORETA_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .25, .5],
    # ['PSF_MNE-LCMV_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .25, .5],
    # ['PSF_dep0-80_011_amplitude_sum_EEGMEG_loo0', 0, .25, .5],
    # ['CTF_MNE-dSPM_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .25, .5],
    # ['CTF_MNE-sLORETA_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .25, .5],
    # ['CTF_MNE-LCMV_011_amplitude_sum_EEGMEG_loo0_dep0', 0, .25, .5],
    # ['CTF_dep0-80_011_amplitude_sum_EEGMEG_loo0', 0, .25, .5]
]


# parameters for mne_make_movie common to all files
subject = '--subject fsaverage'
smooth = '--smooth 5'
pick = '--pick 0'
spm = '--spm'
comments = '--nocomments'

print('Visualising %d files.' % len(filelist))

# Keep track of error from mne_make_movie output
errors = []

for ff in filelist:

    filein = op.join(inpath, ff[0])
    fileout = op.join(outpath, ff[0])

    # parameters for mne_make_movie for individual files
    stcin = '--stcin %s' % filein
    jpg = '--jpg %s' % fileout
    fthresh = '--fthresh %f  --fmid %f  --fmax %f' % (ff[1], ff[2], ff[3])

    # create string for mne_make_movie command
    cmd = 'mne_make_movie %s %s %s %s %s %s %s %s' % (subject, smooth, pick, spm, comments, stcin, jpg, fthresh)

    # execute mne_make_movie command in linux
    osout = os.system(cmd)

    # if error occurs, keep track
    if not osout==0:

        errors.append([ff, osout])

if not errors == []:

    print('\n\nThe following errors occurred while executing mne_make_movie:\n')
    
    for ee in errors:

        print('%s, %s' % (ee[0][0], ee[1]))

else:

    print('No errors occurred while executing mne_make_movie!')

# Done all