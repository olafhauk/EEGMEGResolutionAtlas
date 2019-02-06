"""
=========================================================
Morph STCs of resolution metrics for W&H data set.
e.g.: run WH_Res_MorphSTC.py WH_MNE_Resolution_config 11 SensitivityMaps SensMap RMS
or: run WH_Res_MorphSTC.py WH_MNE_Resolution_config 11 ResolutionMetrics locerr_peak
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

# get subject ID to process
# qsub start at 0, thus +1 here
sbj_ids = [int(sys.argv[2]) + 1]

# hack to have variables via qsub
stc_path, stc_type, metric = '', '', ''

# read variables specified via qsub
if len(sys.argv)>3: # if additional variable specified

    stc_path = sys.argv[3] # pathname (e.g. 'SensitivityMap')

if len(sys.argv) > 4:

    stc_type = sys.argv[4] # beginning of filename (e.g. 'LocErrPeak')

if len(sys.argv) > 5:

    metric = sys.argv[5] # end of filename (e.g. 'RMS')

    print('###\nChosen variables: %s, %s, %s.\n###' % (stc_path, stc_type, metric))

# for filenames
st_duration = C.res_st_duration
origin = C.res_origin

# only one morph_mat per subject needed
morph_mat = []

###
for sbj in sbj_ids:

    subject = 'Sub%02d' % sbj

    print('###\nAbout to morph STCs for %s.\n###' % (subject))

    for modality in ['EEGMEG', 'MEG', 'EEG', 'EEGMEG-MEG', 'EEGMEG-EEG']: # EEG/MEG/EEGMEG

        # contrasts for inverse methods only computed for EEGMEG
        if modality == 'EEGMEG':

            res_inv_types = C.res_inv_types + ['MNE-dSPM', 'MNE-sLORETA', 'dSPM-sLORETA', 'MNE-LCMV', 'dSPM-LCMV', 'sLORETA-LCMV']

        else:

            res_inv_types = C.res_inv_types

        for inv_type in res_inv_types: # 'MNE', 'LCMV' etc.

            # for CTFs and PSFs
            for functype in ['CTF', 'PSF']:

                # iterate over inverse operator types
                for loose in C.inv_loose: # orientation constraint

                    for depth in C.inv_depth: # depth weighting

                        print('\n###\nDoing %ss for %s and %s.\n###\n' % (functype, inv_type, modality))

                        if loose == None: loose = 0

                        loo_str = '_loo%s' % str(int(100*loose))

                        if depth == None: depth = 0

                        dep_str = '_dep%s' % str(int(100*depth))

                        mytext = functype + '_' + inv_type + '_' + stc_type + '_' + modality + loo_str + dep_str

                        if metric != '':

                            mytext = mytext + '_' + metric

                        fname_stc = C.fname_STC(C, stc_path, subject, mytext)

                        fname_morph = C.fname_STC(C, stc_path, subject, mytext + '_mph')

                        # read existing source estimate
                        print('Reading: %s.' % fname_stc)
                        stc = mne.read_source_estimate(fname_stc, subject)

                        # compute morph_mat only once per subject
                        if morph_mat == []:

                            vertices_to = mne.grade_to_vertices(subject=C.stc_morph, grade=5, subjects_dir=C.subjects_dir)

                            morph_mat = mne.compute_morph_matrix(subject_from=subject, subject_to=C.stc_morph,
                                                                vertices_from=stc.vertices, vertices_to=vertices_to,
                                                                subjects_dir=C.subjects_dir)

                        # Morphing to standard brain
                        morphed = mne.morph_data_precomputed(subject_from=subject, subject_to=C.stc_morph, stc_from=stc,
                                                              vertices_to=vertices_to, morph_mat=morph_mat)

                        print('Writing morphed to: %s.' % fname_morph)
                        morphed.save(fname_morph)

    # Depth weighting is separate

    for modality in ['EEGMEG']: # EEG/MEG/EEGMEG

        res_inv_types = ['MNE-MNE40', 'MNE-MNE80']

        for inv_type in res_inv_types: # 'MNE', 'LCMV' etc.

            # for CTFs and PSFs
            for functype in ['CTF', 'PSF']:

                loose = 0
                loo_str = '_loo%s' % str(int(100*loose))

                mytext = functype + '_' + inv_type + '_' + stc_type + '_' + modality + loo_str

                if metric != '':

                    mytext = mytext + '_' + metric

                fname_stc = C.fname_STC(C, stc_path, subject, mytext)

                fname_morph = C.fname_STC(C, stc_path, subject, mytext + '_mph')

                # read existing source estimate
                print('Reading: %s.' % fname_stc)
                stc = mne.read_source_estimate(fname_stc, subject)

                # compute morph_mat only once per subject
                if morph_mat == []:

                    vertices_to = mne.grade_to_vertices(subject=C.stc_morph, grade=5, subjects_dir=C.subjects_dir)

                    morph_mat = mne.compute_morph_matrix(subject_from=subject, subject_to=C.stc_morph,
                                                        vertices_from=stc.vertices, vertices_to=vertices_to,
                                                        subjects_dir=C.subjects_dir)

                # Morphing to standard brain
                morphed = mne.morph_data_precomputed(subject_from=subject, subject_to=C.stc_morph, stc_from=stc,
                                                      vertices_to=vertices_to, morph_mat=morph_mat)

                print('Writing morphed to: %s.' % fname_morph)
                morphed.save(fname_morph)

# Done