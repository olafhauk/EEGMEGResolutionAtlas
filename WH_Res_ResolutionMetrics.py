"""
=========================================================
Compute whole-brain resolution metrics for WH data set.
E.g. run WH_Res_ResolutionMetrics.py WH_MNE_Resolution_config 0 locerr peak
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

# for filenames
st_duration = C.res_st_duration
origin = C.res_origin

# resolution metric to be computed
# read variables specified via qsub
if len(sys.argv)>3: # if additional variable specified

    metric_type = sys.argv[3] # metric type, e.g. 'width' or 'locerr'

    metric = sys.argv[4] # exact metric, e.g. 'RMS'


for sbj in sbj_ids:

    subject = 'Sub%02d' % sbj

    print('###\nWorking hard on %s.\n###' % (subject))

    fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, '')

    # create output path if necessary
    if not os.path.exists(fname_stc):

        os.mkdir(fname_stc)


    fwd_fname = C.fname_ForwardSolution(C, subject, 'EEGMEG')

    print('###\nReading EEGMEG forward solutions: %s .\n###' % (fwd_fname))

    fwd = mne.read_forward_solution(fwd_fname)

    fwd = mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True)

    # used vertex numbers, later assumed to be the same as for invop
    vertno_0 = fwd['src'][0]['vertno']
    vertno_1 = fwd['src'][1]['vertno']
    vertno = [vertno_0, vertno_1]

    # locations of used vertices in forward (and inverse) solution
    locations_0 = fwd['src'][0]['rr'][vertno_0,:]
    locations_1 = fwd['src'][1]['rr'][vertno_1,:]
    locations = np.vstack([locations_0, locations_1])


    # covariance matrix (filter with wildcard)
    fname_cov = C.fname_cov(C, subject, st_duration, origin, C.res_cov_latwin, C.inv_method, '*')

    # method may be underspecified, since it may be ranked differently for different subjects
    fname_cov = glob.glob(fname_cov)[0] # be careful if multiple options present

    print('###\nReading noise covariance matrix from: %s.\n###' % (fname_cov))

    noise_cov = mne.read_cov(fname_cov)


    # initialise dict to collect results
    stcs = {'PSF': {}, 'CTF': {}}
    for functype in stcs.keys():
        for inv_type in C.res_inv_types:
            stcs[functype][inv_type] = {}

    # iterate over different combinations of sensors
    for (eeg,meg,modal) in [(True,True,'EEGMEG'), (False,True,'MEG'), (True,False,'EEG')]:

        fwd_use = mne.pick_types_forward(fwd, meg=meg, eeg=eeg)

        mytext = modal + '_fxd_nodep'
        inv_fname = C.fname_InverseOperator(C, subject, st_duration, origin, C.res_cov_latwin, mytext)

        print('\n###\nReading inverse operator from %s.\n###\n' % (inv_fname))

        invop = mne.minimum_norm.read_inverse_operator(inv_fname)

        info = invop['info']
        ch_names = info['ch_names']
        info['sfreq'] = 1000.
        info['projs'] = noise_cov['projs']
        info['comps'] = ''
                
        # restrict noise covariance to channels in forward solution
        noise_cov_use = mne.cov.pick_channels_cov(noise_cov, ch_names)

        # regularise noise covariance
        if C.inv_method == 'empirical': # if unregularised

            noise_cov_use = mne.cov.regularize(noise_cov_use, info, mag=C.res_lambda_empirical['mag'],
                                                    grad=C.res_lambda_empirical['grad'], eeg=C.res_lambda_empirical['eeg'])


        for inv_type in C.res_inv_types: # 'MNE', 'LCMV' etc.

            print('\n###\nComputing resolution matrix for %s.\n###\n' % (inv_type))

            if inv_type in ['MNE', 'sLORETA', 'dSPM']: # minimum-norm-based methods

                resmat = R.make_resolution_matrix(fwd_use, invop, method=inv_type, lambda2=1./3.**2)

            elif inv_type == 'LCMV': # LCMV beamformer

                # NOTE: here noise and data cov identical
                resmat = R.make_resolution_matrix_lcmv(fwd_use, info, noise_cov, noise_cov)

            
            # turn resolution matrix to STC and save
            for (mat, matname) in [(resmat, 'PSF'), (resmat.T, 'CTF')]:

                stc = R.resmat_to_stc(mat, fwd['src'])

                mytext = matname + 'plot_' + inv_type + '_' + modal

                fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                print('###\nWriting %s as STC file to: %s.\n###' % (matname, fname_stc))

                stc.save(fname_stc)


            # compute CTFs and PSFs
            for (axis, functype) in ((0, 'CTF'), (1,'PSF')):
       
                # NOTE: edit metric, metric_text and metric_mat at the same time
                if metric_type == 'locerr':

                    print('###\nComputing localisation error metric for %s.\n###' % functype)

                    # compute requested resolution metric
                    metric_mat = R.localisation_error(resmat, locations, axis=axis, metric=metric)

                    # convert to cm for visualisation (NOTE: may depend on metric)
                    metric_mat = 100.*metric_mat

                elif metric_type == 'width':

                    print('###\nComputing spatial width metric for %s.\n###' % functype)

                    metric_mat = R.spatial_width(resmat, locations, axis=axis, metric=metric)

                    # convert to cm (NOTE it's sqrt of distance in Molins SD metric, hence 10)
                    metric_mat = 100.*metric_mat

                elif metric_type == 'amplitude':

                    print('###\nComputing amplitude metric for %s.\n###' % functype)

                    metric_mat = R.relative_amplitude(resmat, locations, axis, metric=metric)

                # for filename
                metric_text = '_' + metric_type + '_' + metric + '_'

                # fake multiple time steps
                # cm for visualisation
                metric_mat_rep = np.repeat(metric_mat[:,np.newaxis], 5, axis=1)

                # convert norms to source estimate
                stc_metric = mne.SourceEstimate(metric_mat_rep, vertno, tmin=0., tstep=1.)

                stcs[functype][inv_type][modal] = stc_metric

                mytext = functype + '_' + inv_type + metric_text + modal

                fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                print('###\nWriting STC file to: %s.\n###' % (fname_stc))

                stc_metric.save(fname_stc)

    print('\n###\nContrasting modalities.\n###')

    for functype in stcs.keys(): # for function types specified above (PSF, CTF)

        for inv_type in C.res_inv_types: # 'norm', 'lcmv'

            for (modal1,modal2) in [('EEGMEG', 'MEG'), ('EEGMEG', 'EEG'), ('MEG','EEG')]:
           
                stc_contr = stcs[functype][inv_type][modal1] - stcs[functype][inv_type][modal2]

                modal = modal1 + '-' + modal2

                stcs[functype][inv_type][modal] = stc_contr

                mytext = functype + '_' + inv_type + metric_text + modal

                fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                print('Saving STC to: %s.' % fname_stc)

                stc_contr.save(fname_stc)


    ### Visualisation:

    # clim = {'kind': 'value', 'pos_lims': (0, 0.5, 1)}

    # stc_norm['EEGMEG-MEG'].plot(subject=subject, subjects_dir=C.subjects_dir, hemi='both',
    #                         time_viewer=True, transparent=False, colormap='mne', clim=clim)


# Done