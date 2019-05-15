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
    fname_cov = C.fname_cov(C, subject, st_duration, origin, C.res_cov_latwin, C.inv_cov_method, '*')

    # method may be underspecified, since it may be ranked differently for different subjects
    fname_cov = glob.glob(fname_cov)[0] # be careful if multiple options present

    print('###\nReading noise covariance matrix from: %s.\n###' % (fname_cov))

    noise_cov = mne.read_cov(fname_cov)

    # INITIALISE stcs dict to collect results
    # the more the merrier...
    stcs = {'PSF': {}, 'CTF': {}}

    for functype in stcs.keys(): # PSF/CTF

        # contrasts for different inverse methods and subtractions (also depth weighted)
        for inv_type in C.res_inv_types + [x+'-'+y for (x,y) in C.meth_contr] + ['dep'+str(int(100*x))+'-'+str(int(100*y)) for (x,y) in C.meth_contr_dep]: # inverse methods

            stcs[functype][inv_type] = {}

            for lambda2 in C.res_lambda2s: # regularisation parameters

                lamb2_str = str(lambda2).replace('.', '')
                if len(lamb2_str) > 3:
                    lamb2_str = lamb2_str[:3]

                stcs[functype][inv_type][lamb2_str] = {}

                for modal in C.modalities + [x+'-'+y for (x,y) in C.modal_contr]: # ['EEGMEG', 'MEG', 'EEG', 'EEGMEG-MEG', 'EEGMEG-EEG', 'MEG-EEG']:

                    stcs[functype][inv_type][lamb2_str][modal] = {}

                    for loose in C.inv_loose: # orientation constraint

                        if loose == None: loose = 0
                        loo_str = 'loo%s' % str(int(100*loose))

                        stcs[functype][inv_type][lamb2_str][modal][loo_str] = {}

                        for depth in C.inv_depth: # depth weighting

                            if depth == None: depth = 0
                            dep_str = 'dep%s' % str(int(100*depth))

                            stcs[functype][inv_type][lamb2_str][modal][loo_str][dep_str] = {}


    # iterate over different combinations of sensors
    modal_iter = []
    for modal in C.modalities:
        if modal == 'EEGMEG':
            modal_iter.append([True, True, 'EEGMEG'])
        elif modal == 'MEG':
            modal_iter.append([False, True, 'MEG'])
        elif modal == 'EEG':
            modal_iter.append([True, False, 'EEG'])
        else:
            print('Modality not defined: %s.' % modal)

    for (eeg,meg,modal) in modal_iter:

        fwd_use = mne.pick_types_forward(fwd, meg=meg, eeg=eeg)

        info = fwd_use['info']
        ch_names = info['ch_names']
        info['sfreq'] = 1000.
                
        # restrict noise covariance to channels in forward solution
        noise_cov_use = mne.cov.pick_channels_cov(noise_cov, ch_names)

        info['projs'] = noise_cov_use['projs']
        info['comps'] = ''

        # regularise noise covariance
        if C.inv_cov_method == 'empirical': # if unregularised

            noise_cov_use = mne.cov.regularize(noise_cov_use, info, mag=C.res_lambda_empirical['mag'],
                                                    grad=C.res_lambda_empirical['grad'], eeg=C.res_lambda_empirical['eeg'])

        # Compute some null space components (not sure about it yet)
        nullspace, leadfield, pseudoinv, null_proj = R.leadfield_nullspace(fwd_use, noise_cov_use, regpar=0.)

        stc = R.mat_to_stc(nullspace, fwd['src'])

        mytext = 'Nullspace_' + modal

        fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

        print('###\nWriting nullspace for %s as STC file to: %s.\n###' % (modal, fname_stc))

        stc.save(fname_stc)


        # iterate over inverse operator types
        for loose in C.inv_loose: # orientation constraint

            for depth in C.inv_depth: # depth weighting

                inv_fname = C.fname_InverseOperator(C, subject, st_duration, origin, C.res_cov_latwin,
                                                    modal, loose, depth)

                print('\n###\nReading inverse operator from %s.\n###\n' % (inv_fname))

                invop = mne.minimum_norm.read_inverse_operator(inv_fname)

                for inv_type in C.res_inv_types: # 'MNE', 'LCMV' etc.

                    for lambda2 in C.res_lambda2s: # regularisation parameters

                        lamb2_str = str(lambda2).replace('.', '')
                        if len(lamb2_str) > 3:
                            lamb2_str = lamb2_str[:3]

                        print('\n###\nComputing resolution matrix for %s.\n###\n' % (inv_type))

                        if inv_type in ['MNE', 'sLORETA', 'dSPM']: # minimum-norm-based methods

                            resmat = R.make_resolution_matrix(fwd_use, invop, method=inv_type, lambda2=lambda2)

                        elif inv_type == 'LCMV': # LCMV beamformer

                            # NOTE: here noise and data cov identical
                            resmat = R.make_resolution_matrix_lcmv(fwd_use, info, noise_cov_use, noise_cov_use)

                        
                        # turn resolution matrix to STC and save
                        for (mat, matname) in [(resmat, 'PSF'), (resmat.T, 'CTF')]:

                            stc = R.mat_to_stc(mat, fwd['src'])

                            mytext = '%s_plot_%s_%s_%s' % (matname, inv_type, lamb2_str, modal)

                            fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                            print('###\nWriting %s as STC file to: %s.\n###' % (matname, fname_stc))

                            stc.save(fname_stc)


                        # compute CTFs and PSFs
                        for (axis, functype) in [(0, 'PSF'), (1,'CTF')]:
                   
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

                                # convert to cm
                                metric_mat = 100.*metric_mat

                            elif metric_type == 'amplitude':

                                print('###\nComputing amplitude metric for %s.\n###' % functype)

                                metric_mat = R.relative_amplitude(resmat, locations, axis, metric=metric)

                            # fake multiple time steps
                            # cm for visualisation
                            metric_mat_rep = np.repeat(metric_mat[:,np.newaxis], 5, axis=1)

                            # convert norms to source estimate
                            stc_metric = mne.SourceEstimate(metric_mat_rep, vertno, tmin=0., tstep=1.)

                            # for filename
                            metric_text = metric_type + '_' + metric

                            if loose == None: loose = 0

                            loo_str = 'loo%s' % str(int(100*loose))

                            if depth == None: depth = 0

                            dep_str = 'dep%s' % str(int(100*depth))

                            mytext = '%s_%s_%s_%s_%s_%s_%s' % (functype, inv_type, lamb2_str, metric_text, modal, loo_str, dep_str)

                            fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                            print('###\nWriting STC file to: %s.\n###' % (fname_stc))

                            stc_metric.save(fname_stc)

                            # collect results in dict
                            stcs[functype][inv_type][lamb2_str][modal][loo_str][dep_str] = stc_metric

    
    print('\n###\nContrasting modalities.\n###')

    for functype in stcs.keys(): # for function types specified above (PSF, CTF)

        for inv_type in C.res_inv_types: # 'MNE', 'dSPM', 'lcmv' etc.

            for lambda2 in C.res_lambda2s: # regularisation parameters

                lamb2_str = str(lambda2).replace('.', '')
                if len(lamb2_str) > 3:
                    lamb2_str = lamb2_str[:3]

                for (modal1,modal2) in C.modal_contr: # [('EEGMEG', 'MEG'), ('EEGMEG', 'EEG'), ('MEG','EEG')]:

                    # iterate over inverse operator types
                    for loose in C.inv_loose: # orientation constraint

                        for depth in C.inv_depth: # depth weighting

                            if loose == None: loose = 0

                            loo_str = 'loo%s' % str(int(100*loose))

                            if depth == None: depth = 0

                            dep_str = 'dep%s' % str(int(100*depth))
               
                            stc_contr = stcs[functype][inv_type][lamb2_str][modal1][loo_str][dep_str] - \
                                            stcs[functype][inv_type][lamb2_str][modal2][loo_str][dep_str]

                            modal = modal1 + '-' + modal2

                            stcs[functype][inv_type][lamb2_str][modal][loo_str][dep_str] = stc_contr

                            mytext = '%s_%s_%s_%s_%s_%s_%s' % (functype, inv_type, lamb2_str, metric_text, modal, loo_str, dep_str)

                            fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                            print('Saving STC to: %s.' % fname_stc)

                            stc_contr.save(fname_stc)

                    print "Nothing."


    print('\n###\nContrasting inverse methods.\n###')

    for functype in stcs.keys(): # for function types specified above (PSF, CTF)

        # methods to subtract from each other, meth1-meth2
        for modal in ['EEGMEG']:

            for (meth1,meth2) in C.meth_contr: # [('MNE','dSPM'), ('MNE','sLORETA'), ('dSPM', 'sLORETA'), ('MNE', 'LCMV'), ('dSPM', 'LCMV'), ('sLORETA','LCMV')]:

                for lambda2 in C.res_lambda2s: # regularisation parameters

                    lamb2_str = str(lambda2).replace('.', '')
                    if len(lamb2_str) > 3:
                        lamb2_str = lamb2_str[:3]

                    # iterate over inverse operator types
                    for loose in C.inv_loose: # orientation constraint

                        for depth in C.inv_depth: # depth weighting

                            if loose == None: loose = 0

                            loo_str = 'loo%s' % str(int(100*loose))

                            if depth == None: depth = 0

                            dep_str = 'dep%s' % str(int(100*depth))
               
                            stc_contr = stcs[functype][meth1][lamb2_str][modal][loo_str][dep_str] - \
                                            stcs[functype][meth2][lamb2_str][modal][loo_str][dep_str]

                            meth = meth1 + '-' + meth2

                            stcs[functype][meth][lamb2_str][modal][loo_str][dep_str] = stc_contr

                            mytext = '%s_%s_%s_%s_%s_%s_%s' % (functype, meth, lamb2_str, metric_text, modal, loo_str, dep_str)

                            fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                            print('Saving STC to: %s.' % fname_stc)

                            stc_contr.save(fname_stc)

            
            # Depth weighting is separate
            # Do subtractions if specified for different depth weightings
            if not C.meth_contr_dep == []:

                loose = 0
                loo_str = 'loo%s' % str(int(100*loose))

                # For filenames
                meth1 = 'MNE'
                meth2 = 'MNE'

                for (dep1,dep2) in C.meth_contr_dep: # what to subtract from each other
                  
                    dep_str1 = 'dep%s' % str(int(100*dep1))
                    dep_str2 = 'dep%s' % str(int(100*dep2))

                    for lambda2 in C.res_lambda2s: # regularisation parameters

                        lamb2_str = str(lambda2).replace('.', '')
                        if len(lamb2_str) > 3:
                            lamb2_str = lamb2_str[:3]
                                               
                        # subtract depth-weighted from non-weighted MNE
                        stc_contr = stcs[functype][meth1][lamb2_str][modal][loo_str][dep_str1] - \
                                                    stcs[functype][meth2][lamb2_str][modal][loo_str][dep_str2]

                        # for filename
                        meth = 'dep%s-%s' % (dep_str1[3:], dep_str2[3:])

                        stcs[functype][meth][lamb2_str][modal][loo_str][dep_str2] = stc_contr

                        mytext = '%s_%s_%s_%s_%s_%s' % (functype, meth, lamb2_str, metric_text, modal, loo_str)

                        fname_stc = C.fname_STC(C, 'ResolutionMetrics', subject, mytext)

                        print('Saving STC to: %s.' % fname_stc)

                        stc_contr.save(fname_stc)


    ### Visualisation:

    # clim = {'kind': 'value', 'pos_lims': (0, 0.5, 1)}

    # stc_norm['EEGMEG-MEG'].plot(subject=subject, subjects_dir=C.subjects_dir, hemi='both',
    #                         time_viewer=True, transparent=False, colormap='mne', clim=clim)


# Done