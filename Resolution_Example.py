"""
=========================================================
Compute resolution metrics for MNE sample data set.
Peak Localisation Error (PLE) and Spatial Deviation (SD)
are computed for L2-MNE and sLORETA, respectively.
sLORETA has zero PLE, therefore outperforming L2-MNE.
However, L2-MNE outperforms sLORETA with respect to SD.
The CTFs for L2-MNE and sLORETA have the same shape,
thus localisation error and width metrics are the same.
For L2-MNE, PSFs and CTFs are the same, hence the metrics
are, too.
=========================================================
OH, May 2019
"""
import numpy as np

import mne
from mne.datasets import sample

import importlib

print(__doc__)

# get functions for resolution matrix, metrics etc.
Res = importlib.import_module('Resolution_Example_Functions')
reload(Res)

# get data from MNE sample dataset
data_path = sample.data_path()
subjects_dir = data_path + '/subjects/'

fname_fwd = data_path + '/MEG/sample/sample_audvis-meg-oct-6-fwd.fif'

fname_inv = data_path + '/MEG/sample/sample_audvis-meg-oct-6-meg-inv.fif'

fname_evoked = data_path + '/MEG/sample/sample_audvis-ave.fif'

fname_cov = data_path + '/MEG/sample/sample_audvis-shrunk-cov.fif'

# read forward solution
forward = mne.read_forward_solution(fname_fwd)

# only use normal components in forward solution
forward = mne.convert_forward_solution(forward, surf_ori=True,
                                        force_fixed=True)

# noise covariance matrix
noise_cov = mne.read_cov(fname_cov)

# evoked data for info
evoked = mne.evoked.read_evokeds(fname_evoked, 0)

# make inverse operator from forward solution
inverse_operator = mne.minimum_norm.make_inverse_operator(info=evoked.info,
        forward=forward, noise_cov=noise_cov, loose=0., depth=None, fixed=True)

# regularisation parameter based on SNR
snr = 3.0
lambda2 = 1.0 / snr ** 2

# vertices used in forward and inverse operator
vertno_lh = forward['src'][0]['vertno']
vertno_rh = forward['src'][1]['vertno']
vertno = [vertno_lh, vertno_rh]

# locations corresponding to vertices for both hemispheres
locations_lh = forward['src'][0]['rr'][vertno_lh,:]
locations_rh = forward['src'][1]['rr'][vertno_rh,:]
locations = np.vstack([locations_lh, locations_rh])

for (func,axis) in zip(['PSF','CTF'], [0,1]):

    # initialise lists for resolution metrics
    locerr, width = [], []

    for method in ['MNE', 'sLORETA']:

        # compute resolution matrix
        RM = Res.make_resolution_matrix(forward, inverse_operator, method=method,
                                        lambda2=lambda2)

        # compute peak localisation error (cm)
        # axis=0: PSF (columns of RM); axis=1: CTF (rows)
        locerr.append(100.*Res.localisation_error(RM, locations, axis=axis,
                                                metric='peak')[:,np.newaxis])

        # convert to source estimate
        stc = mne.SourceEstimate(locerr[-1], vertno, tmin=0., tstep=1.)

        stc.save('%s_locerr_%s' % (func, method))

        ### compute spatial deviation

        # compute spatial deviation (cm)
        width.append(100.*Res.spatial_width(RM, locations, axis=axis,
                                            metric='sd')[:,np.newaxis])

        # convert to source estimate
        stc = mne.SourceEstimate(width[-1], vertno, tmin=0., tstep=1.)

        stc.save('%s_width_%s' % (func, method))

    # Compute differences between methods
    stc = mne.SourceEstimate(locerr[0]-locerr[1], vertno, tmin=0., tstep=1.)

    stc.save('%s_locerr_MNE-sLORETA' % func)

    stc = mne.SourceEstimate(width[0]-width[1], vertno, tmin=0., tstep=1.)

    stc.save('%s_width_MNE-sLORETA' % func)