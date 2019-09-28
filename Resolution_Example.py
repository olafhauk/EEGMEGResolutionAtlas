"""
=============================================================
Compute resolution metrics for MNE sample data set.
Peak Localisation Error (PLE) is computed for L2-MNE
and sLORETA, and for point-spread and cross-talk functions
(PSFs and CTFs), respectively.
sLORETA has zero PLE for PSFs, therefore outperforming L2-MNE.
However, the CTFs for L2-MNE and sLORETA have the same shape,
thus PLE for CTFs is the same (and non-zero) for both
MNE and sLORETA.
==============================================================
OH, Sep 2019
"""
import numpy as np

import mne
from mne.datasets import sample

import importlib

print(__doc__)

# get functions for resolution matrix, metrics etc.
Res = importlib.import_module('Resolution_Example_Functions')
importlib.reload(Res)

R = importlib.import_module('WH_Resolution_Functions')
importlib.reload(R)

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

# compute resolution matrix for MNE
RM_MNE = Res.make_resolution_matrix(forward, inverse_operator, method='MNE',
                                lambda2=lambda2)

# compute resolution matrix for sLORETA
RM_LOR = Res.make_resolution_matrix(forward, inverse_operator, method='sLORETA',
                                lambda2=lambda2)

# compute peak localisation error (cm) for MNE's PSF
# axis=0: PSF (columns of RM); axis=1: CTF (rows)
locerr_PSF_MNE = 100.*Res.localisation_error(RM_MNE, locations, axis=0,
                                                metric='peak')[:,np.newaxis]

# convert to source estimate
stc_PSF_MNE = mne.SourceEstimate(locerr_PSF_MNE, vertno, tmin=0., tstep=1.)


# compute peak localisation error (cm) for sLORETA's PSF
# axis=0: PSF (columns of RM); axis=1: CTF (rows)
locerr_PSF_LOR = 100.*Res.localisation_error(RM_LOR, locations, axis=0,
                                                metric='peak')[:,np.newaxis]
# convert to source estimate
stc_PSF_LOR = mne.SourceEstimate(locerr_PSF_LOR, vertno, tmin=0., tstep=1.)


# compute peak localisation error (cm) for MNE's CTF
# axis=1: PSF (columns of RM); axis=1: CTF (rows)
locerr_CTF_MNE = 100.*Res.localisation_error(RM_MNE, locations, axis=1,
                                                metric='peak')[:,np.newaxis]

# convert to source estimate
stc_CTF_MNE = mne.SourceEstimate(locerr_CTF_MNE, vertno, tmin=0., tstep=1.)


# compute peak localisation error (cm) for sLORETA's CTF
# axis=1: PSF (columns of RM); axis=1: CTF (rows)
locerr_CTF_LOR = 100.*Res.localisation_error(RM_LOR, locations, axis=1,
                                                metric='peak')[:,np.newaxis]

# convert to source estimate
stc_CTF_LOR = mne.SourceEstimate(locerr_CTF_LOR, vertno, tmin=0., tstep=1.)

# Visualise
from mayavi import mlab


brain_PSF_MNE = stc_PSF_MNE.plot('sample', 'inflated', 'both',
                                    subjects_dir=subjects_dir, 
                                    clim=dict(kind='value', lims=(0, 2, 4)),
                                    title='PLE PSF MNE', figure=1)
mlab.title('PLE PSF MNE')

brain_PSF_LOR = stc_PSF_LOR.plot('sample', 'inflated', 'both',
                                    subjects_dir=subjects_dir, 
                                    clim=dict(kind='value', lims=(0, 2, 4)),
                                    title='PLE PSF sLORETA', figure=2)
mlab.title('PLE PSF sLORETA')

brain_CTF_MNE = stc_CTF_MNE.plot('sample', 'inflated', 'both',
                                    subjects_dir=subjects_dir, 
                                    clim=dict(kind='value', lims=(0, 2, 4)),
                                    title='PLE CTF MNE', figure=3)
mlab.title('PLE CTF MNE')

brain_CTF_LOR = stc_CTF_LOR.plot('sample', 'inflated', 'both',
                                    subjects_dir=subjects_dir, 
                                    clim=dict(kind='value', lims=(0, 2, 4)),
                                    title='PLE CTF sLORETA', figure=4)
mlab.title('PLE CTF sLORETA')