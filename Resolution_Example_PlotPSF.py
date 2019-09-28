"""
=============================================================
Compute Point-Spread and Cross-Talk Functions (PSF and CTF)
for sLORETA. 
The PSF has zero peak localisation error (i.e. the peak occurs
at the correct vertex), but the CTF does not.
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

# R = importlib.import_module('WH_Resolution_Functions')
# importlib.reload(R)

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

# compute resolution matrix for sLORETA
RM_LOR = Res.make_resolution_matrix(forward, inverse_operator, method='sLORETA',
                                lambda2=lambda2)

# which PSF/CTF to plot
idx = 1000

# column of resolution matrix
PSF = RM_LOR[:,idx]
# row of resolution matrix
CTF = RM_LOR[idx,:]


# convert to source estimate
stc_PSF = mne.SourceEstimate(PSF, vertno, tmin=0., tstep=1.)
stc_CTF = mne.SourceEstimate(CTF, vertno, tmin=0., tstep=1.)

# maximum value for scaling
mp = np.abs(PSF).max()
mc = np.abs(CTF).max()

# vertex corresponding to PSF/CTF
vertex = vertno_lh[idx]

vert_max_psf = vertno_lh[PSF.argmax()]
vert_max_ctf = vertno_lh[CTF.argmax()]

# Visualise
from mayavi import mlab

brain_PSF = stc_PSF.plot('sample', 'inflated', 'lh',
                                    subjects_dir=subjects_dir, 
                                    clim=dict(kind='value', lims=(0, mp/2, mp)),
                                    title='MNE', figure=1)

brain_PSF.add_foci([vertex], coords_as_verts=True, scale_factor=1.,
                   hemi='lh', color='green')
brain_PSF.add_foci([vert_max_psf], coords_as_verts=True, scale_factor=1.,
                   hemi='lh', color='black')

mlab.title('PSF sLORETA')

brain_CTF = stc_CTF.plot('sample', 'inflated', 'lh',
                                    subjects_dir=subjects_dir, 
                                    clim=dict(kind='value', lims=(0, mc/2, mc)),
                                    title='LOR', figure=2)

brain_CTF.add_foci([vertex], coords_as_verts=True, scale_factor=1.,
                   hemi='lh', color='green')
brain_CTF.add_foci([vert_max_ctf], coords_as_verts=True, scale_factor=1.,
                   hemi='lh', color='black')

mlab.title('CTF sLORETA')