"""
=========================================================
Functions for EEG/MEG Resolution
=========================================================
"""
# OH July 2018

import glob

from copy import deepcopy

import numpy as np

from scipy.linalg import svd as sci_svd

import mne

import mne.minimum_norm.WH_psf_ctf as psf_ctf


# TO DO: Still need to find good example for source distribution from which to compute null space.
def leadfield_nullspace(fwd, noise_cov, regpar=0.):
    """ Compute null space distributions for leadfield.

    Parameters
    ----------
    fwd: forward solution
        Used to get leadfield matrix.
    noise_cov: noise covariance matrix
        Used to (diagonally) whiten leadfield.
    regpar: float
        Regularisation parameter (Tiknonov) for pseudoinverse.
        Default: 0.

    Returns
    -------
    nullspace: 2D numpy array (n_ch distributions).
        Distributions in leadfield nullspace.
        Null projector will be applied to point sources.
    """

    # get leadfield matrix from forward solution
    leadfield = psf_ctf._pick_leadfield(fwd['sol']['data'], fwd, fwd['info']['ch_names'])

    # number of channels and sources
    n_c, n_s = leadfield.shape

    # covariance matrix will be modified
    cov = deepcopy(noise_cov)

    # use diagonal noise covariance matrix
    cov = cov.as_diag()

    # compute whitener from noise covariance matrix
    info = fwd['info']
    info['projs'] = noise_cov['projs']
    info['comps'] = ''
    whitener = mne.cov.compute_whitener(cov, info)[0]

    # whiten leadfield
    leadfield = whitener.dot(leadfield)

    # scale leadfield, maybe better for stability
    leadfield = leadfield / leadfield.max()

    # compute pseudoinverse the "matrix way", for intuitive regularisation

    # Gram matrix
    gram = leadfield.dot(leadfield.T)

    trace_gram = np.trace(gram)

    eye_gram = np.eye(n_c)

    # invert Gram matrix using Tikhonov regularisation
    # (don't use np.linalg.inv here)
    gram_inv = np.linalg.pinv(gram + regpar*(trace_gram/n_c)*eye_gram)

    # full pseudoinverse
    pseudoinv = leadfield.T.dot(gram_inv)

    # null space projector
    null_proj = np.eye(n_s) - pseudoinv.dot(leadfield)

    sources = np.empty([n_s,5])

    # point sources, every tenth
    sources = np.eye(n_s)[:,::10]

    # null-space projection applied to sources
    nullspace = null_proj.dot(sources)
    
    return nullspace, leadfield, pseudoinv, null_proj


def make_resolution_matrix(fwd, invop, method, lambda2):
    """ Compute resolution matrix for linear inverse operator.

    Parameters
    ----------
    fwd: forward solution
        Used to get leadfield matrix.
    invop: inverse operator
        Inverse operator to get inverse matrix.
        pick_ori='normal' will be selected.
    method: string
        Inverse method to use (MNE, dSPM, sLORETA).
    lambda2: float
        The regularisation parameter.

    Returns
    -------
        resmat: 2D numpy array.
        Resolution matrix (inverse matrix times leadfield).
    """

    info = fwd['info']

    # get leadfield matrix from forward solution
    leadfield = psf_ctf._pick_leadfield(fwd['sol']['data'], fwd, fwd['info']['ch_names'])


    invmat, _ = psf_ctf._get_matrix_from_inverse_operator(invop, fwd, method=method, lambda2=lambda2,
                                                pick_ori='normal')

    resmat = invmat.dot(leadfield)

    return resmat


def make_resolution_matrix_lcmv(fwd, info, noise_cov, data_cov):
    """ Compute resolution matrix for linear inverse operator.

    Parameters
    ----------
    forward : dict
        The forward operator.
    info: instance of Info
        Should contain measurement information, e.g. sfreq, projs.
    noise_cov: noise covariance matrix
        Used to compute whitener. Should be regularised.
    data_cov: data covariance matrix
        Used to compute LCMV beamformer. Should be regularised.

    Returns
    -------
        resmat: 2D numpy array.
        Resolution matrix (inverse matrix times leadfield).
    """    

    # get leadfield matrix from forward solution
    leadfield = psf_ctf._pick_leadfield(fwd['sol']['data'], fwd, info['ch_names'])


    invmat = _get_matrix_from_LCMV_beamformer(fwd, info, noise_cov, noise_cov)

    resmat = invmat.dot(leadfield)

    return resmat



def mat_to_stc(mat, src):
    """ Turn columns of matrix to STCs.

    Parameters
    ----------
    mat: 2D numpy array
        For example resolution matrix (or its transpose).
    src: SourceSpace
        Source space with vertex numbers, e.g. forward['src'].

    Returns
    -------
    stc: Instance of SourceEstimate
        PSFs or CTFs as source estimate object.
    """

    vertno = [src[0]['vertno'], src[1]['vertno']]

    # convert norms to source estimate
    stc = mne.SourceEstimate(mat, vertno, tmin=0., tstep=1.)

    return stc


def relative_amplitude(resmat, locations, axis, metric='peak'):
    """ Compute relative amplitude metrics for resolution matrix.

    Parameters
    ----------
    resmat: 2D numpy array
        The resolution matrix (nloc-by-nloc).
    locations: 2D (nloc-by-3) numpy array
        Locations (in m) to be used for resolution metrics (distances etc.).
    axis: integer (0 or 1)
        Whether to compute metrics for columns (=0, PSFs) or rows (=1, CTFs).
    metric: string ('peak')
        Which amplitudes to use.
        'peak': Ratio between absolute maximum amplitudes of peaks per location
                and maximum peak across locations.
        'sum': Ratio between sums of absolute amplitudes.

    Returns
    -------
        relamp: 1D numpy array.
        Relative amplitude metric per location.
    """

    # NOTE: locations needed?

    # only use absolute values
    resmat = np.absolute(resmat)

    # Ratio between amplitude at peak and global peak maximum
    if metric.lower() == 'peak':

        # maximum amplitudes along specified axis
        maxamps = resmat.max(axis=axis)

        # global absolute maximum
        maxmaxamps = maxamps.max()

        relamp = maxamps / maxmaxamps

    # ratio between sums of absolute amplitudes
    elif metric.lower() == 'sum':

        # sum of amplitudes per location
        sumamps = np.sum(resmat, axis=axis)

        # maximum of summed amplitudes
        sumampsmax = sumamps.max()

        relamp = sumamps / sumampsmax
       
    return relamp


def spatial_width(resmat, locations, axis, metric='sd'):
    """ Compute spatial width metrics for resolution matrix.

    Parameters
    ----------
    resmat: 2D numpy array
        The resolution matrix (nloc-by-nloc).
    locations: 2D (nloc-by-3) numpy array
        Locations (in m) to be used for resolution metrics (distances etc.).
    axis: integer (0 or 1)
        Whether to compute metrics for columns (=0, PSFs) or rows (=1, CTFs).
    metric: string ('sd' | 'rad')
        What type of width metric to compute.
        'sd': spatial deviation (e.g. Molins et al.).
        'maxrad': maximum radius to 50% of max amplitude.

    Returns
    -------
        width: 1D numpy array.
        Spatial width metric per location.
    """

    # only use absolute values
    resmat = np.absolute(resmat)

    # The below will operate on columns
    if axis == 1:

        resmat = resmat.T

    # find indices of maxima along rows
    resmax = resmat.argmax(axis=0)

    # initialise output array
    width = np.empty(len(resmax))

    # spatial deviation as in Molins et al.
    if metric.lower() == 'sd':

        for ii in range(0, locations.shape[0]):

            # locations relative to true source
            diffloc = locations - locations[ii,:]

            # squared Euclidean distances to true source
            locerr = np.sum(diffloc**2,1)

            # pick current row
            resvec = resmat[:,ii]**2

            # spatial deviation (Molins et al, NI 2008, eq. 12)
            width[ii] = np.sqrt(np.sum(np.multiply(locerr, resvec))/np.sum(resvec))

    # maximum radius to 50% of max amplitude
    elif metric.lower() == 'maxrad':

        # peak amplitudes per location across columns
        maxamp = resmat.max(axis=0)

        for (ii,aa) in enumerate(maxamp): # for all locations

            # pick current column
            resvec = resmat[:,ii]

            # indices of elements where values are larger than 50% of peak amplitude
            amps50idx = np.where(resvec > 0.5*aa)[0]

            # get distances for those indices from true source position
            locs50 = locations[amps50idx,:] - locations[ii,:]

            # get maximum distance
            width[ii] = np.sqrt(np.sum(locs50**2, 1).max())
       
    return width


def localisation_error(resmat, locations, axis, metric='peak'):
    """ Compute localisation error metrics for resolution matrix.

    Parameters
    ----------
    resmat: 2D numpy array
        The resolution matrix (nloc-by-nloc).
    locations: 2D (nloc-by-3) numpy array
        Locations (in m) to be used for resolution metrics (distances etc.).
    axis: integer (0 or 1)
        Whether to compute metrics for columns (=0, PSFs) or rows (=1, CTFs).
    metric: string ('peak')
        What type of localisation error to compute.
        'peak': peak localisation error, Euclidean distance.

    Returns
    -------
        locerr: 1D numpy array.
        Localisation error per location (m).
    """

    # only use absolute values
    resmat = np.absolute(resmat)

    # The below will operate on columns
    if axis == 1:

        resmat = resmat.T    

    # Euclidean distance between true location and maximum
    if metric.lower() == 'peak':

        # find indices of maxima along columns
        resmax = resmat.argmax(axis=0)

        # locations of maxima
        maxloc = locations[resmax,:]

        # difference between locations of maxima and true locations
        diffloc = locations - maxloc

        # Euclidean distance
        locerr = np.sqrt(np.sum(diffloc**2,1))

    # centre of gravity
    elif metric.lower() == 'cog':

        # initialise result array
        locerr = np.empty(locations.shape[0])

        for (ii, rr) in enumerate(locations): # for every vertex

            # differences to true location
            difflocs = locations - rr

            # corresponding column of resmat
            resvec = resmat[:,ii].T

            # centre of gravity
            cog = resvec.dot(locations) / np.sum(resvec)

            # centre of gravity
            locerr[ii] = np.sqrt( np.sum( (rr - cog)**2 ) )

    return locerr


def sensitivity_map(fwd, noise_cov, diag=True, metric='norm', maxnorm=False):
    """ Compute sensitivity maps for EEG/MEG (norms of leadfield columns).

    Parameters
    ----------
    fwd: forward solution
        Used to get leadfield matrix.
    noise_cov: noise covariance matrix
        Used to whiten leadfield. Should already be regularised.
        Diagonal will be used.
    diag: Boolean
        Whether to use only the diagonal (True) or whole
        matrix for whitening. Default: True (diagonal).
    metric: string
        Whether 'SNR' (Goldenholz) or 'RMS' of columns to be computed.
        Default 'SNR'.
    maxnorm: Boolean
        Whether to normalise sensitivity map to maximum or not.

    Returns
    -------
        stc_metric: Source estimate.
        Distribution of leadfield column metric.
    """

    info = fwd['info']

    # covariance matrix will be modified
    cov = deepcopy(noise_cov)

    # use diagonal noise covariance matrix
    cov = cov.as_diag()

    whitener = mne.cov.compute_whitener(cov, info)[0]
    
    print('###\nGetting leadfield, computing column norms.\n###')

    vertno = [fwd['src'][0]['vertno'], fwd['src'][1]['vertno']]

    leadfield = psf_ctf._pick_leadfield(fwd['sol']['data'], fwd, fwd['info']['ch_names'])

    print('Leadfield has dimensions (%d, %d).' % leadfield.shape)

    print('Whitening leadfield.')
    lfd_white = whitener.dot(leadfield)

    if metric.lower() == 'snr':
        # SNR in decibel units, as in Goldenholz et al., HBM 2009 (eq. 1)
        ldf_metric = 10.*np.log10(np.average(lfd_white**2, axis=0))

    elif metric.lower() == 'rms':

        # compute norm per column
        ldf_metric = np.sqrt(np.average(lfd_white**2, axis=0))

    # if specified, normalise metrics to absolute maximum across vertices
    if maxnorm:

        ldf_metric = ldf_metric / np.abs(ldf_metric).max()

    # fake multiple time steps
    ldf_metric_rep = np.repeat(ldf_metric[:,np.newaxis],5, axis=1)

    # convert norms to source estimate
    stc_metric = mne.SourceEstimate(ldf_metric_rep, vertno, tmin=0., tstep=1.)

    return stc_metric


def leadfield_svd(fwd, noise_cov, diag=True, norm=True):
    """ Compute singular values of whitened leadfield matrix.

    Parameters
    ----------
    fwd: forward solution
        Used to get leadfield matrix.
    noise_cov: noise covariance matrix
        Used to whiten leadfield. Should already be regularised.
        Diagonal will be used.
    diag: Boolean
        Whether to use only the diagonal (True) or whole
        matrix for whitening. Default: True (diagonal).
    norm: Boolean
        Whether to normalise to maximum singular value or not.
        Default: True.

    Returns
    -------
        sing_vals: numpy array.
        Singular values of leadfield.
    """

    info = fwd['info']

    # covariance matrix will be modified
    cov = deepcopy(noise_cov)

    # use diagonal noise covariance matrix
    cov = cov.as_diag()

    whitener = mne.cov.compute_whitener(cov, info)[0]
    
    print('###\nGetting and whitening leadfield.\n###')

    vertno = [fwd['src'][0]['vertno'], fwd['src'][1]['vertno']]

    leadfield = psf_ctf._pick_leadfield(fwd['sol']['data'], fwd, fwd['info']['ch_names'])

    print('Leadfield has dimensions (%d, %d).' % leadfield.shape)

    print('Whitening leadfield.')
    lfd_white = whitener.dot(leadfield)

    print('Compute SVD of whitened leadfield.')
    sing_vals = sci_svd(lfd_white, compute_uv=False)

    if norm: # normalise to maximum

        sing_vals = sing_vals / sing_vals[0]

    return sing_vals


def normalise_stc(stc):
    """ Normalise data in STC object to absolute maximum.

    Parameters
    ----------
    stc: SourceEstimate
        The data to normalise.

    Returns
    -------
        stc_norm: Source estimate.
        STC with data normalised to absolute maximum.
    """

    data = stc.data

    data = data / np.absolute(data).max()

    # convert normalised data to source estimate
    stc_norm = mne.SourceEstimate(data, stc.vertices, tmin=stc.tmin, tstep=stc.tstep)

    return stc_norm


def find_nearest_vertex(src, vertex, hemi=0):
    """ Find nearest used vertex (e.g. fwd, invop) from Freesurfer surface.

    Parameters
    ----------
    src: Source space
        The source space with vertex information.
    vertex: int
        The vertex to look for (e.g. from mne_analyze).
    hemi: 0 | 1
        Hemisphere from which vertex to choose.
        0: LH, 1: RH.

    Returns
    -------
    nearest_vert: int.
        Nearest used vertex.
    """

    # coordinate of target vertex
    loc_vertex = src[hemi]['rr'][vertex]

    # vertices used in source estimate
    verts_used = src[hemi]['vertno']

    # coordinates of used vertices
    locs_used = src[hemi]['rr'][verts_used]

    # difference to target vertex
    locs_diff = locs_used - loc_vertex

    # squared distance
    locs_dist = np.sum(locs_diff**2, axis=1)

    min_idx = locs_dist.argmin()

    nearest_vert = verts_used[min_idx]

    return nearest_vert



def _get_matrix_from_LCMV_beamformer(forward, info, noise_cov, data_cov):
    """Get inverse matrix for LCMV beamformer.    
    Parameters
    ----------
    forward : dict
        The forward operator.
    info: instance of Info
        Should contain measurement information, e.g. sfreq.
    noise_cov: noise covariance matrix
        Used to compute whitener. Should be regularised.
    data_cov: data covariance matrix
        Used to compute LCMV beamformer. Should be regularised.
    
    Returns
    -------
    invmat : ndarray
        Inverse matrix associated with LCMV beamformer.
    """
    # based on _get_matrix_from_inverse_operator() from psf_ctf module.

    # number of channels for identity matrix
    n_chs = len(info['ch_names'])

    # create identity matrix as input for inverse operator
    # set elements to zero for non-selected channels
    id_mat = np.eye(n_chs)

    # convert identity matrix to evoked data type (pretending it's an epoch)
    evo_ident = mne.EvokedArray(id_mat, info=info, tmin=0.)

    ### apply beamformer to identity matrix
    stc_lcmv = mne.beamformer.lcmv(evo_ident, forward, noise_cov, data_cov,
                                    max_ori_out='signed')

    # turn source estimate into numpsy array
    invmat = stc_lcmv.data

    return invmat