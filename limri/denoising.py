# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Denoising tools.
"""

# Import
from dipy.denoise.nlmeans import nlmeans
from dipy.denoise.noise_estimate import estimate_sigma


def nlm_denoising(arr, n_coils=0):
    """ Non-local means for denoising 3D and 4D images, using blockwise
    averaging approach.

    Parameters
    ----------
    arr: 3D or 4D ndarray
        The array to be denoised
    n_coils: int, default 0
        Number of coils of the receiver array. Use N = 1 in case of a SENSE
        reconstruction (Philips scanners) or the number of coils for a GRAPPA
        reconstruction (Siemens and GE). Use 0 to disable the correction
        factor, as for example if the noise is Gaussian distributed.

    Returns
    -------
    denoised_arr: ndarray
        the denoised ``arr`` which has the same shape as ``arr``.

    References
    ----------
    .. [Coupe08] P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C.
                 Barillot, An Optimized Blockwise Non Local Means Denoising
                 Filter for 3D Magnetic Resonance Images, IEEE Transactions on
                 Medical Imaging, 27(4):425-441, 2008
    .. [Coupe11] Pierrick Coupe, Jose Manjon, Montserrat Robles, Louis Collins.
                Adaptive Multiresolution Non-Local Means Filter for 3D MR Image
                Denoising IET Image Processing, Institution of Engineering and
                Technology, 2011
    """
    sigma = estimate_sigma(arr, N=n_coils)
    denoised_arr = nlmeans(arr, sigma=sigma, patch_radius=1, block_radius=2,
                           rician=True)
    return denoised_arr
