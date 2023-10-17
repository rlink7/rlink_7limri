# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
MinMax matching.
"""

# Imports
import numpy as np
from sklearn import mixture


def minmax_matching(source, template, mask, concentration=2.):
    """ Adjust the pixel values of a grayscale image such that its dynamic
    matches that of a target image.

    Parameters
    ----------
    source: np.ndarray
        the image to transform.
    template: np.ndarray
        the template image: same dimensions as the source image. In this case
        the template is a phantom with 1 compartment.
    mask: np.ndarray
        the mask image: same dimensions as the source image.
    concentration: float, default 2.
        the compartment concentration in milli mols / litre.

    Returns
    -------
    matched: np.ndarray
        the transformed source image.
    """
    # Segment the compartment
    clf = mixture.GaussianMixture(n_components=2, covariance_type="full")
    clf.fit(template.reshape(-1, 1))
    m1, m2 = clf.means_
    thr = (m1 + m2) / 2.
    template[template < thr] = 0.

    # Get the reference value
    ref_val = np.mean(template[template > 0])

    # Normalize data
    matched = (source * concentration) / ref_val

    return matched


def norm(source, ref_val, concentration=2.):
    """ Simply normalize the intensities from a reference intensity value and
    the associated compartment concentration.

    Parameters
    ----------
    source: np.ndarray
        the image to transform.
    ref_val: int
        reference value of phantom intensity for the corresponding site.
    concentration: float, default 2.
        the compartment concentration in milli mols / litre.

    Returns
    -------
    matched: np.ndarray
        the transformed source image.
    """
    matched = (source * concentration) / ref_val
    return matched
