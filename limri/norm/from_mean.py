# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
From mean send by Newcastle Team matching.
"""

# Imports
import numpy as np


def from_mean_matching(source, ref_val, concentration=2.):
    """ Adjust the pixel values of a grayscale image.

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
    # Normalize data
    matched = (source * concentration) / ref_val
    return matched
