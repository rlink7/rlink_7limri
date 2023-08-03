# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Histogram matching.
"""

# Imports
import numpy as np


def hist_matching(source, template, mask, plot=False):
    """ Adjust the pixel values of a grayscale image such that its histogram
    matches that of a target image.

    Parameters
    ----------
    source: np.ndarray
        the image to transform: the histogram is computed over the flattened
        array.
    template: np.ndarray
        the template image: same dimensions as the source image.
    mask: np.ndarray
        the mask image: same dimensions as the source image.
    plot: bool, default False
        plot the matched histograms.

    Returns
    -------
    matched: np.ndarray
        the transformed source image.
    """
    # Compute the source and template histograms
    mask_indices = np.where(mask == 1)
    h1, c1 = np.histogram(source[mask_indices], bins=65536)
    h1t, c1t = np.histogram(template[mask_indices], bins=65536)
    vec_size = len(h1)
    vec_size_f = float(vec_size)

    # Compute the normalized cumulated density functions
    cdf1 = h1.cumsum().astype(np.float64) / np.sum(h1)
    cdf2 = h1t.cumsum().astype(np.float64) / np.sum(h1t)

    # Compute the mapping
    M = np.zeros((vec_size, ))
    for idx in range(vec_size):
        ind = np.argmin(np.abs(cdf1[idx] - cdf2))
        M[idx] = ind
    M /= vec_size_f
    B = np.asarray(range(vec_size)) / vec_size_f

    # Interpolate the data by fitting a piecewise linear interpolant
    data1x = (source - c1.min()) / c1.max()
    data1x[np.where(data1x < 0)] = 0
    data1x[np.where(data1x > 1)] = 1
    indices = np.where((data1x > 0) & (data1x <= 1))
    xnew = data1x[indices]
    ynew = np.interp(xnew, B, M)

    # Plot if verbose
    if plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(B, M, 'k-')
        ax.plot(xnew[::10000], ynew[::10000], 'ro')
        plt.show()

    # Reshape the result
    matched = np.zeros(source.shape, dtype=source.dtype)
    matched[indices] = ynew * (c1t.max() - c1t.min()) + c1t.min()

    return matched


def _hist_matching(source, template, mask, plot=False):
    """ Adjust the pixel values of a grayscale image such that its histogram
    matches that of a target image.

    Parameters
    ----------
    source: np.ndarray
        the image to transform: the histogram is computed over the flattened
        array.
    template: np.ndarray
        the template image: same dimensions as the source image.
    mask: np.ndarray
        the mask image: same dimensions as the source image.
    plot: bool, default False
        plot the matched histograms.

    Returns
    -------
    matched: np.ndarray
        the transformed source image.
    """
    # Image information
    shape = source.shape
    dtype = source.dtype
    mask_indices = np.where(mask == 1)
    source = source[mask_indices]
    template = template[mask_indices]

    # Get the set of unique pixel values and their corresponding indices and
    # counts
    s_values, bin_idx, s_counts = np.unique(source, return_inverse=True,
                                            return_counts=True)
    t_values, t_counts = np.unique(template, return_counts=True)

    # Take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template images (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    # Interpolate linearly to find the pixel values in the template image
    # that correspond most closely to the quantiles in the source image
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)
    matched = np.zeros(shape, dtype=dtype)
    matched[mask_indices] = interp_t_values[bin_idx]

    return matched
