# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Eyes detection workflows definition.
"""

# Imports
import os
import nibabel
import numpy as np
from skimage import measure
from sklearn import mixture
from scipy.stats import norm
from scipy import ndimage
from collections import Counter
import matplotlib.pyplot as plt
import limri
from limri.regtools import apply_translation
from limri.color_utils import print_title, print_subtitle, print_result


def li2mnieyes(li2mni_file, outdir, bins=300):
    """ Detect the eyes in a Lithium MRI image in the MNI space and determine
    a potential shift as a translation.

    Parameters
    ----------
    li2mni_file: str
        path to the Li image.
    outdir: str
        path to the destination folder.
    bins: int, default 300
        the number of bins in the histogram.
    """
    print_title("Load data...")
    im = nibabel.load(li2mni_file)
    arr = im.get_fdata()
    ref_file = os.path.join(os.path.dirname(limri.__file__), "resources",
                            "MNI152_T1_2mm_eye_mask.nii.gz")
    ref_im = nibabel.load(ref_file)
    ref_arr = ref_im.get_fdata()

    print_title("Last peak extraction: GMM...")
    data = arr[arr > 0]
    data.shape += (1, )
    mode = get_last_mode(data, bins=bins, snapdir=outdir)
    print_result(f"last mode: {mode}")

    print_title("Extract eyes...")
    arr[arr < 2 * mode] = 0
    arr[arr > 0] = 1
    mask_im = nibabel.Nifti1Image(arr, im.affine)
    nibabel.save(mask_im, os.path.join(outdir, "li2mnieyes.nii.gz"))
    arr = ndimage.binary_opening(arr, iterations=3)
    li_labels = measure.label(arr, background=0)
    label_im = nibabel.Nifti1Image(li_labels, im.affine)
    nibabel.save(label_im, os.path.join(outdir, "li2mnilabels.nii.gz"))
    assert li_labels.max() >= 2, "Assume at least 2 CCs."
    count = Counter(li_labels[li_labels != 0])
    li_centroids = []
    for key in list(count.keys())[:2]:
        li_centroids.append(
            np.asarray(np.where(li_labels == key)).mean(axis=1))
    li_centroids = np.array(li_centroids)
    print_result(f"li eyes centroids: {li_centroids}")
    ref_arr = ndimage.binary_erosion(ref_arr, iterations=5).astype(int)
    ref_labels = measure.label(ref_arr, background=0)
    assert ref_labels.max() >= 2, "Assume at least 2 CCs."
    count = Counter(ref_labels[ref_labels != 0])
    ref_centroids = []
    for key in list(count.keys())[:2]:
        ref_centroids.append(
            np.asarray(np.where(ref_labels == key)).mean(axis=1))
    ref_centroids = np.array(ref_centroids)
    print_result(f"ref eyes centroids: {ref_centroids}")

    print_title("Compute translation from barycenters...")
    li_bary = np.mean(li_centroids, axis=0)
    ref_bary = np.mean(ref_centroids, axis=0)
    li2ref_translation = ref_bary - li_bary
    print_result(f"li2ref estimated translation: {li2ref_translation}")

    print_title("Apply translation...")
    shiftedli2mni_file = os.path.join(outdir, "shiftedli2mni.nii.gz")
    apply_translation(image_file=li2mni_file, translation=li2ref_translation,
                      filename=shiftedli2mni_file)
    print_result(shiftedli2mni_file)


def get_last_mode(data, bins=300, snapdir=None):
    """ Grabs the last peak or shoulder.

    Parameters
    ----------
    data: array
        the values of the histogram. See density and weights for a
    bins: int, default 300
        the number of bins in the histogram.
    snapdir: str, default None
        a folder where QC images will be generated.

    Returns
    -------
    last_mode: int
        the last mode in the histogram.
    """
    clf = mixture.GaussianMixture(n_components=2, covariance_type="full")
    clf.fit(data)
    m1, m2 = clf.means_
    w1, w2 = clf.weights_
    c1, c2 = clf.covariances_
    last_mode = max(clf.means_)
    if snapdir is not None:
        fig = plt.figure()
        x = data.copy().ravel()
        x.sort()
        gauss1 = w1 * norm.pdf(x, m1, np.sqrt(c1)).ravel()
        gauss2 = w2 * norm.pdf(x, m2, np.sqrt(c2)).ravel()
        plt.hist(data, bins=bins, density=True, histtype="step", facecolor="g",
                 alpha=0.5, lw=2)
        plt.plot(x, gauss1, c="C0")
        plt.plot(x, gauss2, c="C1")
        plt.plot(x, gauss1 + gauss2, lw=3, c="C2", ls="dashed")
        plt.axvline(x=last_mode, c="r")
        plt.title("Last peak")
        filename = os.path.join(snapdir, "last_peak.png")
        fig.savefig(filename)
    return last_mode
