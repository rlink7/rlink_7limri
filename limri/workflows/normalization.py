# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Intensity normalization.
"""

# Imports
import os
import nibabel
from limri.norm import hist_matching, minmax_matching
from limri.color_utils import print_title, print_subtitle, print_result

# Global parameters
NORM_MAP = {
    "hist": hist_matching,
    "minmax": minmax_matching
}


def li2mninorm(li2mni_file, li2mniref_file, mask_file, outdir, norm="hist"):
    """ Normalize intensities using histogram matching.

    Parameters
    ----------
    li2mni_file: str
        path to the Li image.
    li2mniref_file: str
        path to the reference Li image.
    mask_file: str
        the brain mask image.
    outdir: str
        path to the destination folder.
    norm: str, default 'hist'
        the normalization method.
    """
    print_title("Load data...")
    li2mni = nibabel.load(li2mni_file)
    li2mni_arr = li2mni.get_fdata()
    li2mniref = nibabel.load(li2mniref_file)
    li2mniref_arr = li2mniref.get_fdata()
    mask = nibabel.load(mask_file)
    mask_arr = mask.get_fdata()

    print_title("Normalization...")
    norm_fn = NORM_MAP.get(norm)
    if norm_fn is None:
        raise ValueError("Normalization method not defined.")
    norm_arr = norm_fn(li2mni_arr, li2mniref_arr, mask_arr)
    norm = nibabel.Nifti1Image(norm_arr, li2mni.affine)
    norm_file = os.path.join(outdir, "li2mninorm.nii.gz")
    nibabel.save(norm, norm_file)
    print_result(norm_file)
