# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Workflows definition.
"""

import os
import limri
from .registration import li2mni, applytrf
from .maskeyes import li2mnieyes
from .normalization import li2mninorm


def li2mni_all(li_file, lianat_file, hanat_file, outdir, thr_factor=2,
               bins=300):
    """ Transform the Lithium (Li) data to the MNI space by using intermediate
    Hydrogene (H) data: l2mni, li2mnieyes, applytrf.

    Parameters
    ----------
    li_file: str
        path to the Li image.
    lianat_file: str
        path to the anat image acquired with the Li coil.
    hanat_file: str
        path of the anat image acquired with the H coil.
    outdir: str
        path to the destination folder.
    thr_factor: float, default 2
        multiply the mean of the second mode in the histogram to get a
        threshold to detect the eyes in the Lithium image.
    bins: int, default 300
        the number of bins in the histogram.
    """
    li2mni(li_file, lianat_file, hanat_file, outdir)
    li2mni_file = os.path.join(outdir, "li2mni.nii.gz")
    li2mnieyes(li2mni_file, outdir, thr_factor=2, bins=300)
    ref_file = os.path.join(os.path.dirname(limri.__file__), "resources",
                            "MNI152_T1_2mm.nii.gz")
    shiftedli2mni_file = os.path.join(outdir, "shiftedli2mni.nii.gz")
    transformlist = [
        os.path.join(outdir, "h2mni1Warp.nii.gz"),
        os.path.join(outdir, "h2mni0GenericAffine.mat"),
        os.path.join(outdir, "li2h0GenericAffine.mat"),
        os.path.join(outdir, "li2lianat0GenericAffine.mat")]
    applytrf(ref_file, li_file, transformlist, shiftedli2mni_file)
