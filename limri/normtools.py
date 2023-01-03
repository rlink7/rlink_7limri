# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2022
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Normalization tools.
"""

# Imports
import os
import subprocess
import numpy as np
import nibabel
from .regtools import flirt2aff


def gzfile(input_image, output_image):
    """ Gzip a file if necessary.

    Parameters
    ----------
    input_image: str
        the image to reorient.
    output_image: str
        the gzip image.

    Returns
    -------
    gzip_image: str
        the gzip file.
    """
    im = nibabel.load(input_image)
    nibabel.save(im, output_image)
    return output_image


def fslreorient2std(input_image, output_image, save_trf=True):
    """ Reorient an image to match the approximate orientation of the standard
    template image (MNI152).

    Parameters
    ----------
    input_image: str
        the image to reorient.
    output_image: str
        the reoriented image.
    save_trf: bool, default True
        opttionally save the reorientation matrix.
    """
    cmd1 = ["fslreorient2std", input_image, output_image]
    cmd2 = ["fslreorient2std", input_image]
    subprocess.check_call(cmd1)
    if save_trf:
        stdout = subprocess.check_output(cmd2)
        fsl_trf_file = output_image.split(".")[0] + ".fsl.trf"
        with open(fsl_trf_file, "wt") as open_file:
            open_file.write(stdout.decode("utf8"))
        trf_file = output_image.split(".")[0] + ".trf"
        np.savetxt(trf_file, flirt2aff(fsl_trf_file, output_image,
                                       input_image))
    return output_image


def fast(input_file, out_fileroot, klass=3, im_type=1, segments=False,
         bias_field=True, bias_corrected_im=True, probabilities=False):
    """ FAST (FMRIB's Automated Segmentation Tool) segments a 3D image of
    the brain into different tissue types (Grey Matter, White Matter, CSF,
    etc.), whilst also correcting for spatial intensity variations (also
    known as bias field or RF inhomogeneities).

    Parameters
    ----------
    input_file: str
        the image to be segmented.
    out_fileroot: str
        output basename.
    klass: int, default 3
        number of tissue-type classes.
    im_type: int, default 1
        type of image 1=T1, 2=T2, 3=PD.
    segments: bool, default False
        outputs a separate binary image for each tissue type.
    bias_field, default True
        output estimated bias field.
    bias_corrected_im: bool, default True
        output bias-corrected image.
    probabilities: bool, default False
        outputs individual probability maps.

    Returns
    -------
    biascorrected_file: str
        the bias corrected input image.
    """
    bool_params = {
        "-g": segments,
        "-b": bias_field,
        "-B": bias_corrected_im,
        "-p": probabilities
    }
    cmd = ["fast", "-o", out_fileroot, "-n", str(klass), "-t", str(im_type)]
    for name, value in bool_params.items():
        if value:
            cmd.append(name)
    cmd.append(input_file)
    subprocess.check_call(cmd)
    image_ext = ".nii.gz"
    biascorrected_file = out_fileroot + "_restore" + image_ext
    if not os.path.isfile(biascorrected_file):
        biascorrected_file = None
    return biascorrected_file
