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

# Imports
import os
import glob
import limri
from limri.normtools import fslreorient2std, fast, gzfile
from limri.regtools import antsregister, apply_transforms, apply_translation
from limri.color_utils import print_title, print_result, print_warning


def li2mni(li_file, lianat_file, hanat_file, outdir, li2lianat=None):
    """ Transform the Lithium (Li) data to the MNI space by using intermediate
    Hydrogene (H) data.

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
    li2lianat: 3-uplet, default None
        the translation applied on the Li image to compensate for different
        field of view between the Li and Li anat images (in mm).
    """
    print_title("Reorient images...")
    lianat_reo_file = os.path.join(outdir, "lianat.nii.gz")
    if not os.path.isfile(lianat_reo_file):
        gzfile(lianat_file, lianat_reo_file)
        fslreorient2std(lianat_reo_file, lianat_reo_file, save_trf=True)
    else:
        print_warning("lianat already reoriented")
    print_result(lianat_reo_file)
    hanat_reo_file = os.path.join(outdir, "hanat.nii.gz")
    if not os.path.isfile(hanat_reo_file):
        gzfile(hanat_file, hanat_reo_file)
        fslreorient2std(hanat_reo_file, hanat_reo_file, save_trf=True)
    else:
        print_warning("hanat already reoriented")
    print_result(hanat_reo_file)

    print_title("Bias field correction...")
    lianat_bcorr_file = os.path.join(outdir, "lianat_restore.nii.gz")
    if not os.path.isfile(lianat_bcorr_file):
        lianat_bcorr_file = fast(
            lianat_reo_file, lianat_reo_file.replace(".nii.gz", ""))
    else:
        print_warning("lianat already bias corrected")
    print_result(lianat_bcorr_file)
    hanat_bcorr_file = os.path.join(outdir, "hanat_restore.nii.gz")
    if not os.path.isfile(hanat_bcorr_file):
        hanat_bcorr_file = fast(
            hanat_reo_file, hanat_reo_file.replace(".nii.gz", ""))
    else:
        print_warning("hanat already bias corrected")
    for key1 in ("lianat", "hanat"):
        for key2 in ("pve", "mixeltype", "seg"):
            regex = os.path.join(outdir, f"{key1}_{key2}*.nii.gz")
            for path in glob.glob(regex):
                os.remove(path)
    print_result(hanat_bcorr_file)

    print_title("Coregistration & normalization...")
    rigid_transforms = [os.path.join(outdir, "li2h0GenericAffine.mat")]
    deform_transforms = [
        os.path.join(outdir, "h2mni1Warp.nii.gz"),
        os.path.join(outdir, "h2mni0GenericAffine.mat")]
    is_generated = True
    for path in rigid_transforms + deform_transforms:
        if not os.path.isfile(path):
            is_generated = False
            break
    ref_file = os.path.join(os.path.dirname(limri.__file__), "resources",
                            "MNI152_T1_2mm.nii.gz")
    if not is_generated:
        antsregister(
            template_file=ref_file, lianat_file=lianat_bcorr_file,
            hanat_file=hanat_bcorr_file, outdir=outdir)
    else:
        print_warning("li2mni transformation already computed")
    print_result(deform_transforms + rigid_transforms)

    print_title("Li image to MNI space...")
    li2mni_file = os.path.join(outdir, "li2mni.nii.gz")
    if not os.path.isfile(li2mni_file):
        li2lianat = li2lianat or (0, 0, 0)
        apply_translation(image_file=li_file, translation=li2lianat,
                          filename=li2mni_file)
        apply_transforms(
            fixed_file=ref_file, moving_file=li2mni_file,
            transformlist=deform_transforms + rigid_transforms,
            filename=li2mni_file)
    else:
        print_warning("li2mni transformation already applied")
    print_result(li2mni_file)
