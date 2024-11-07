# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Registration tools.
"""

# Imports
import os
import subprocess
import numpy as np
import nibabel
import scipy.io as sio
from limri.color_utils import print_subtitle, print_result


def flirt(in_file, ref_file, out, omat=None, init=None, cost="corratio",
          usesqform=False, displayinit=False, anglerep="euler", bins=256,
          interp="trilinear", dof=12, applyxfm=False, applyisoxfm=None,
          nosearch=False, wmseg=None, verbose=0):
    """ Affine registration with FSL flirt.

    Parameters
    ----------
    in_file: str
        input volume.
    ref_file: str
        reference volume.
    out: str
        output volume
    omat: str, default None
        matrix filename - output in 4x4 ascii format.
    init: str, default None
        input 4x4 affine matrix
    cost: str, default 'corratio'
        choose the most appropriate option: 'mutualinfo', 'corratio',
        'normcorr', 'normmi', 'leastsq', 'labeldiff', 'bbr'.
    usesqform: bool, default False
        initialise using appropriate sform or qform.
    displayinit: bool, default False
        display initial matrix.
    anglerep: str, default 'euler'
        choose the most appropriate option: 'quaternion', 'euler'.
    bins: int, default 256
        number of histogram bins
    interp: str, default 'trilinear'
        choose the most appropriate option: 'trilinear', 'nearestneighbour',
        'sinc', 'spline'.
    dof: int, default 12
        number of transform dofs.
    applyxfm: bool, default False
        ppplies transform (no optimisation) - requires -init.
    applyisoxfm: float, default None
        the integer defines the scale - as applyxfm but forces isotropic
        resampling.
    verbose: int, default 0
        0 is least and default.
    nosearch: bool, default False
        if set perform no search to initialize the optimization.
    wmseg: str, default None
        white matter segmentation volume needed by BBR cost function.

    Returns
    -------
    out: str
        output volume.
    omat: str
        output matrix filename - output in 4x4 ascii format.
    """
    cmd = ["flirt",
           "-in", in_file,
           "-ref", ref_file,
           "-cost", cost,
           "-searchcost", cost,
           "-anglerep", anglerep,
           "-bins", str(bins),
           "-interp", interp,
           "-dof", str(dof),
           "-out", out,
           "-verbose", str(verbose)]
    if usesqform:
        cmd += ["-usesqform"]
    if displayinit:
        cmd += ["-displayinit"]
    if applyxfm:
        cmd += ["-applyxfm"]
    if nosearch:
        cmd += ["-nosearch"]
    if init is not None:
        cmd += ["-init", init]
    if applyisoxfm is not None:
        cmd += ["-applyisoxfm", str(applyisoxfm)]
    if cost == "bbr":
        cmd += ["-wmseg", wmseg]
    if not applyxfm:
        cmd += ["-omat", omat]
    subprocess.check_call(cmd)
    return out, omat


def applywarp(in_file, ref_file, out_file, warp_file, pre_affine_file=None,
              post_affine_file=None, interp="trilinear", verbose=0):
    """ Apply SPM deformation field.

    Parameters
    ----------
    in_file: str
        filename of input image (to be warped).
    ref_file: str
        filename for reference image.
    out_file: str
        filename for output (warped) image.
    warp_file: str
        filename for warp/coefficient (volume).
    pre_affine_file: str, default None
        filename for pre-transform (affine matrix).
    post_affine_file: str, default None
        filename for post-transform (affine matrix).
    interp: str, default 'trilinear'
        interpolation method {nn,trilinear,sinc,spline}
    verbose: int, default 0
        the verbosity level.
    """
    cmd = ["applywarp",
           "-i", in_file,
           "-r", ref_file,
           "-o", out_file,
           "-w", warp_file,
           "--abs",
           "--interp={0}".format(interp),
           "--verbose={0}".format(verbose)]
    if pre_affine_file is not None:
        cmd.append("--premat={0}".format(pre_affine_file))
    if post_affine_file is not None:
        cmd.append("--postmat={0}".format(post_affine_file))
    subprocess.check_call(cmd)


def flirt2aff(mat_file, in_file, ref_file):
    """ Map from 'in_file' image voxels to 'ref_file' voxels given `mat_file`
    omat FSL affine transformation.

    Parameters
    ------------
    mat_file: str
        filename of output '-omat' transformation file from FSL flirt.
    in_file: str
        filename of the image passed to flirt as the '-in' image.
    ref_file: str
        filename of the image passed to flirt as the '-ref' image.

    Returns
    -------
    omat: array (4, 4)
        array containing the transform from voxel coordinates in image
        for 'in_file' to voxel coordinates in image for 'ref_file'.
    """
    flirt_affine = np.loadtxt(mat_file)
    in_img = nibabel.load(in_file)
    ref_img = nibabel.load(ref_file)
    in_hdr = in_img.header
    ref_hdr = ref_img.header

    def _x_flipper(n):
        flipr = np.diag([-1, 1, 1, 1])
        flipr[0, 3] = n - 1
        return flipr

    inspace = np.diag(in_hdr.get_zooms()[:3] + (1, ))
    refspace = np.diag(ref_hdr.get_zooms()[:3] + (1, ))
    if np.linalg.det(in_img.affine) >= 0:
        inspace = np.dot(inspace, _x_flipper(in_hdr.get_data_shape()[0]))
    if np.linalg.det(ref_img.affine) >= 0:
        refspace = np.dot(refspace, _x_flipper(ref_hdr.get_data_shape()[0]))

    omat = np.dot(np.linalg.inv(refspace), np.dot(flirt_affine, inspace))

    return omat


def normalize2field():
    """ Transform a SPM normalization field to a FSL field.

    For the deformation fields generated in SPM, the values in each of the
    three components (ie (:,:,:,1,1), (:,:,:,1,2) and (:,:,:,1,3) ) encode
    the x,y,z coordinates of the corresponding location (in units of mm).

    SPM conventions:
    1) The file encodes the mapping rather than displacement fields for a
    couple of reasons:
    a) Even in 2012, many people still seem to think that displacements
    can be combined by addition or subtraction, rather than by composing
    or inverting the mappings.  (Attempting to compose deformations by
    adding the displacement fields would be analogous to saying 2*2=3 or
    3*3=5). Using the mapping instead of the displacements would (I hope)
    make this abuse a bit less likely.
    b) It saves needing to encode an additional affine transform matrix
    for dealing with the identity transform to which the displacements are
    added. This makes life simpler.
    2) It encodes coordinates in mm, rather than voxels. If the mapping
    was to voxels, it would not work so well for images that are in
    alignment via the matrices encoded by the sform or qform fields.
    Also, if the mapping is to voxels, the values would depend on whether
    the first voxel in the file is denoted as 1,1,1 (as in MATLAB, Fortran
    or SPM) or as 0,0,0 (as in C or FSL). The mm coordinates are
    unambiguous - providing the sform (or sform) fields of the image that
    the deformation points to is filled in.

    The bad news here is that you do need to use the sform (qform)
    information to convert from SPM to FSL conventions. In FSL we can use
    either absolute coordinates (as is used in SPM) or relative ones.
    However, our mm coordinate system is not the same, so you'll need to map
    between FSL's and SPM's if you want to get the warp correct.

    The problem you are experiencing will be about coordinate system
    conventions. Our mm coordinate system is just a scalar multiple of the
    voxel coordinates if the image has a negative sform (or qform)
    determinant. If the determinant of the sform (or qform) matrix is
    positive then is it the same except that the x-voxel-coordinate is
    flipped first (0 to N-1 becomes N-1 to 0). Note that this means that
    our origin is in the 0,0,0 voxel (the centre of the voxel) and that we
    do not use the qform or sform information (except for the sign of the
    determinant). If you can figure out what the SPM conventions are then you
    should be able to adjust for this, although it may not be easy.

    References
    ----------
    https://github.com/nipy/nitransforms
    https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;ed6a4473.1208
    """
    raise NotImplementedError


def antsregister(template_file, li_file, lianat_file, hanat_file, outdir,
                 mask_file=None):
    """ Compute the deformation field with Ants from a T1w image to a template.
    """
    try:
        import ants
    except:
        raise ImportError("You will need to install AntsPy to execute this "
                          "function.")

    print_subtitle("Load data...")
    li = ants.image_read(li_file)
    print_result(f"li spacing: {li.spacing}")
    print_result(f"li origin: {li.origin}")
    print_result(f"li direction: {li.direction}")
    filename = os.path.join(outdir, "li.png")
    li.plot_ortho(
        flat=True, xyz_lines=False, orient_labels=False,
        title="li", filename=filename)
    lianat = ants.image_read(lianat_file)
    print_result(f"lianat spacing: {lianat.spacing}")
    print_result(f"lianat origin: {lianat.origin}")
    print_result(f"lianat direction: {lianat.direction}")
    filename = os.path.join(outdir, "lianat.png")
    lianat.plot_ortho(
        flat=True, xyz_lines=False, orient_labels=False,
        title="lianat", filename=filename)
    hanat = ants.image_read(hanat_file)
    print_result(f"hanat spacing: {hanat.spacing}")
    print_result(f"hanat origin: {hanat.origin}")
    print_result(f"hanat direction: {hanat.direction}")
    filename = os.path.join(outdir, "hanat.png")
    hanat.plot_ortho(
        flat=True, xyz_lines=False, orient_labels=False,
        title="hanat", filename=filename)
    template = ants.image_read(template_file)
    print_result(f"template spacing: {template.spacing}")
    print_result(f"template origin: {template.origin}")
    print_result(f"template direction: {template.direction}")
    filename = os.path.join(outdir, "template.png")
    template.plot_ortho(
        flat=True, xyz_lines=False, orient_labels=False,
        title="template", filename=filename)

    print_subtitle("Normalize...")
    lianat = ants.iMath_normalize(lianat)
    hanat = ants.iMath_normalize(hanat)
    template = ants.iMath_normalize(template)

    print_subtitle("Rigid: lianat -> hanat...")
    lianat2h = ants.registration(
        fixed=hanat, moving=lianat, type_of_transform="Rigid",
        outprefix=os.path.join(outdir, "lianat2h"))
    print_result(f"rigid transforms: {lianat2h['fwdtransforms']}")
    lianat2hanat = ants.apply_transforms(
        fixed=hanat, moving=lianat, transformlist=lianat2h["fwdtransforms"],
        interpolator="bSpline")
    filename = os.path.join(outdir, "lianat2hanat.nii.gz")
    lianat2hanat.to_filename(filename)
    print_result(f"lianat2h T1: {filename}")
    filename = os.path.join(outdir, "lianat2hanat.png")
    lianat2hanat.plot_ortho(
        hanat, flat=True, xyz_lines=False, orient_labels=False,
        title="lianat2hanat", filename=filename, overlay_alpha=0.5)
    li2hanat = ants.apply_transforms(
        fixed=hanat, moving=li, transformlist=lianat2h["fwdtransforms"],
        interpolator="bSpline")
    filename = os.path.join(outdir, "li2hanat.nii.gz")
    li2hanat.to_filename(filename)
    print_result(f"li2h T1: {filename}")
    filename = os.path.join(outdir, "li2hanat.png")
    li2hanat.plot_ortho(
        hanat, flat=True, xyz_lines=False, orient_labels=False,
        title="li2hanat", filename=filename, overlay_alpha=0.5)

    print_subtitle("Rigid + Affine + deformation field: hanat -> template...")
    if mask_file is None:
        h2mni = ants.registration(
            fixed=template, moving=hanat, type_of_transform="SyNRA",
            outprefix=os.path.join(outdir, "h2mni"))
    else:
        h2mni = ants.registration(
            fixed=template, moving=hanat, type_of_transform="Affine",
            outprefix=os.path.join(outdir, "_h2mni"))
        mask = ants.image_read(mask_file)
        print_result(f"mask spacing: {mask.spacing}")
        print_result(f"mask origin: {mask.origin}")
        print_result(f"mask direction: {mask.direction}")
        h2mni = ants.registration(
            fixed=template, moving=hanat, type_of_transform="SyNOnly",
            mask=mask, initial_transform=h2mni["fwdtransforms"][0],
            outprefix=os.path.join(outdir, "h2mni"))

    print_result(f"deform transforms: {h2mni['fwdtransforms']}")
    jac = ants.create_jacobian_determinant_image(
        domain_image=hanat, tx=h2mni["fwdtransforms"][0])
    jac -= 1
    hanat2mni = ants.apply_transforms(
        fixed=template, moving=hanat, transformlist=h2mni["fwdtransforms"],
        interpolator="bSpline")
    h2mnijac = ants.apply_transforms(
        fixed=template, moving=jac, transformlist=h2mni["fwdtransforms"],
        interpolator="bSpline")
    lianat2mni = ants.apply_transforms(
        fixed=template, moving=lianat, interpolator="bSpline",
        transformlist=h2mni["fwdtransforms"] + lianat2h["fwdtransforms"])
    filename = os.path.join(outdir, "hjac.nii.gz")
    jac.to_filename(filename)
    print_result(f"h jacobian: {filename}")
    filename = os.path.join(outdir, "h2mnijac.nii.gz")
    h2mnijac.to_filename(filename)
    print_result(f"h2mni jacobian: {filename}")
    filename = os.path.join(outdir, "lianat2mni.nii.gz")
    lianat2mni.to_filename(filename)
    print_result(f"li2mni T1: {filename}")
    filename = os.path.join(outdir, "hanat2mni.nii.gz")
    hanat2mni.to_filename(filename)
    print_result(f"h2mni T1: {filename}")
    filename = os.path.join(outdir, "lianat2mni.png")
    lianat2mni.plot_ortho(
        template, flat=True, xyz_lines=False, orient_labels=False,
        title="lianat2mni", filename=filename, overlay_alpha=0.5)
    filename = os.path.join(outdir, "hanat2mni.png")
    hanat2mni.plot_ortho(
        template, flat=True, xyz_lines=False, orient_labels=False,
        title="hanat2mni", filename=filename, overlay_alpha=0.5)


def apply_transforms(fixed_file, moving_file, transformlist, filename):
    """ Apply a transform list to map an image from one domain to another.

    Parameters
    ----------
    fixed_file: str
        fixed image defining domain into which the moving image is transformed.
    moving_file: str
        moving image to be mapped to fixed space.
    transformlist: list of str
        list of transforms generated by ants.registration where each transform
        is a filename.
    filename: str
        the name of the transformed image.
    """
    try:
        import ants
    except:
        raise ImportError("You will need to install AntsPy to execute this "
                          "function.")

    fixed = ants.image_read(fixed_file)
    moving = ants.image_read(moving_file)
    li2mni = ants.apply_transforms(
        fixed=fixed, moving=moving, interpolator="bSpline",
        transformlist=transformlist)
    li2mni.to_filename(filename)


def apply_translation(image_file, translation, filename):
    """ Apply a translation to an image.

    Parameters
    ----------
    image_file: str
        input image.
    translation: 3-uplet
        the translation in mm.
    filename: str
        the name of the transformed image.
    """
    im = nibabel.load(image_file)
    affine = im.affine
    affine[:3, 3] += translation
    im = nibabel.Nifti1Image(im.get_fdata(), affine)
    nibabel.save(im, filename)


def save_translation(translation, filename):
    """ Save a translation in Ants format.

    Parameters
    ----------
    translation: 3-uplet
        the translation in mm.
    filename: str
        the name of the transformed image.
    """
    try:
        import ants
    except:
        raise ImportError("You will need to install AntsPy to execute this "
                          "function.")

    ras2lps = np.array([-1, -1, 1])
    tx = ants.create_ants_transform(
        transform_type="AffineTransform", translation=(ras2lps * translation))
    ants.write_transform(tx, filename)


def ants2affine(mat_file):
    """ Map from 'in_file' image voxels to 'ref_file' voxels given `mat_file`
    omat Ants affine transformation.

    Parameters
    ------------
    mat_file: str
        filename of affine transformation file from Ants.

    Returns
    -------
    omat: array (4, 4)
        array containing the transform from voxel coordinates in image
        for 'in_file' to voxel coordinates in image for 'ref_file'.
    """
    transfo_dict = sio.loadmat(mat_file)
    lps2ras = np.diag([-1, -1, 1])
    key = list(transfo_dict.keys()).remove("fixed")[0]
    rot = transfo_dict[key][0:9].reshape((3, 3))
    trans = transfo_dict[key][9:12]
    offset = transfo_dict["fixed"]
    r_trans = (np.dot(rot, offset) - offset - trans).T * [1, 1, -1]
    omat = np.eye(4)
    omat[0: 3, 3] = r_trans
    omat[:3, :3] = np.dot(np.dot(lps2ras, rot), lps2ras)
    return omat
