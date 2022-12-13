import nibabel
import os
import scipy.io
import nipype
import nipype.interfaces.spm as spm
import nipype.interfaces.spm.utils as spmu
import numpy as np
from skimage.restoration import denoise_nl_means, estimate_sigma
from nilearn import plotting
from nipype.interfaces import cat12
import subprocess
import rlink_7Limri_2022
import argparse


def denoising_nlm(image, output):
    image = nibabel.load(image)
    arr_image = image.get_fdata()
    sigma_est = np.mean(estimate_sigma(arr_image))
    arr_denoised_image = denoise_nl_means(arr_image, h=1.15 * sigma_est,
                                          patch_size=9,
                                          patch_distance=5)
    denoised_image = nibabel.Nifti1Image(arr_denoised_image, image.affine)
    nibabel.save(denoised_image, output)
    return 0


def denoising_cat12(registered_Li):
    # non-local means (SANLM)
    # denoising filter (Manjón et al., 2010)
    c = cat12.CAT12SANLMDenoising()
    c.inputs.in_files = registered_Li
    c.run()
    return 0


def coregistration(target_file, moving_file, coreg_path):
    coreg = spmu.CalcCoregAffine()
    coreg.inputs.target = target_file
    coreg.inputs.moving = moving_file
    coreg.inputs.mat = coreg_path
    coreg.inputs.use_mcr = True
    coreg.run()
    return 0


def apply_defor(defor, moving_file):
    # apply H coil to MNI
    norm12 = spm.Normalize12()
    norm12.inputs.apply_to_files = moving_file
    norm12.inputs.deformation_file = defor
    # norm12.inputs.write_bounding_box = [[-78, -112, -70], [78, 76, 85]]
    norm12.inputs.write_bounding_box = [[-90, -126, -72], [90, 90, 108]]
    norm12.inputs.write_voxel_sizes = [1.5, 1.5, 1.5]
    norm12.inputs.write_interp = 4
    norm12.inputs.jobtype = 'write'
    norm12.inputs.use_mcr = True
    # norm12.inputs.out_prefix = "MNI_"
    norm12.run()
    return 0


def apply_transform(input, matrix, output):
    applymat = spmu.ApplyTransform()
    applymat.inputs.in_file = input
    applymat.inputs.mat = matrix  # .mat
    applymat.inputs.out_file = output
    applymat.run()
    return 0


def resclice2(input, target):
    r2ref = spmu.ResliceToReference()
    r2ref.inputs.in_files = input
    r2ref.inputs.target = target
    r2ref.run()
    return 0


def resclice(input, target, output):
    r2ref = spmu.Reslice()
    r2ref.inputs.in_file = input
    r2ref.inputs.space_defining = target
    r2ref.inputs.out_file = output
    r2ref.run()
    return 0


def scale(imfile, scaledfile, scale):
    """ Scale the MRI image.
    .. note:: This function is based on FSL.
    Parameters
    ----------
    imfile: str
        the input image.
    scaledfile: str
        the path to the scaled input image.
    scale: int
        the scale factor in all directions.
    Returns
    -------
    scaledfile, trffile: str
        the generated files.
    """
    trffile = scaledfile.split(".")[0] + ".txt"
    cmd = ["flirt", "-in", imfile, "-ref", imfile, "-out",
           scaledfile, "-applyisoxfm", str(scale), "-omat", trffile]
    subprocess.check_call(cmd)
    return scaledfile, trffile


def dicom_to_nii(liste_filename):
    di = spmu.DicomImport()
    di.inputs.in_files = liste_filename
    di.run()
    return 0


def write_matlabbatch(template, nii_file, tpm_file, darteltpm_file, outfile):
    """ Complete matlab batch from template.

    Parameters
    ----------
    template: str
        path to template batch to be completed.
    nii_files: list
        the Nifti image to be processed.
    tpm_file: str
        path to the SPM TPM file.
    darteltpm_file: str
        path to the CAT12 tempalte file.
    outfile: str
        path to the generated matlab batch file that can be used to launch
        CAT12 VBM preprocessing.
    """
    nii_file_str = ""
    for i in nii_file:
        nii_file_str += "'{0}' \n".format(i)
    with open(template, "r") as of:
        stream = of.read()
    stream = stream.format(anat_file=nii_file_str, tpm_file=tpm_file,
                           darteltpm_file=darteltpm_file)
    with open(outfile, "w") as of:
        of.write(stream)
    return 0


def plot_anat_Li(li_mni_denoised, anat_MNI, threshold, figname):
    li_img = nibabel.load(li_mni_denoised)
    bg_img = nibabel.load(anat_MNI)
    arr = li_img.get_fdata()
    arr[arr < threshold] = 0
    li_img = nibabel.Nifti1Image(arr, li_img.affine)
    display = plotting.plot_stat_map(li_img, bg_img, cut_coords=(-35, 54, -44),
                                     cmap=plotting.cm.black_blue)
    display.savefig(figname)
    # plotting.show()
    display.close()
    return 0


def find_threshold(mask_MNI, li_mni_denoised):
    li_img = nibabel.load(li_mni_denoised)
    mask_img = nibabel.load(mask_MNI)
    arr = li_img.get_fdata()
    arr[mask_img.get_fdata() != 0] = 0
    threshold = np.max(arr)
    print(threshold)
    return threshold


def pipeline_lithium(target_anatLi, target_anat, moving_file_Li,
                     transfo_folder,
                     executable_cat12, standalone_cat12,
                     mcr_matlab, matlabbatch, tpm_file, darteltpm_file,
                     threshold=None, bfc=True):
    print("pipeline launch")
    if bfc:
        ''' biasfield correction'''
        # anat
        bfcfile_anat = target_anat.replace(".nii", "_bfc.nii")
        bfcfile_anat = os.path.join(transfo_folder,
                                    os.path.basename(bfcfile_anat))
        if not os.path.exists(bfcfile_anat):
            bfcfile_anat, _ = biasfield(target_anat, bfcfile_anat,
                                        maskfile=None,
                                        nb_iterations=50,
                                        convergence_threshold=0.001,
                                        bspline_grid=(1, 1, 1),
                                        shrink_factor=1,
                                        bspline_order=3,
                                        histogram_sharpening=(0.15, 0.01, 200),
                                        check_pkg_version=False)
            target_anat = bfcfile_anat
        print("anat biasfield correction ok")
        # anat Li
        bfcfile_anatLi = target_anat.replace(".nii", "_bfc.nii")
        bfcfile_anatLi = os.path.join(transfo_folder,
                                      os.path.basename(bfcfile_anatLi))
        if not os.path.exists(bfcfile_anatLi):
            bfcfile_anatLi, _ = biasfield(target_anat, bfcfile_anatLi,
                                          maskfile=None,
                                          nb_iterations=50,
                                          convergence_threshold=0.001,
                                          bspline_grid=(1, 1, 1),
                                          shrink_factor=1,
                                          bspline_order=3,
                                          histogram_sharpening=(0.15,
                                                                0.01,
                                                                200),
                                          check_pkg_version=False)
            target_anatLi = bfcfile_anatLi
        print("anatLi biasfield correction ok")

    ''' transformations estimations'''
    # # Li to anat Li # does not work
    # coreg_path_Li = os.path.join(transfo_folder, "Li_to_Lianat.mat")
    # if not os.path.exists(coreg_path_Li):
    #     coregistration(target_anatLi, moving_file_Li, coreg_path_Li)
    # print("coreg Li Lianat ok")
    # anat Li to anat
    coreg_path_anat = os.path.join(transfo_folder, "anatLi_to_anat.mat")
    if not os.path.exists(coreg_path_anat):
        coregistration(target_anat, target_anatLi, coreg_path_anat)
    print("coreg anatLi anat ok")

    # Non Linear registration
    deformfiel_nl = target_anat.split(os.sep)
    deformfiel_nl.insert(-1, "mri")
    deformfiel_nl[-1] = "y_{0}".format(deformfiel_nl[-1])
    deformfiel_nl = os.sep.join(deformfiel_nl)

    if not os.path.exists(deformfiel_nl):
        # write matlabbatch
        print("matlabbatch")
        batch_name = os.path.join(transfo_folder, "matlabbatch.m")
        write_matlabbatch(matlabbatch, [target_anat], tpm_file, darteltpm_file,
                          batch_name)
        # anat to MNI
        print("launch car12vbm")
        subprocess.check_call([executable_cat12,
                               "-s", standalone_cat12,
                               "-m", mcr_matlab,
                               "-b", batch_name])
    print("cat12vbm ok")

    '''Transformations applications'''
    # trans_Lianat = scipy.io.loadmat(os.path.join(transfo_folder,
    #                                 "inverse_Li_to_Lianat.mat"))
    # trans_Lianat_mat = trans_Lianat['M']
    trans_anat = scipy.io.loadmat(os.path.join(transfo_folder,
                                  "inverse_anatLi_to_anat.mat"))
    trans_anat_mat = trans_anat['M']
    # Linear registrations combinaison
    # mixte_mat = np.dot(trans_anat_mat, trans_Lianat_mat)
    mixte_mat = trans_anat_mat
    # Modify Li affine, linear registration
    li_modified_affine = os.path.join(transfo_folder, "li_modified_affine.nii")
    if not os.path.exists(li_modified_affine):
        img = nibabel.load(moving_file_Li)
        new_affine = np.dot(mixte_mat, img.affine)
        normalized = nibabel.Nifti1Image(img.get_fdata(), new_affine)
        nibabel.save(normalized, li_modified_affine)
    print("Modified Li affine ok")
    # Modify Li affine, linear registration
    anatli_modified_affine = os.path.join(transfo_folder,
                                          "anatli_modified_affine.nii")
    if not os.path.exists(anatli_modified_affine):
        img = nibabel.load(target_anatLi)
        new_affine = np.dot(mixte_mat, img.affine)
        normalized = nibabel.Nifti1Image(img.get_fdata(), new_affine)
        nibabel.save(normalized, anatli_modified_affine)
    print("Modified anatLi affine ok")
    # Apply NL deformation on Li
    Li_MNI = os.path.join(transfo_folder, "wli_modified_affine.nii")
    if not os.path.exists(Li_MNI):
        apply_defor(deformfiel_nl, li_modified_affine)
    print("NL deformation applied on 7Li")
    # Apply NL deformation on anat
    anat_MNI_auto = os.path.join(os.path.dirname(target_anat),
                                 "w{0}".format(os.path.basename(target_anat)))
    anat_MNI = os.path.join(transfo_folder,
                            "w{0}".format(os.path.basename(target_anat)))
    if not os.path.exists(anat_MNI):
        apply_defor(deformfiel_nl, target_anat)
        if anat_MNI_auto != anat_MNI:
            subprocess.check_call(["mv", anat_MNI_auto, anat_MNI])
    print("NL deformation applied on anat 3T")
    # Apply NL deformation on Lianat
    Lianat_MNI_auto = os.path.join(os.path.dirname(anatli_modified_affine),
                                   "w{0}".format(os.path.basename
                                                 (anatli_modified_affine)))
    if not os.path.exists(Lianat_MNI_auto):
        apply_defor(deformfiel_nl, anatli_modified_affine)
    print("NL deformation applied on anat Li")
    # Apply denoising on Li MNI
    name_denoised_Li = os.path.join(transfo_folder,
                                    "sanlm_wli_modified_affine.nii")
    if not os.path.exists(name_denoised_Li):
        denoising_cat12(Li_MNI)
    print("denoising finished : ", name_denoised_Li)
    # MNI space shape checks
    li_img = nibabel.load(name_denoised_Li)
    bg_img = nibabel.load(anat_MNI)
    print("Li MNI shape : ", li_img.get_fdata().shape)
    print("Anat MNI shape : ", bg_img.get_fdata().shape)
    assert li_img.get_fdata().shape ==\
           (121, 145, 121), li_img.get_fdata().shape
    assert li_img.get_fdata().shape ==\
           (121, 145, 121), li_img.get_fdata().shape
    # Threshold
    if not threshold:
        arr = li_img.get_fdata()
        threshold = np.max(arr)*0.25

    ''' Plot results'''
    figname_no_denoi = os.path.join(transfo_folder, "figure_li_anat_MNI")
    fignamedenoi = os.path.join(transfo_folder, "figure_denoised_li_anat_MNI")
    plot_anat_Li(name_denoised_Li, anat_MNI, threshold, fignamedenoi)
    plot_anat_Li(Li_MNI, anat_MNI, threshold, figname_no_denoi)

    return 0


def biasfield(imfile, bfcfile, maskfile=None, nb_iterations=50,
              convergence_threshold=0.001, bspline_grid=(1, 1, 1),
              shrink_factor=1, bspline_order=3,
              histogram_sharpening=(0.15, 0.01, 200), check_pkg_version=False):
    """ Perform MRI bias field correction using N4 algorithm.
    .. note:: This function is based on ANTS.
    Parameters
    ----------
    imfile: str
        the input image.
    bfcfile: str
        the bias fieled corrected file.
    maskfile: str, default None
        the brain mask image.
    nb_iterations: int, default 50
        Maximum number of iterations at each level of resolution. Larger
        values will increase execution time, but may lead to better results.
    convergence_threshold: float, default 0.001
        Stopping criterion for the iterative bias estimation. Larger values
        will lead to smaller execution time.
    bspline_grid: int, default (1, 1, 1)
        Resolution of the initial bspline grid defined as a sequence of three
        numbers. The actual resolution will be defined by adding the bspline
        order (default is 3) to the resolution in each dimension specified
        here. For example, 1,1,1 will result in a 4x4x4 grid of control points.
        This parameter may need to be adjusted based on your input image.
        In the multi-resolution N4 framework, the resolution of the bspline
        grid at subsequent iterations will be doubled. The number of
        resolutions is implicitly defined by Number of iterations parameter
        (the size of this list is the number of resolutions).
    shrink_factor: int, default 1
        Defines how much the image should be upsampled before estimating the
        inhomogeneity field. Increase if you want to reduce the execution
        time. 1 corresponds to the original resolution. Larger values will
        significantly reduce the computation time.
    bspline_order: int, default 3
        Order of B-spline used in the approximation. Larger values will lead
        to longer execution times, may result in overfitting and poor result.
    histogram_sharpening: 3-uplate, default (0.15, 0.01, 200)
        A vector of up to three values. Non-zero values correspond to Bias
        Field Full Width at Half Maximum, Wiener filter noise, and Number of
        histogram bins.
    check_pkg_version: bool, default False
        optionally check the package version using dpkg.
    Returns
    -------
    bfcfile, bffile: str
        the generatedd files.
    """
    ndim = 3
    bspline_grid = [str(e) for e in bspline_grid]
    histogram_sharpening = [str(e) for e in histogram_sharpening]
    bffile = bfcfile.split(".")[0] + "_field.nii.gz"
    cmd = [
        "N4BiasFieldCorrection",
        "-d", str(ndim),
        "-i", imfile,
        "-s", str(shrink_factor),
        "-b", "[{0}, {1}]".format("x".join(bspline_grid), bspline_order),
        "-c", "[{0}, {1}]".format(
            "x".join([str(nb_iterations)] * 4), convergence_threshold),
        "-t", "[{0}]".format(", ".join(histogram_sharpening)),
        "-o", "[{0}, {1}]".format(bfcfile, bffile),
        "-v"]
    if maskfile is not None:
        cmd += ["-x", maskfile]
    subprocess.check_call(cmd)
    return bfcfile, bffile


###############################################################################
def main():

    # PARSER
    parser = argparse.ArgumentParser("Launch 7Li preprocessing")
    # required
    parser.add_argument('--Li',
                        help='path to Li image',
                        required=True,
                        type=str)
    parser.add_argument('--Lianat',
                        help='path of anat image'
                             ' acquired with Li coil',
                        required=True,
                        type=str)
    parser.add_argument('--anat',
                        help='path of anat image'
                             ' acquired with H coil',
                        required=True,
                        type=str)
    parser.add_argument('--output',
                        help='path the preprocessing outputs',
                        required=True,
                        type=str)
   
    # optionnal
    parser.add_argument('--executable_cat12',
                        help='Path to the executable_cat12.\n'
                             'example : [...]/cat12-standalone/'
                             'standalone/cat_standalone.sh',
                        default='/i2bm/local/cat12-standalone/'
                                'standalone/cat_standalone.sh',
                        type=str)
    parser.add_argument('--standalone_cat12',
                        help='cat12 standalone folder.\n'
                             'example : [...]/cat12-standalone',
                        default='/i2bm/local/cat12-standalone',
                        type=str,)
    parser.add_argument('--mcr_matlab',
                        help='Path to the mcr_matlab.\n'
                             'example : [...]/cat12-standalone/mcr/v93',
                        default='/i2bm/local/cat12-standalone/mcr/v93',
                        type=str)
    parser.add_argument('--tpm_file',
                        help='Path to the tpm_file cat12vbm file.\n'
                             'example : cat12-standalone/spm12_mcr/'
                             'home/gaser/gaser/spm/spm12/tpm/TPM.nii',
                        default='/i2bm/local/cat12-standalone/spm12_mcr/'
                                'home/gaser/gaser/spm/spm12/tpm/TPM.nii',
                        type=str)
    parser.add_argument('--darteltpm_file',
                        help='Path to the dartel tpm file.\n'
                             'example : [...]/cat12-standalone/spm12_mcr/'
                             'home/gaser/gaser/spm/spm12/toolbox/cat12/'
                             'templates_volumes/Template_1_IXI555_MNI152.nii',
                        default='/i2bm/local/cat12-standalone/spm12_mcr/'
                                'home/gaser/gaser/spm/spm12/toolbox/cat12/tem'
                                'plates_volumes/Template_1_IXI555_MNI152.nii',
                        type=str)
    parser.add_argument('--threshold',
                        help="threshold",
                        default=None,
                        type=int)
    options = parser.parse_args()
    # required
    Li = options.Li
    Lianat = options.Lianat
    anat = options.anat
    output = options.output
    # initialization
    matlab_cmd = "/i2bm/local/cat12-standalone/run_spm12.sh"\
                 " /i2bm/local/cat12-standalone/mcr/v93/ script"
    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
    print("spm mcr version", spm.SPMCommand().version)
    print(nipype.__version__)
    # standalone cat12vbm matlab config
    executable_cat12 = options.executable_cat12
    standalone_cat12 = options.standalone_cat12
    mcr_matlab = options.mcr_matlab
    # templates
    tpm_file = options.tpm_file
    darteltpm_file = options.darteltpm_file
    # resources dir
    resource_dir = os.path.join(
            os.path.dirname(rlink_7Limri_2022.__file__), "resources")
    matlabbatch = os.path.join(resource_dir, "cat12vbm_matlabbatch.m")

    # threshold 20 or 200
    threshold = options.threshold

    # computation
    sub = os.path.basename(Li).split("_")[0]
    print("launch {0}".format(sub))
    transfo_folder = os.path.join(output, "preproc-lithium_{0}".format(sub))
    subprocess.check_call(["mkdir", "-p", transfo_folder])
    pipeline_lithium(Lianat, anat, Li,
                     transfo_folder,
                     executable_cat12, standalone_cat12,
                     mcr_matlab, matlabbatch, tpm_file,
                     darteltpm_file, threshold, bfc=True)


if __name__ == "__main__":
    main()
