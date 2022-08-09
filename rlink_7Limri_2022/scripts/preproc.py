import nibabel
import os
import scipy.io
import nipype.interfaces.spm as spm
import nipype.interfaces.spm.utils as spmu
import numpy as np
from skimage.restoration import denoise_nl_means, estimate_sigma
from nilearn import plotting
from nipype.interfaces import cat12
import subprocess
import rlink_7Limri_2022
import glob


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
    coreg.run()
    return 0


def apply_defor(defor, moving_file):
    # apply H coil to MNI
    norm12 = spm.Normalize12()
    norm12.inputs.apply_to_files = moving_file
    norm12.inputs.deformation_file = defor
    # norm12.inputs.write_bounding_box = [[-78, -112, -70], [78, 76, 85]]
    # norm12.inputs.write_bounding_box = [[-42, -68, -35], [78, 76, 85]]
    # norm12.inputs.write_voxel_sizes = [1, 1, 1]
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
                           cmap=plotting.cm.black_blue)  # , output_file=".png"
    display.savefig(figname)
    plotting.show()
    display.close()
    
    # plotting.show()
    return 0


def find_threshold(mask_MNI, li_mni_denoised):
    li_img = nibabel.load(li_mni_denoised)
    mask_img = nibabel.load(mask_MNI)
    arr = li_img.get_fdata()
    arr[mask_img.get_fdata() != 0] = 0
    threshold = np.max(arr)
    print(threshold)
    return threshold


def bids_naming(path, liregex, anatliregex, anatregex):
    target_anatLi = os.path.join(path, anatliregex)
    list_anat_Li = glob.glob(target_anatLi)
    assert len(list_anat_Li) != 0, target_anatLi
    target_anat = os.path.join(path, anatregex)
    list_anat = glob.glob(target_anat)
    assert len(list_anat) != 0, target_anat
    moving_file_Li = os.path.join(path, liregex)
    list_Li = glob.glob(moving_file_Li)
    assert len(list_Li) != 0, moving_file_Li
    list_sub = [i.split(os.sep)[-4] for i in list_Li]
    dico = {}
    for c, i in enumerate(list_sub):
        dico[i] = []
        dico[i].append(list_Li[c])
    for i in list_anat_Li:
        sub = i.split(os.sep)[-4]
        if sub in dico:
            dico[sub].append(i)
    for i in list_anat:
        sub = i.split(os.sep)[-4]
        if sub in dico:
            dico[sub].append(i)
    for i in dico:
        if len(dico[i]) < 3:
            del dico[i]
    print("number of subjects to preprocess : ", len(dico))
        
    return dico



def pipeline_lithium(target_anatLi, target_anat, moving_file_Li,
                     transfo_folder,
                     executable_cat12, standalone_cat12,
                     mcr_matlab, matlabbatch, tpm_file, darteltpm_file,
                     threshold):
    print("pipeline launch")
    # Li to anat Li
    coreg_path_Li = os.path.join(transfo_folder, "Li_to_Lianat.mat")
    if not os.path.exists(coreg_path_Li):
        # print("coreg Li Lianat")
        coregistration(target_anatLi, moving_file_Li, coreg_path_Li)
    print("coreg Li Lianat ok")
    # anat Li to anat
    coreg_path_anat = os.path.join(transfo_folder, "anatLi_to_anat.mat")
    if not os.path.exists(coreg_path_anat):
        # print("coreg anatLi anat")
        coregistration(target_anat, target_anatLi, coreg_path_anat)
    print("coreg anatLi anat ok")

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
    ############################################################################
    # create combine transformation
    trans_Lianat = scipy.io.loadmat(os.path.join(transfo_folder,
                                    "inverse_Li_to_Lianat.mat"))
    trans_Lianat_mat = trans_Lianat['M']

    trans_anat = scipy.io.loadmat(os.path.join(transfo_folder,
                                  "inverse_anatLi_to_anat.mat"))
    trans_anat_mat = trans_anat['M']

    # mixte_mat = np.dot(trans_anat_mat, trans_Lianat_mat) # to determine
    mixte_mat = trans_anat_mat
    #### test no combine
    img = nibabel.load(moving_file_Li)
    new_affine = np.dot(mixte_mat, img.affine)
    normalized = nibabel.Nifti1Image(img.get_fdata(), new_affine)

    # img = nibabel.load(deformfiel_nl)
    # new_affine = np.dot(mixte_mat, img.affine)
    # normalized = nibabel.Nifti1Image(img.get_fdata(), new_affine)
    nibabel.save(normalized, os.path.join(transfo_folder, "li_linearapply.nii"))
    # defor = os.path.join(transfo_folder, "y_simple.nii")
    defor = deformfiel_nl
  
    
    new_name = os.path.join(transfo_folder, "li_linearapply.nii")
    apply_defor(defor, new_name) 
    # check point
    Li_MNI = new_name.split(os.sep)
    Li_MNI[-1] = "w{0}".format(Li_MNI[-1]) 
    Li_MNI = os.sep.join(Li_MNI)
    img_Li = nibabel.load(Li_MNI)
    # assert img_Li.shape == (121, 145, 121), img_Li.shape  # (157, 189, 156)

    apply_defor(defor, target_anat)

    denoising_cat12(Li_MNI)
    print("sanlm_wli_linearapply.nii")
    name_denoised_Li = os.path.join(transfo_folder, "sanlm_wli_linearapply.nii")

    # # find threshold
    # threshold = find_threshold(mask_path_MNI, Li_MNI_denoised)
        
    # plot results
    anat_MNI = target_anat.split(os.sep)
    anat_MNI[-1] = "w{0}".format(anat_MNI[-1])
    anat_MNI = os.sep.join(anat_MNI)
    fignamenodenoi = os.path.join(transfo_folder, "figureliandanat")
    fignamedenoi = os.path.join(transfo_folder, "figure_denoised_liandanat")
    plot_anat_Li(name_denoised_Li, anat_MNI, threshold, fignamenodenoi)
    plot_anat_Li(Li_MNI, anat_MNI, threshold, fignamedenoi)


################################################################################
# LAUNCH
# initialization
matlab_cmd = '/i2bm/local/cat12-standalone/run_spm12.sh /i2bm/local/cat12-standalone/mcr/v93/ script'
# matlab_cmd = '/i2bm/local/spm12-standalone/run_spm12.sh /i2bm/local/spm12-standalone/mcr/v713/ script'
# matlab_cmd = '/neurospin/local/bin/matlab-R2019a'
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True) 

# standalone cat12vbm matlab config
executable_cat12 = "/i2bm/local/cat12-standalone/standalone/cat_standalone.sh"
standalone_cat12 = "/i2bm/local/cat12-standalone"
mcr_matlab = "/i2bm/local/cat12-standalone/mcr/v93"
# templates
tpm_file = "/i2bm/local/cat12-standalone/spm12_mcr/home/gaser/gaser/spm/spm12/tpm/TPM.nii"
darteltpm_file = "/i2bm/local/cat12-standalone/spm12_mcr/home/gaser/gaser/spm/spm12/toolbox/"\
                "cat12/templates_volumes/Template_1_IXI555_MNI152.nii"

# resources dir
resource_dir = os.path.join(
        os.path.dirname(rlink_7Limri_2022.__file__), "resources")
matlabbatch = os.path.join(resource_dir, "cat12vbm_matlabbatch.m")

# options
# method : one or list
method = "list"
# launch: local or cluster
launch = "local"
# threshold 20 or 200
threshold = 200

if method == "list":
    path = "/neurospin/ciclops/projects/BIPLi7/BIDS/derivatives/Limri_preproc"
    liregex = "sub-*/ses-*/lithium/sub-*_ses-*_acq-trufi_limri.nii"
    anatliregex = "sub-*/ses-*/anat/sub-*_ses-*_acq-7T_T1w.nii"
    anatregex = "sub-*/ses-*/anat/sub-*_ses-*_acq-3T_T1w.nii"
    dico = bids_naming(path, liregex, anatliregex, anatregex)
    
    if launch == "local":
        for i in dico:
            if i == "sub-20181207":
                print("launch {0}".format(i))
                transfo_folder = os.path.join(os.path.dirname(dico[i][0]), "transfo")
                subprocess.check_call(["mkdir", "-p", transfo_folder])
                target_anatLi = dico[i][0]
                target_anat = dico[i][2]
                moving_file_Li = dico[i][1]
                name_denoised_Li = os.path.join(transfo_folder, "sanlm_wli_linearapply.nii")
                if not os.path.exists(name_denoised_Li):  # do not redoo if it s already done
                    pipeline_lithium(target_anatLi, target_anat, moving_file_Li,
                                    transfo_folder,
                                    executable_cat12, standalone_cat12,
                                    mcr_matlab, matlabbatch, tpm_file, darteltpm_file,
                                    threshold)
    elif launch == "cluster":
        pass

elif method == "one":
    # inputs
    # target_anatLi = "/volatile/7Li/data_test/08_09_21/data_rlink/anat_3TLi/01004RL20210430M033DT1Li_noPN_DIS2D_S007.nii"
    # target_anat = "/volatile/7Li/data_test/08_09_21/data_rlink/anat_3T/01004RL20210430M033DT1_noPN_DIS2D_S009.nii"
    # moving_file_Li = "/volatile/7Li/data_test/08_09_21/data_rlink/Li/01004RL20210430M03trufi_S005.nii"
    # transfo_folder = "/volatile/7Li/data_test/08_09_21/transfo_rlink"

    target_anatLi = "/neurospin/psy_sbox/temp_julie/Fawzi/example/sub-20181116/Anatomy7T/t1_mpr_tra_iso1_0mm.nii"
    target_anat = "/neurospin/psy_sbox/temp_julie/Fawzi/example/sub-20181116/Anatomy3T/t1_weighted_sagittal_1_0iso.nii"
    moving_file_Li = "/neurospin/psy_sbox/temp_julie/Fawzi/example/sub-20181116/Trufi/01-Raw/TRUFI_1000_1.nii"
    transfo_folder = "/neurospin/psy_sbox/temp_julie/Fawzi/example/transfo"

    pipeline_lithium(target_anatLi, target_anat, moving_file_Li,
                    transfo_folder,
                    executable_cat12, standalone_cat12,
                    mcr_matlab, matlabbatch, tpm_file, darteltpm_file,
                    threshold)

elif method == "plot":
    transfo_folder = "/neurospin/psy_sbox/temp_julie/Fawzi/example/transfo"
    target_anat = "/neurospin/psy_sbox/temp_julie/Fawzi/example/sub-20181116/Anatomy3T/t1_weighted_sagittal_1_0iso.nii"
    li_mni_denoised = os.path.join(transfo_folder, "sanlm_wli_linearapply.nii")
    anat_MNI = target_anat.split(os.sep)
    anat_MNI[-1] = "w{0}".format(anat_MNI[-1])
    anat_MNI = os.sep.join(anat_MNI)
    plot_anat_Li(li_mni_denoised, anat_MNI, threshold)
