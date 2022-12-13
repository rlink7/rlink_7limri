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
import glob
import argparse


def apply_transform(moving_file, matrix, output):
    applymat = spmu.ApplyTransform()
    applymat.inputs.in_file = moving_file
    applymat.inputs.mat = matrix  # .mat
    applymat.inputs.out_file = output
    applymat.run()
    return 0


def coreg(target_file, moving_file):
    coreg = spm.Coregister()
    coreg.inputs.target = target_file
    coreg.inputs.source = moving_file
    coreg.inputs.use_mcr = True
    coreg.run()


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


# initialization
matlab_cmd = "/i2bm/local/cat12-standalone/run_spm12.sh"\
            " /i2bm/local/cat12-standalone/mcr/v93/ script"
spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
print("spm mcr version", spm.SPMCommand().version)
print(nipype.__version__)
# standalone cat12vbm matlab config
executable_cat12 = '/i2bm/local/cat12-standalone/standalone/cat_standalone.sh'
standalone_cat12 = '/i2bm/local/cat12-standalone'
mcr_matlab = '/i2bm/local/cat12-standalone/mcr/v93'
# templates
tpm_file = '/i2bm/local/cat12-standalone/spm12_mcr/home/gaser/gaser/'\
           'spm/spm12/tpm/TPM.nii'
darteltpm_file = '/i2bm/local/cat12-standalone/spm12_mcr/'\
                 'home/gaser/gaser/spm/spm12/toolbox/cat12/tem'\
                 'plates_volumes/Template_1_IXI555_MNI152.nii'

# paths
# example2
# root = "/neurospin/psy_sbox/temp_julie/Fawzi/example/written/example2"
# target_file = os.path.join(root,
#               "sub-01004_ses-M03_acq-3DT1_rec-yBCyGC_run-1_T1w.nii")
# moving_file = os.path.join(root,
#               "sub-01004_ses-M03Li_acq-3DT1Li_rec-nBCyGC_run-1_T1w.nii")
# output = os.path.join(root, "anatli_anatspace.nii")
# matrix = os.path.join(root, "anatLi_to_anat.mat")
# moving_file2 = os.path.join(root, "rsub-01004_ses-M03Li_acq-3DT1Li_rec-nBCyGC_run-1_T1w.nii")
# defor = os.path.join(root, "y_sub-01004_ses-M03_acq-3DT1_rec-yBCyGC_run-1_T1w.nii")

# # example3
# root = "/neurospin/psy_sbox/temp_julie/Fawzi/example/written/example3"
# target_file = os.path.join(root,
#               "sub-20181116_t1_weighted_sagittal_1_0iso.nii")
# moving_file = os.path.join(root,
#               "sub-20181116_t1_mpr_tra_iso1_0mm.nii")
# output = os.path.join(root, "anatli_anatspace.nii")
# matrix = os.path.join(root, "anatLi_to_anat.mat")
# moving_file2 = os.path.join(root, "rsub-20181116_t1_mpr_tra_iso1_0mm.nii")
# defor = os.path.join(root, "y_sub-20181116_t1_weighted_sagittal_1_0iso.nii")

# # example4
root = "/neurospin/psy_sbox/temp_julie/Fawzi/example/written/example4"
target_file = os.path.join(root,
              "sub-20181221_ses-1_acq-3T_T1w.nii")
moving_file = os.path.join(root,
              "sub-20181221_ses-1_acq-7T_T1w.nii")
output = os.path.join(root, "anatli_anatspace.nii")
matrix = os.path.join(root, "anatLi_to_anat.mat")
moving_file2 = os.path.join(root, "rsub-20181221_ses-1_acq-7T_T1w.nii")
defor = os.path.join(root, "y_sub-20181221_ses-1_acq-3T_T1w.nii")

apply_transform(moving_file, matrix, output)
coreg(target_file, moving_file)
apply_defor(defor, moving_file2)