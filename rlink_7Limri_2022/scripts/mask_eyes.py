import subprocess
import sys

def execute_command(command, cwd=None):
    """ Execute a command.

    Parameters
    ----------
    command: list of str
        the command to be executed.
    """
    print(" ".join(command))
    proc = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    output, error = proc.communicate()
    if proc.returncode != 0:
        raise ValueError(
            "\nCommand {0} failed:\n\n- output:\n{1}\n\n- error: "
            "{2}\n\n".format(" ".join(command),
                             output.decode("utf8"),
                             error.decode("utf8")))


def check_command(command):
    """ Check if a command is installed.

    .. note:: This function is based on which linux command.

    Parameters
    ----------
    command: str
        the name of the command to locate.
    """
    if sys.platform != "linux":
        raise ValueError("This code works only on a linux machine.")
    process = subprocess.Popen(
        ["which", command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode("utf8")
    stderr = stderr.decode("utf8")
    exitcode = process.returncode
    if exitcode != 0:
        print("Command {0}: {1}".format(command, stderr))
        raise ValueError("Impossible to locate command '{0}'.".format(command))


def mask_eyes(img_nii, cwd=None):
    ## superbet2 https://github.com/saslan-7/super-bet2/blob/master/ICV_t1.sh
    # std space
    check_command("standard_space_roi")
    command = ["standard_space_roi", img_nii, "premask", "-b"]
    execute_command(command, cwd=cwd)
    # bet brain extraction
    check_command("bet2")
    command = ["bet2", "premask", "bet2_result", "-e", "-f", "0.05"]
    execute_command(command, cwd=cwd)
    # register natif to std, keep transfo matrix
    check_command("flirt")
    command = ["flirt",
               "-ref", "/i2bm/local/fsl/data/standard/MNI152_T1_1mm_brain.nii",
               "-in", "premask",
               "-out", "flirt_result", "-omat", "flirt_result.mat"]
    execute_command(command, cwd=cwd)
    # rename
    check_command("mv")
    command = ["mv", "bet2_result_mesh.vtk", "bet2_result_mesh.off"]
    execute_command(command, cwd=cwd)
    # betsurf generate skin, brain, squelette masks in native space
    check_command("betsurf")
    command = ["betsurf", "-1", "-m", "-s", img_nii, "bet2_result_mesh.off",
               "flirt_result.mat", "final_result"]
    execute_command(command, cwd=cwd)

    ## masks
    # generate skull masked inversed 1->0 0->1
    check_command("fslmaths")
    command = ["fslmaths", "final_result_skull_mask.nii.gz", "-mul", "-1",
               "-add", "1", "skull_mask_i"]
    execute_command(command, cwd=cwd)
    # generate outskull masked inversed 1->0 0->1
    check_command("fslmaths")
    command = ["fslmaths", "final_result_outskull_mask.nii.gz", "-mul", "-1",
               "-add", "1", "outskull_mask_i"]
    execute_command(command, cwd=cwd)
    # multiply skull mask i and outskull mask i (volumes to remove)
    check_command("fslmaths")
    command = ["fslmaths", "skull_mask_i", "-mul", "outskull_mask_i", "internvol_i"]
    execute_command(command, cwd=cwd)

    ## eyes
    # raw T1w to fsl std space, keep matrix
    check_command("standard_space_roi")
    command = ["standard_space_roi", img_nii, "tmp_premask", "-b", "-d"]
    execute_command(command, cwd=cwd)
    # register average mask eye (std space) on T1w(native space)
    check_command("flirt")
    command = ["flirt", "-ref", img_nii, "-in",
               "/i2bm/local/fsl/data/standard/MNI152_T1_2mm_eye_mask",
               "-applyxfm", "-init", "tmp_premask_tmp_to_std_inv.mat",
               "-datatype", "float", "-out", "tmp_eye_mask_nativspace"]
    execute_command(command, cwd=cwd)
    # binarise eye mask
    check_command("fslmaths")
    command = ["fslmaths", "tmp_eye_mask_nativspace.nii.gz", "-bin",
               "eye_mask_nativspace1"]
    execute_command(command, cwd=cwd)

    ## apply masks
    # remove intern volumes
    check_command("fslmaths")
    command = ["fslmaths", "eye_mask_nativspace1", "-mul", "internvol_i",
               "eye_mask_nativspace2"]
    execute_command(command, cwd=cwd)
    # remove extern volumes
    check_command("fslmaths")
    command = ["fslmaths", "eye_mask_nativspace2", "-mul",
               "final_result_outskin_mask.nii.gz",
               "eye_mask_nativspace_final"]
    execute_command(command, cwd=cwd)
    return 0


img_nii = "/neurospin/psy_sbox/temp_julie/Fawzi/example/sub-20181116/Anatomy3T/t1_weighted_sagittal_1_0iso.nii"
woking_directory = "/neurospin/psy_sbox/temp_julie/Fawzi/mask_eyes"
mask_eyes(img_nii, cwd=woking_directory)
