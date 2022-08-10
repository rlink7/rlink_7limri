|Python27|_ |Python35|_ |PyPi|_

.. |Python27| image:: https://img.shields.io/badge/python-2.7-blue.svg
.. _Python27: https://badge.fury.io/py/pycaravel

.. |Python35| image:: https://img.shields.io/badge/python-3.5-blue.svg
.. _Python35: https://badge.fury.io/py/pycaravel

.. |PyPi| image:: https://badge.fury.io/py/pycaravel.svg
.. _PyPi: https://badge.fury.io/py/pycaravel


rlink_7Limri_2022
=================

This module has been created to do some preprocessings on lithium mri data.


Important links
===============


Dependencies
============

The required dependencies to use the software are:



Install
=======

Check the official documentation dedicated page.


How to use
==========

Neurospin specifics:
Activate the env :
.. code::

    source /neurospin/tmp/jv261711/env/li/bin/activate
Copy/Paste the Rawdata to /derivatives/preproc_limri/ (this will be your input folder)

Options
=======
# method : "one" or "list"
You can launch the preprocessing on one subject, with the "one" option.
You can launch the preprocessing on several subjects, with the "list" option.

# launch: "local" or "cluster"
You can launch the preprocessing on your local computer.
You can launch the preprocessing on the alambic cluster. (neurospin specifics)

# file.txt : 1 or 0
If you know exactly which subjetcs you want to preprocess, you can launch the preprocessing from txt files, set the option to 1.
If you want to use regex or to use the one method, set the option to 0.

# threshold
You can define a specific threshold for the lithium imaging.

Inputs
======

If method = "one" :
* target_anatLi : path to the anatomic mri acquired with the Li coil. (T1w.nii)
* target_anat : path to the anatomic mri acquired with the H coil. (T1w.nii)
* moving_file_Li : path to the Li mri mri acquired with the Li coil. (lithium.nii)
* transfo_folder = path to the output folder.

If method = "list" :
    If file_txt = 0:
        * path = path to the input folder. (common path before the regex)
        * liregex = regex of Li mri mri acquired with the Li coil. (lithium.nii)
        * anatliregex = regex of the anatomic mri acquired with the Li coil. (T1w.nii)
        * anatregex = regex of the anatomic mri acquired with the H coil. (T1w.nii)
    If file_txt = 1:
        * path_Li = Path to the Li.txt file : content : one row per subject.
        * path_anat_Li = Path to the anatLi.txt file : content : one row per subject.
        * path_anat = Path to the anatH.txt file : content : one row per subject.

Preprocessing Pipeline Description
==================================

Estimate and apply transformation
* linear transformation : anat Li coil => anat H coil.
* non linear transformation : anat H coil => MNI Template.
* Apply combinaison of transfomations.
* Apply denoising on Li.
* Save Li and anat H in the MNI space.
* Save plots.

Outputs
=======
* anatLi_to_anat.mat => linear transfo matrix from anat mri Li coil to anat mri Hcoil
* inverse_anatLi_to_anat.mat		 
* Li_to_Lianat.mat => linear transfo matrix from Li mri Li coil to anat mri Licoil (head motions)
* inverse_Li_to_Lianat.mat
* li_modified_affine.nii => Li with modified affine (from linear registration above)
* wli_modified_affine.nii => Li mri Li coil in the MNI space.
* wt1_weighted_sagittal_1_0iso.nii => Anat mri H coil in the MNI space.
* sanlm_wli_modified_affine.nii => Denoised Li mri Li coil in the MNI space.
* figure_li_anat_MNI.png => plot of Li mri Li coil (threshold) and anat mri H coil
* figure_denoised_li_anat_MNI.png => plot of denoised Li mri Li coil (threshold) and anat mri H coil 	 
	 
	 






