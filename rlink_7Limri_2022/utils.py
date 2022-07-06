# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2019
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
A module with common functions.
"""

# System import
import logging
import os
import gzip
import re
import subprocess

# Third party imports


logger = logging.getLogger("Limri")


def logo():
    """ http://patorjk.com/software/taag/#p=display&f=Big&t=Limri

    Returns
    -------
    logo: str
        the logo.
    """
    logo = r"""
  _      _                _ 
 | |    (_)              (_)
 | |     _ _ __ ___  _ __ _ 
 | |    | | '_ ` _ \| '__| |
 | |____| | | | | | | |  | |
 |______|_|_| |_| |_|_|  |_|
                            """
    return logo


def ungzip_file(zfile, prefix="u", outdir=None):
    """ Copy and ungzip the input file.

    Parameters
    ----------
    zfile: str
        input file to ungzip.
    prefix: str, default 'u'
        the prefix of the result file.
    outdir: str, default None)
        the output directory where ungzip file is saved. If not set use the
        input image directory.

    Returns
    -------
    unzfile: str
        the ungzip file.
    """
    # Checks
    if not os.path.isfile(zfile):
        raise ValueError("'{0}' is not a valid filename.".format(zfile))
    if outdir is not None:
        if not os.path.isdir(outdir):
            raise ValueError("'{0}' is not a valid directory.".format(outdir))
    else:
        outdir = os.path.dirname(zfile)

    # Get the file descriptors
    base, extension = os.path.splitext(zfile)
    basename = os.path.basename(base)

    # Ungzip only known extension
    if extension in [".gz"]:
        basename = prefix + basename
        unzfile = os.path.join(outdir, basename)
        with gzip.open(zfile, "rb") as gzfobj:
            data = gzfobj.read()
        with open(unzfile, "wb") as openfile:
            openfile.write(data)

    # Default, unknown compression extension: the input file is returned
    else:
        unzfile = zfile

    return unzfile


def create_bids_filetree(filepath, dest_path):
    # write the results in bids format
    path_tmp = filepath.split(os.sep)
    filename = os.path.basename(filepath)
    if not re.search("ses-", filename):
        raise ValueError("session key is needed")
    ses = path_tmp[-3]
    sub = path_tmp[-4]
    sub_dest = os.path.join(dest_path, sub)
    ses_dest = os.path.join(sub_dest, ses)
    anatdest = os.path.join(ses_dest, "anat")
    # create filetree
    subprocess.check_call(['mkdir', '-p', anatdest])
    return 0


def create_symbolic_link(filepath, dest_path):
    newpath = os.path.dirname(filepath)
    newpath = newpath.split(os.sep)[-3:]
    newpath = os.sep.join(dest_path, newpath)
    cmd = ["ln", "-s", filepath, newpath]
    subprocess.check_call(cmd)
    return 0
