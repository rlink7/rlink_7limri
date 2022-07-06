# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2021
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# Module current version
version_major = 0
version_minor = 0
version_micro = 0

# Expected by setup.py: string of form "X.Y.Z"
__version__ = "{0}.{1}.{2}".format(version_major, version_minor, version_micro)

# Expected by setup.py: the status of the project
CLASSIFIERS = ["Development Status :: 5 - Production/Stable",
               "Environment :: Console",
               "Environment :: X11 Applications :: Qt",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering",
               "Topic :: Utilities"]

# Project descriptions
description = """
Package that provides tools for brain Lithium MRI PreProcessing.
"""
SUMMARY = """
.. container:: summary-carousel

    `Limri` is a toolbox that provides common lithium brain
    MRI PreProcessing scripts:

    * preproc
    * normalisation
    * quantification
"""
long_description = (
    "Package that provides tools for brain Lithium MRI PreProcessing.\n")

# Main setup parameters
NAME = "Limri"
ORGANISATION = "CEA"
MAINTAINER = "Julie Victor"
MAINTAINER_EMAIL = "julie.victor@cea.fr"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
EXTRANAME = "NeuroSpin webPage"
EXTRAURL = (
    "https://joliot.cea.fr/drf/joliot/Pages/Entites_de_recherche/"
    "NeuroSpin.aspx")
LINKS = {"scripts": "https://github.com/neurospin/scripts"}
URL = "https://github.com/neurospin/scripts/2017_bipli_lithium_imaging"
DOWNLOAD_URL = \
    "https://github.com/neurospin/scripts/2017_bipli_lithium_imaging"
LICENSE = "CeCILL-B"
AUTHOR = """
Limri developers
"""
AUTHOR_EMAIL = "julie.victor@cea.fr"
PLATFORMS = "OS Independent"
ISRELEASE = True
VERSION = __version__
PROVIDES = ["Limri"]
REQUIRES = [
    "numpy>=1.17.1",
    "nibabel>=2.4.0"
]
SCRIPTS = [
    "Limri/scripts/preproc",
    "Limri/scripts/normalisation",
    "Limri/scripts/quantification"
]
