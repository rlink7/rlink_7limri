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
from .info import __version__, LICENSE, AUTHOR, DESCRIPTION


def logo():
    """ Module logo.

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


def info():
    """ Dispaly some usefull information about the package.

    Returns
    -------
    info: str
        package information.
    """
    rchar = "\n"
    version = f"Package version: {__version__.strip(rchar)}\n"
    license = f"License: {LICENSE.strip(rchar)}\n"
    authors = f"Authors: {AUTHOR.strip(rchar)}\n"
    desc = f"Description: {DESCRIPTION.strip(rchar)}\n"
    return logo() + "\n" + desc + version + license + authors
