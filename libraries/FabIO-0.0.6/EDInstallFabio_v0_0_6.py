#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Project: EDNA Libraries
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
#
#    Principal authors:   Olof Svensson (svensson@esrf.fr)
#                         Jerome Kieffer (Jerome.Kieffer@ESRF.eu)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#

"""EDNA installer for Fabio v0.0.6"""

__authors__ = ["Olof Svensson", "Jerome Kieffer"]
__contact__ = "svensson@esrf.fr"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "20110722"
import os, sys
if "EDNA_HOME" not in os.environ:
    EDNA_HOME = os.path.dirname(os.path.dirname(__file__))
    os.environ["EDNA_HOME"] = EDNA_HOME
else:
    EDNA_HOME = os.environ["EDNA_HOME"]
kernel_src = os.path.join(EDNA_HOME, "kernel", "src")
if kernel_src not in sys.path:
    sys.path.append(kernel_src)

from EDVerbose                  import EDVerbose
from EDUtilsPlatform            import EDUtilsPlatform
from EDUtilsLibraryInstaller    import EDUtilsLibraryInstaller, installLibrary
from EDFactoryPluginStatic      import EDFactoryPluginStatic

moduleName = "fabio"
modulePath = os.path.join(os.environ["EDNA_HOME"], "libraries", "FabIO-0.0.6", EDUtilsPlatform.architecture)
moduleVersion = "0.0.6"
################################################################################
# Import the good version of FabIO ... 
################################################################################

oModule = EDFactoryPluginStatic.preImport(moduleName)
if not oModule:
    oModule = EDFactoryPluginStatic.preImport(moduleName, modulePath)
    if oModule is None:
        installLibrary(modulePath)
        oModule = EDFactoryPluginStatic.preImport(moduleName, modulePath)
try:
    version = oModule.version
except AttributeError:
    version = "0.0.0"

if version.split(".") < moduleVersion.split("."):
    EDVerbose.screen("Wrong %s library:  %s " % (moduleName, version))
    EDFactoryPluginStatic.unImport(moduleName)
    oModule = EDFactoryPluginStatic.preImport(moduleName, modulePath, moduleVersion)

if oModule is None:
    EDVerbose.ERROR("Unable to download, compile or install module %s" % moduleName)
else:
    EDVerbose.screen("Version of %s: %s from %s" % (moduleName, oModule.version, oModule.__file__))
