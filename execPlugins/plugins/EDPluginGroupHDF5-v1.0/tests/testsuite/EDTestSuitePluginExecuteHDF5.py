#
#    Project: PROJECT
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) ESRF Grenoble
#
#    Principal author:       Jerome Kieffer
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

__author__ = "Jerome Kieffer"
__license__ = "GPLv3+"
__copyright__ = "ESRF Grenoble"

import os, sys
from EDTestSuite  import EDTestSuite
from EDVerbose                  import EDVerbose
from EDUtilsLibraryInstaller    import EDUtilsLibraryInstaller, installLibrary
from EDFactoryPluginStatic      import EDFactoryPluginStatic

EDFactoryPluginStatic.loadModule("EDInstallNumpyv1_3")
EDFactoryPluginStatic.loadModule("EDInstallPILv1_1_7")
EDFactoryPluginStatic.loadModule("EDInstallFabio_v0_0_6")
EDFactoryPluginStatic.loadModule("EDInstallH5Pyv1_3_0")

################################################################################
# EDNA_SITE is not needed for this plugin so why bother !
################################################################################
if not "EDNA_SITE" in os.environ:
    os.environ["EDNA_SITE"] = "edna"



class EDTestSuitePluginExecuteHDF5(EDTestSuite):
    """
    This is the test suite for EDNA plugin HDF5StackImagesv10 
    It will run subsequently all unit tests and execution tests.     
    """

    def process(self):
        self.addTestSuiteFromName("EDTestSuitePluginExecuteHDF5StackImagesv10")
        self.addTestCaseFromName("EDTestCasePluginExecuteHDF5MapSpectrav10")

#EDTestSuitePluginExecuteHDF5MapSpectra
#EDTestSuitePluginExecuteHDF5StackImagesv10

if __name__ == '__main__':

    edTestSuitePluginExecuteHDF5 = EDTestSuitePluginExecuteHDF5("EDTestSuitePluginHDF5StackImagesv10")
    edTestSuitePluginExecuteHDF5.execute()

