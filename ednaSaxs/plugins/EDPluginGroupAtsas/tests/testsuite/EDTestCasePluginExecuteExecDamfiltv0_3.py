#
#    Project: PROJECT
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) DLS
#
#    Principal author:       irakli
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

__author__ = "irakli & Jerome Kieffer"
__license__ = "GPLv3+"
__copyright__ = "2010 DLS, 2014 ESRF"

import os

from EDTestCasePluginExecute             import EDTestCasePluginExecute
from XSDataCommon                    import XSDataString
from parse_atsas import get_ATSAS_version


class EDTestCasePluginExecuteExecDamfiltv0_3(EDTestCasePluginExecute):
    """
    Those are all execution tests for the EDNA Exec plugin Damfiltv0_3
    """

    def __init__(self, _strTestName=None):
        """
        """
        EDTestCasePluginExecute.__init__(self, "EDPluginExecDamfiltv0_3")
#        self.setConfigurationFile(os.path.join(self.getPluginTestsDataHome(),
#                                               "XSConfiguration_Damfilt.xml"))
        self.setConfigurationFile(self.getRefConfigFile())
        self.setDataInputFile(os.path.join(self.getPluginTestsDataHome(), \
                                           "XSDataInputDamfilt_reference.xml"))
        self.setReferenceDataOutputFile(os.path.join(self.getPluginTestsDataHome(), \
                                                     "XSDataResultDamfilt_reference.xml"))


    def testExecute(self):
        """
        """
        inputPdbPath = XSDataString()

        # Add path to the input data file
        inputPdbName = self.getPlugin().getDataInput().getInputPdbFile().getPath().getValue()
        inputPdbPath.setValue(os.path.join(self.getPluginTestsDataHome(), inputPdbName))
        self.getPlugin().getDataInput().getInputPdbFile().setPath(inputPdbPath)

        self.run()


    def process(self):
        """
        """
        self.addTestMethod(self.testExecute)



if __name__ == '__main__':

    testDamfiltv0_3instance = EDTestCasePluginExecuteExecDamfiltv0_3("EDTestCasePluginExecuteExecDamfiltv0_3")
    testDamfiltv0_3instance.execute()
