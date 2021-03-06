#
#    Project: PROJECT
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
#
#    Principal author:       Jerome Kieffer (Jerome.Kieffer@ESRF.eu)
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
__author__ = "Jerome Kieffer (Jerome.Kieffer@ESRF.eu)"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"


from EDTestSuite  import EDTestSuite

class EDTestSuitePluginExecCommandLinev10(EDTestSuite):
    """
    This is the test suite for EDNA plugin CommandLinev10 
    It will run subsequently all unit tests and execution tests. 
    
    """

    def process(self):
        """
        """
        self.addTestCaseFromName("EDTestCasePluginUnitExecCommandLinev10")
        self.addTestCaseFromName("EDTestCasePluginUnitExecRsync")
        self.addTestCaseFromName("EDTestCasePluginExecuteExecCommandLinev10")
        self.addTestCaseFromName("EDTestCasePluginExecuteExecCommandLinev10_fireAndForget")



##############################################################################
if __name__ == '__main__':

    edTestSuitePluginExecCommandLinev10 = EDTestSuitePluginExecCommandLinev10("EDTestSuitePluginExecCommandLinev10")
    edTestSuitePluginExecCommandLinev10.execute()

##############################################################################
