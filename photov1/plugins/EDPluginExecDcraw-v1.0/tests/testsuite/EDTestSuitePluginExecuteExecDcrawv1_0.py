#
#    Project: Photo-v1.0
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) ESRF, Grenoble
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
__copyright__ = "ESRF, Grenoble"

from EDTestSuite  import EDTestSuite

class EDTestSuitePluginExecuteExecDcrawv1_0(EDTestSuite):
    """
    This is the test suite for execution for EDNA plugin Dcrawv1_0 
    It will run subsequently all unit tests and execution tests.     
    """

    def process(self):
        self.addTestCaseFromName("EDTestCasePluginExecuteExecDcrawv1_0")
        self.addTestCaseFromName("EDTestCasePluginExecuteExecDcrawv1_0_thumbnail")
        self.addTestCaseFromName("EDTestCasePluginExecuteExecDcrawv1_0_autoWB")
        self.addTestCaseFromName("EDTestCasePluginExecuteExecDcrawv1_0_tiff")
        self.addTestCaseFromName("EDTestCasePluginExecuteExecDcrawv1_0_interpol")

if __name__ == '__main__':

    edTestSuitePluginExecuteExecDcrawv1_0 = EDTestSuitePluginExecuteExecDcrawv1_0("EDTestSuitePluginExecuteExecDcrawv1_0")
    edTestSuitePluginExecuteExecDcrawv1_0.execute()

