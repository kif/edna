# coding: utf-8
#
#    Project: PROJECT
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) EMBL + ESRF + DLS
#
#    Principal author:       Al, Irakli, Alun, Jerome, Olof, Peter and Claudio
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

__author__ = "Al, Irakli, Alun, Jerome, Olof, Peter and Claudio"
__license__ = "GPLv3+"
__copyright__ = "EMBL + ESRF + DLS"

from EDVerbose import EDVerbose
from EDTestCasePluginUnit import EDTestCasePluginUnit

from XSDataEdnaSaxs import XSDataInputAutoRg, XSDataSaxsSample, XSDataFile

class EDTestCasePluginUnitExecAutoRgv1_0(EDTestCasePluginUnit):
    """
    Those are all units tests for the EDNA Exec plugin AutoRgv1_0
    """

    def __init__(self, _strTestName=None):
        """
        """
        EDTestCasePluginUnit.__init__(self, "EDPluginExecAutoRgv1_0")


    def testCheckParameters(self):
        xsDataInput = XSDataInputAutoRg()
        xsDataInput.sample = XSDataSaxsSample()
        xsDataInput.inputCurve = [XSDataFile()]
        edPluginExecAutoRg = self.createPlugin()
        edPluginExecAutoRg.setDataInput(xsDataInput)
        edPluginExecAutoRg.checkParameters()



    def process(self):
        self.addTestMethod(self.testCheckParameters)



if __name__ == '__main__':

    edTestCasePluginUnitExecAutoRgv1_0 = EDTestCasePluginUnitExecAutoRgv1_0("EDTestCasePluginUnitExecAutoRgv1_0")
    edTestCasePluginUnitExecAutoRgv1_0.execute()
