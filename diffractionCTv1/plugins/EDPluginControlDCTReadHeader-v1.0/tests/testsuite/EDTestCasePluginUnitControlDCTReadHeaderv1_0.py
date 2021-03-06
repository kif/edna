#
#    Project: DiffractionCTv1
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) 2008-2009 European Synchrotron Radiation Facility
#                            Grenoble, France
#
#    Principal authors: Jerome Kieffer (kieffer@esrf.fr)
#                       Olof Svensson (svensson@esrf.fr) 
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


from EDVerbose import EDVerbose
from EDTestCasePluginUnit import EDTestCasePluginUnit
from XSDataDiffractionCTv1 import XSDataInputReadHeader

class EDTestCasePluginUnitControlDCTReadHeaderv1_0(EDTestCasePluginUnit):
    """
    """

    def __init__(self, _edStringTestName=None):
        """
        """
        EDTestCasePluginUnit.__init__(self, "EDPluginControlDCTReadHeaderv1_0", "EDPluginControlDCTReadHeader-v1.0", _edStringTestName)


    def testCheckParameters(self):
        xsDataInput = XSDataInputReadHeader()
        edPluginExecDCTReadHeader = self.createPlugin()
        edPluginExecDCTReadHeader.setDataInput(xsDataInput)
        edPluginExecDCTReadHeader.checkParameters()



    def process(self):
        self.addTestMethod(self.testCheckParameters)



##############################################################################


if __name__ == '__main__':

    EDTestCasePluginUnitControlDCTReadHeaderv1_0 = EDTestCasePluginUnitControlDCTReadHeaderv1_0("EDTestCasePluginUnitControlDCTReadHeaderv1_0")
    EDTestCasePluginUnitControlDCTReadHeaderv1_0.execute()
