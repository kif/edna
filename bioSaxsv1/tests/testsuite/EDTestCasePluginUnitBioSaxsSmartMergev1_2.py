# coding: utf-8
#
#    Project: templatev1
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) 2011, ESRF Grenoble
#
#    Principal author:        Jérôme Kieffer
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

__author__ = "Jérôme Kieffer"
__license__ = "GPLv3+"
__copyright__ = "2011, ESRF Grenoble"


from EDVerbose import EDVerbose
from EDTestCasePluginUnit import EDTestCasePluginUnit
from EDFactoryPluginStatic import EDFactoryPluginStatic
EDFactoryPluginStatic.loadModule("XSDataEdnaSaxs")
from XSDataBioSaxsv1_0 import XSDataInputBioSaxsSmartMergev1_0, XSDataFile, XSDataDouble

class EDTestCasePluginUnitBioSaxsSmartMergev1_2(EDTestCasePluginUnit):


    def __init__(self, _edStringTestName=None):
        EDTestCasePluginUnit.__init__(self, "EDPluginBioSaxsSmartMergev1_2")


    def testCheckParameters(self):
        """
      complex type XSDataInputBioSaxsSmartMergev1_2 extends XSDataInput  {
    inputFile: XSDataFile []
    absoluteFidelity: XSDataDouble optional 
    relativeFidelity: XSDataDouble optional
    mergedCurve: XSDataFile
        """
        xsDataInput = XSDataInputBioSaxsSmartMergev1_0()
        xsDataInput.inputCurves = [XSDataFile()]
        xsDataInput.absoluteFidelity = XSDataDouble()
        xsDataInput.relativeFidelity = XSDataDouble()
        xsDataInput.mergedCurve = XSDataFile()
        edPluginExecBioSaxsSmartMerge = self.createPlugin()
        edPluginExecBioSaxsSmartMerge.setDataInput(xsDataInput)
        edPluginExecBioSaxsSmartMerge.checkParameters()



    def process(self):
        self.addTestMethod(self.testCheckParameters)




if __name__ == '__main__':

    EDTestCasePluginUnitBioSaxsSmartMergev1_2 = EDTestCasePluginUnitBioSaxsSmartMergev1_2("EDTestCasePluginUnitBioSaxsSmartMergev1_2")
    EDTestCasePluginUnitBioSaxsSmartMergev1_2.execute()
