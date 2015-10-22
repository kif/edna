# coding: utf-8
#
#    Project: BioSaxs
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) 2011, ESRF
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
__copyright__ = "2011, ESRF"


from EDVerbose import EDVerbose
from EDTestCasePluginUnit import EDTestCasePluginUnit

from XSDataBioSaxsv1_0 import XSDataInputBioSaxsSingleSamplev1_0, \
    XSDataFile, XSDataFileSeries

class EDTestCasePluginUnitBioSaxsSingleSamplev1_0(EDTestCasePluginUnit):


    def __init__(self, _edStringTestName=None):
        EDTestCasePluginUnit.__init__(self, "EDPluginBioSaxsSingleSamplev1_0")


    def testCheckParameters(self):
        xsDataInput = XSDataInputBioSaxsSingleSamplev1_0()
        xsDataInput.directory1D = XSDataFile()
        xsDataInput.directory2D = XSDataFile()
        xsDataInput.directoryMisc = XSDataFile()
        xsDataInput.bufferSeries = [XSDataFileSeries()]
        xsDataInput.sampleSeries = [XSDataFileSeries()]
        edPluginExecBioSaxsSingleSamplev1_0 = self.createPlugin()
        edPluginExecBioSaxsSingleSamplev1_0.setDataInput(xsDataInput)
        edPluginExecBioSaxsSingleSamplev1_0.checkParameters()



    def process(self):
        self.addTestMethod(self.testCheckParameters)




if __name__ == '__main__':

    EDTestCasePluginUnitBioSaxsSingleSamplev1_0 = EDTestCasePluginUnitBioSaxsSingleSamplev1_0("EDTestCasePluginUnitBioSaxsSingleSamplev1_0")
    EDTestCasePluginUnitBioSaxsSingleSamplev1_0.execute()
