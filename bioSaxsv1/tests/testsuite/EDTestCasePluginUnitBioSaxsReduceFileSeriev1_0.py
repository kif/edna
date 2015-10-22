# coding: utf-8
#
#    Project: BioSaxs
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

from XSDataBioSaxsv1_0 import XSDataInputBioSaxsReduceFileSeriev1_0, XSDataFileSeries, \
        XSDataFile, XSDataBioSaxsSample, XSDataBioSaxsExperimentSetup

class EDTestCasePluginUnitBioSaxsReduceFileSeriev1_0(EDTestCasePluginUnit):


    def __init__(self, _edStringTestName=None):
        EDTestCasePluginUnit.__init__(self, "EDPluginBioSaxsReduceFileSeriev1_0")


    def testCheckParameters(self):
        xsDataInput = XSDataInputBioSaxsReduceFileSeriev1_0()
        xsDataInput.fileSerie = XSDataFileSeries()
        xsDataInput.sample = XSDataBioSaxsSample()
        xsDataInput.experimentSetup = XSDataBioSaxsExperimentSetup()
        xsDataInput.directory1D = XSDataFile()
        xsDataInput.directory2D = XSDataFile()
        xsDataInput.directoryMisc = XSDataFile()
        edPluginExecBioSaxsReduceFileSerie = self.createPlugin()

        edPluginExecBioSaxsReduceFileSerie.setDataInput(xsDataInput)
        edPluginExecBioSaxsReduceFileSerie.checkParameters()



    def process(self):
        self.addTestMethod(self.testCheckParameters)




if __name__ == '__main__':

    EDTestCasePluginUnitBioSaxsReduceFileSeriev1_0 = EDTestCasePluginUnitBioSaxsReduceFileSeriev1_0("EDTestCasePluginUnitBioSaxsReduceFileSeriev1_0")
    EDTestCasePluginUnitBioSaxsReduceFileSeriev1_0.execute()
