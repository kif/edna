# coding: utf-8
# 
#
#    Project: ExecPlugins
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) European Synchrotron Radiation Facility, Grenoble 2010
#
#    Principal author:       Jérôme Kieffer
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
__contact__ = "jerome.kieffer@esrf.fr"
__license__ = "GPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, 2010"

import os

from EDVerbose                           import EDVerbose
from EDAssert                            import EDAssert
from EDTestCasePluginExecute             import EDTestCasePluginExecute
from XSDataEDFv1_0 import XSDataResult1DPowderEDF


class EDTestCasePluginExecuteExportAsciiPowderEDF_edf2spr(EDTestCasePluginExecute):
    """
    Those are all execution tests for the EDNA Exec plugin <pluginName>
    """

    def __init__(self, _strTestName=None):
        """
        """
        EDTestCasePluginExecute.__init__(self, "EDPluginExportAsciiPowderv1_0")
#        self.setConfigurationFile(os.path.join(self.getPluginTestsDataHome(),
#                                               "XSConfiguration_<basePluginName>.xml"))
        self.setDataInputFile(os.path.join(self.getPluginTestsDataHome(), \
                                           "XSDataInput1DPowder_edf2spr.xml"))
        self.setReferenceDataOutputFile(os.path.join(self.getPluginTestsDataHome(), \
                                                     "XSDataResult1DPowder_edf2spr.xml"))
        self.refOutput = "moke-azim.spr"


    def preProcess(self):
        """
        """
        EDTestCasePluginExecute.preProcess(self)
        self.loadTestImage(["moke-azim.edf" , self.refOutput])
        tmpdir = "/tmp/edna-%s" % os.environ["USER"]
        if not os.path.isdir(tmpdir):
            os.mkdir(tmpdir)
        if os.path.isfile(os.path.join(tmpdir, self.refOutput)):
            os.remove(os.path.join(tmpdir, self.refOutput))


    def testExecute(self):
        """
        """
        self.run()
        plugin = self.getPlugin()

        # Checking obtained results
        xsDataResult = plugin.getDataOutput()
        EDAssert.strAlmostEqual(XSDataResult1DPowderEDF.parseString(self.readAndParseFile(self.getReferenceDataOutputFile())).marshal(), xsDataResult.marshal(), _strComment="XML structures are the same")

        outputData = open(xsDataResult.getOutputFile().getPath().getValue(), "rb").read()
        referenceData = open(os.path.join(self.getTestsDataImagesHome(), self.refOutput), "rb").read()

        EDAssert.strAlmostEqual(referenceData, outputData, "SpreadSheet files are the same")



    def process(self):
        """
        """
        self.addTestMethod(self.testExecute)



if __name__ == '__main__':

    testinstance = EDTestCasePluginExecuteExportAsciiPowderEDF_edf2spr("EDTestCasePluginExecuteExportAsciiPowderEDF_edf2spr")
    testinstance.execute()
