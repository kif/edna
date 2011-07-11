# coding: utf8
#
#    Project: The EDNA Kernel
#             http://www.edna-site.org
#
#    File: "$Id: EDTestCaseParallelExecute.py 2548 2010-12-01 09:14:01Z kieffer $"
#
#    Copyright (C) 2011-2011 European Synchrotron Radiation Facility
#                            Grenoble, France
#
#    Principal authors: Jérôme Kieffer (Jerome.Kieffer@esrf.eu)
#                        
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    and the GNU Lesser General Public License  along with this program.  
#    If not, see <http://www.gnu.org/licenses/>.
#

__authors__ = ["Jérôme Kieffer"]
__contact__ = "Jerome.Kieffer@esrf.eu"
__license__ = "LGPLv3+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

"""Test cases for testing EDShare"""

import  time, tempfile, os
from EDVerbose                           import EDVerbose
from EDTestCase                          import EDTestCase
from EDAssert                            import EDAssert
from EDShare                             import EDShare
if "USER" not in os.environ:
    os.environ["USER"] = "ednatester"

class EDTestCaseEDShare(EDTestCase):
    """
    Unit test for EDNA-share for sharing large objects between plugins
    """

    def unitTestInitialState(self):
        """
        Check the initialization of the share:
        """
        EDAssert.equal(False, EDShare.isInitialized(), "Check that EDShare is uninitialized")
        EDShare.initialize(os.path.join(tempfile.gettempdir(), "edna-%s" % os.environ["USER"]))
        EDShare["test1"] = range(10)
        EDVerbose.screen("Backend used is %s" % EDShare.backend)
        EDAssert.equal(1, len(EDShare.items()))
        EDAssert.equal(True, "test1" in EDShare, "list is actually present")
        for i, j in zip(range(10), EDShare["test1"]):
            EDAssert.equal(i, j, "elements are the same")
        EDShare.close()
        EDShare.initialize(os.path.join(tempfile.gettempdir(), "edna-%s" % os.environ["USER"]))
        EDAssert.equal(True, "test1" in EDShare, "list is still present")
        for i, j in zip(range(10), EDShare["test1"]):
            EDAssert.equal(i, j, "elements are still the same")
        filename = EDShare.filename
        EDShare.close(remove=True)
        EDAssert.equal(False, os.path.isfile(filename), "dump-file has been removed")


#    def unitTestRunning(self):
#        """
#        check the status after a job creation
#        """
#        EDVerbose.DEBUG("EDTestCaseEDShare.unitTestRunning")
#        EDShare.tellRunning(self.strPluginName)
#        EDVerbose.DEBUG("Success Plugins: " + ",".join(EDShare.getSuccess()))
#        EDVerbose.DEBUG("Running Plugins: " + ",".join(EDShare.getRunning()))
#        EDVerbose.DEBUG("Failed Plugins: " + ",".join(EDShare.getFailure()))
#        EDAssert.equal(True, self.strPluginName in EDShare.getRunning(), "Plugin  running")
#        EDAssert.equal(False, self.strPluginName in EDShare.getSuccess(), "Plugin not yet Finished")
#        EDAssert.equal(False, self.strPluginName in EDShare.getFailure(), "Plugin not yet Finished")
#
#    def unitTestFailed(self):
#        """
#        check the failure of a plugin is registerd 
#        """
#        EDVerbose.DEBUG("EDTestCaseEDShare.unitTestFailed")
#        EDShare.tellFailure(self.strPluginName)
#        EDVerbose.DEBUG("Success Plugins: " + ",".join(EDShare.getSuccess()))
#        EDVerbose.DEBUG("Running Plugins: " + ",".join(EDShare.getRunning()))
#        EDVerbose.DEBUG("Failed Plugins: " + ",".join(EDShare.getFailure()))
#        EDAssert.equal(False, self.strPluginName in EDShare.getRunning(), "Plugin not yet running")
#        EDAssert.equal(False, self.strPluginName in EDShare.getSuccess(), "Plugin not yet Finished")
#        EDAssert.equal(True, self.strPluginName in EDShare.getFailure(), "Plugin Failed as expected")

    def process(self):
        self.addTestMethod(self.unitTestInitialState)
#        self.addTestMethod(self.unitTestRunning)
#        self.addTestMethod(self.unitTestFailed)


if __name__ == '__main__':

    edTestCase = EDTestCaseEDShare("EDTestCaseEDShare")
    edTestCase.execute()