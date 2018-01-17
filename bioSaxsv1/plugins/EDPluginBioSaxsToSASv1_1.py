# coding: utf-8
#
#    Project: BioSaxs
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) 2011 ESRF
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
from __future__ import with_statement
__author__ = "Jérôme Kieffer"
__license__ = "GPLv3+"
__copyright__ = "2013 ESRF"
__status__ = "stable"
__date__ = "27/08/2015"

import os, shutil
from EDPluginControl        import EDPluginControl
from EDFactoryPlugin        import edFactoryPlugin
from EDUtilsPath            import EDUtilsPath
from EDUtilsPlatform        import EDUtilsPlatform
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
edFactoryPlugin.loadModule("XSDataWaitFilev1_0")
edFactoryPlugin.loadModule("XSDataExecCommandLine")
edFactoryPlugin.loadModule("XSDataEdnaSaxs")
from XSDataBioSaxsv1_0      import XSDataInputBioSaxsToSASv1_0, XSDataResultBioSaxsToSASv1_0, XSDataInputBioSaxsISPyBModellingv1_0
from XSDataWaitFilev1_0     import XSDataInputWaitFile
from XSDataExecCommandLine  import XSDataInputRsync
from XSDataCommon           import XSDataInteger, XSDataString, XSDataFile, XSPluginItem, XSDataStatus
from XSDataEdnaSaxs         import XSDataInputSaxsAnalysisModeling, XSDataInputSaxsModeling

architecture = EDUtilsPlatform.architecture
numpyPath = os.path.join(EDUtilsPath.EDNA_HOME, "libraries", "20090405-Numpy-1.3", architecture)
numpy = edFactoryPlugin.preImport("numpy", numpyPath)


class EDPluginBioSaxsToSASv1_1(EDPluginControl):
    """
	Replacement of Irakly's pipeline for SAS with a new implementation
    """
    CONF_MAX_THREAD = "maxThreads"
    CONF_FILE_SIZE = "fileSize"
    size = None
    maxThreads = None

    cpWait = "EDPluginWaitFile"
    cpModeling = "EDPluginControlSaxsModelingv1_1"
    cpAnalysisModeling = "EDPluginControlSaxsAnalysisModelingv1_0"
    cpISPyB = "EDPluginBioSaxsISPyBModellingv1_0"
    cpRsync = "EDPluginExecRsync"

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsToSASv1_0)
        self.dataOutput = XSDataResultBioSaxsToSASv1_0()
        self.pluginWait = None
        self.pluginModeling = None
        self.pluginRsync = None
        self.pluginISPyB = None
        self.inputFile = None
        self.strInFile = None
        self.gnomFile = None
        self.outFile = None
        self.wd = None
        self.xsdIspybInput = XSDataInputBioSaxsISPyBModellingv1_0()

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsToSASv1_1.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")

    def configure(self):
        """
        Configures the plugin from the configuration file with the following parameters:
        - CONF_MAX_THREAD max number of threads in the action cluster
        - CONF_FILE_SIZE minimum file size 
        """
        EDPluginControl.configure(self)
        self.DEBUG("EDPluginBioSaxsToSASv1_1.configure")

        if (self.__class__.size is None) or (self.__class__.maxThreads is None):
            xsPluginItem = self.getConfiguration()
            if (xsPluginItem == None):
                self.warning("EDPluginBioSaxsToSASv1_1.configure: No plugin item defined.")
                xsPluginItem = XSPluginItem()
            conf_size = self.config.get(self.CONF_FILE_SIZE, None)
            if conf_size is None:
                strMessage = 'EDPluginBioSaxsToSASv1_1.configure: Configuration parameter missing: \
    %s for %s, defaulting to "1000"' % (self.CONF_FILE_SIZE, EDUtilsPath.EDNA_SITE)
                self.WARNING(strMessage)
                self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                self.__class__.size = 1000
            else:
                self.__class__.size = int(conf_size)
            maxThreads = self.config.get(self.CONF_MAX_THREAD, None)
            if maxThreads is None:
                strMessage = 'EDPluginBioSaxsToSASv1_1.configure: Configuration parameter missing: \
    %s for %s, defaulting to "max"' % (self.CONF_MAX_THREAD, EDUtilsPath.EDNA_SITE)
                self.WARNING(strMessage)
                self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                self.__class__.maxThreads = None
            else:
                self.__class__.maxThreads = int(maxThreads)

    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginBioSaxsToSASv1_1.preProcess")
        if self.dataInput.gnomFile is not None:
            self.gnomFile = self.dataInput.gnomFile.path.value
            self.inputFile = self.dataInput.gnomFile
        elif self.dataInput.subtractedCurve is not None:
            self.strInFile = self.dataInput.subtractedCurve.path.value
            self.inputFile = self.dataInput.subtractedCurve
        else:
            self.error("Neither gnomFile not subtractedCurve are present in the input datastructure")
            self.setFailure()
            return
        if self.dataInput.destinationDirectory is not None:
            self.strWorkingDirectory = self.dataInput.destinationDirectory.path.value
        if self.dataInput.subtractedCurve is not None:
            self.strInFile = self.dataInput.subtractedCurve.path.value

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsToSASv1_1.process")
        self.pluginWait = self.loadPlugin(self.cpWait)
        self.pluginWait.dataInput = XSDataInputWaitFile(expectedSize=XSDataInteger(self.__class__.size),
                                                                expectedFile=self.inputFile)
        self.pluginWait.connectSUCCESS(self.doSuccessExecWait)
        self.pluginWait.connectFAILURE(self.doFailureExecWait)
        self.pluginWait.executeSynchronous()

        if self.isFailure():
            return

        if self.gnomFile and os.path.exists(self.gnomFile):
            self.pluginModeling = self.loadPlugin(self.cpModeling)
            self.pluginModeling.dataInput = XSDataInputSaxsModeling(graphFormat=XSDataString("png"),
                                                                    gnomFile=XSDataFile(XSDataString(self.gnomFile)))
        elif self.strInFile and os.path.exists(self.strInFile):
            self.pluginModeling = self.loadPlugin(self.cpAnalysisModeling)
            self.pluginModeling.dataInput = XSDataInputSaxsAnalysisModeling(graphFormat=XSDataString("png"),
                                                                            scatterCurve=XSDataFile(XSDataString(self.strInFile)))
        else:
            self.error("Neither gnomFile not subtractedCurve are present in the input datastructure")
            self.setFailure()
            return

        self.pluginModeling.connectSUCCESS(self.doSuccessExecSAS)
        self.pluginModeling.connectFAILURE(self.doFailureExecSAS)
        self.pluginModeling.executeSynchronous()

        ########################################################################
        # Send to ISPyB
        ########################################################################

        if self.dataInput.sample and self.dataInput.sample.login and \
                self.dataInput.sample.passwd and self.dataInput.sample.measurementID and \
                self.xsdIspybInput:
            self.addExecutiveSummaryLine("Registering to ISPyB")
            self.pluginISPyB = self.loadPlugin(self.cpISPyB)
            self.xsdIspybInput.sample = self.dataInput.sample
            self.pluginISPyB.dataInput = self.xsdIspybInput
            self.pluginISPyB.connectSUCCESS(self.doSuccessExecISPyB)
            self.pluginISPyB.connectFAILURE(self.doFailureExecISPyB)
            self.pluginISPyB.executeSynchronous()
        ########################################################################
        # Move results
        ########################################################################

        if self.dataInput.destinationDirectory is None:
            if self.strInFile is not None:
                outdir = os.path.join(os.path.dirname(os.path.dirname(self.strInFile)), "ednaSAS")
            else:
                #nothing specified
                outdir = os.path.join(os.path.dirname(os.getcwd()), "ednaSAS")
        else:
            outdir = self.dataInput.destinationDirectory.path.value
        if self.strInFile is not None:
            outdir = os.path.join(outdir, os.path.basename(os.path.splitext(self.strInFile)[0]))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        self.outFile = os.path.join(outdir, "NoResults.html")

        self.pluginRsync = self.loadPlugin(self.cpRsync)
        self.pluginRsync.dataInput = XSDataInputRsync(source=XSDataFile(XSDataString(self.wd)),
                                                      destination=XSDataFile(XSDataString(outdir)),
                                                      options=XSDataString("-avx"))

        self.pluginRsync.connectSUCCESS(self.doSuccessExecRsync)
        self.pluginRsync.connectFAILURE(self.doFailureExecRsync)
        self.pluginRsync.executeSynchronous()

        # if no errors up to now, clean up scratch disk
        if not self.isFailure():
            to_remove = self.pluginModeling.getWorkingDirectory()
            if os.path.isdir(to_remove):
                shutil.rmtree(to_remove)

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsToSASv1_1.postProcess")
        # Create some output data
        if self.outFile:
            self.dataOutput.htmlPage = XSDataFile(XSDataString(self.outFile))

    def finallyProcess(self, _edObject=None):
        EDPluginControl.finallyProcess(self, _edObject=_edObject)
        self.dataOutput.status = XSDataStatus(message=self.getXSDataMessage(),
                                          executiveSummary=XSDataString(os.linesep.join(self.getListExecutiveSummaryLines())))

    def doSuccessExecWait(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doSuccessExecWait")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doSuccessExecWait")
        self.retrieveMessages(_edPlugin)

    def doFailureExecWait(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doFailureExecWait")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doFailureExecWait")
        self.retrieveMessages(_edPlugin)
        self.setFailure()

    def doSuccessExecSAS(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doSuccessExecSAS")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doSuccessExecSAS")
        self.retrieveMessages(_edPlugin)
        self.xsdIspybInput.dammifModels = _edPlugin.dataOutput.dammifModels
        self.xsdIspybInput.damaverModel = _edPlugin.dataOutput.damaverModel
        self.xsdIspybInput.damfiltModel = _edPlugin.dataOutput.damfiltModel
        self.xsdIspybInput.damstartModel = _edPlugin.dataOutput.damstartModel
        self.xsdIspybInput.damminModel = _edPlugin.dataOutput.damminModel
        self.xsdIspybInput.fitFile = _edPlugin.dataOutput.fitFile
        self.xsdIspybInput.logFile = _edPlugin.dataOutput.logFile
        self.xsdIspybInput.pdbMoleculeFile = _edPlugin.dataOutput.pdbMoleculeFile
        self.xsdIspybInput.pdbSolventFile = _edPlugin.dataOutput.pdbSolventFile
        self.xsdIspybInput.chiRfactorPlot = _edPlugin.dataOutput.chiRfactorPlot
        self.xsdIspybInput.nsdPlot = _edPlugin.dataOutput.nsdPlot
#         dammifModels: XSDataSaxsModel [] optional
#         damaverModel: XSDataSaxsModel  optional
#         damfiltModel: XSDataSaxsModel  optional
#         damstartModel: XSDataSaxsModel  optional
#         damminModel: XSDataSaxsModel  optional
#         fitFile: XSDataFile optional
#         logFile: XSDataFile optional
#         pdbMoleculeFile: XSDataFile optional
#         pdbSolventFile: XSDataFile optional
#         chiRfactorPlot: XSDataFile optional
#         nsdPlot: XSDataFile optional
        self.wd = os.path.join(_edPlugin.getWorkingDirectory(), "")

    def doFailureExecSAS(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doFailureExecSAS")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doFailureExecSAS")
        self.retrieveMessages(_edPlugin)
        self.wd = os.path.join(_edPlugin.getWorkingDirectory(), "")
        self.setFailure()

    def doSuccessExecISPyB(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doSuccessExecISPyB")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doSuccessExecISPyB")
        self.retrieveMessages(_edPlugin)

    def doFailureExecISPyB(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doFailureExecISPyB")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doFailureExecISPyB")
        self.retrieveMessages(_edPlugin)
        self.setFailure()

    def doSuccessExecRsync(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doSuccessExecRsync")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doSuccessExecRsync")
        self.retrieveMessages(_edPlugin)

    def doFailureExecRsync(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsToSASv1_1.doFailureExecRsync")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsToSASv1_1.doFailureExecRsync")
        self.retrieveMessages(_edPlugin)
        self.setFailure()
