# coding: utf-8
#
#    Project: BioSaxs
#             http://www.edna-site.org
#
#    Copyright (C) 2012 ESRF
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
__copyright__ = "2012 ESRF"
__date__ = "01/03/2018"
__status__ = "development"


import os
import traceback

from EDPluginControl import EDPluginControl
from EDFactoryPlugin import edFactoryPlugin
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
edFactoryPlugin.loadModule("XSDataEdnaSaxs")

from XSDataBioSaxsv1_0 import XSDataInputBioSaxsHPLCv1_0, XSDataResultBioSaxsHPLCv1_0, \
                            XSDataInputBioSaxsISPyB_HPLCv1_0, XSDataInputBioSaxsToSASv1_0, \
                            XSDataInputBioSaxsISPyBv1_0, XSDataInputBioSaxsISPyBHPLCv1_0
from XSDataEdnaSaxs import XSDataInputDataver, XSDataInputSaxsAnalysis

from XSDataCommon import XSDataString, XSDataStatus, XSDataFile
from EDPluginBioSaxsHPLCv1_5 import EDPluginBioSaxsHPLCv1_5


class EDPluginBioSaxsFlushHPLCv1_5 (EDPluginControl):
    """
    plugin that just flushes the HPLC data to disk
    
    News: 
    v1.1: 
    * adapt for version HPLC plugin v1.1
    * ISPyB
    
    v1.2:
    * adapt for version HPLC plugin v1.2
    
    v1.5: Switch autoRG: atsas -> freesas
    """
    strControlledPluginDatAver = "EDPluginExecDataverv2_0"
    strControlledPluginISPyB = "EDPluginBioSaxsISPyB_HPLCv1_0"
    __strControlledPluginSaxsAnalysis = "EDPluginControlSaxsAnalysisv1_0"
    __strControlledPluginSaxsModeling = "EDPluginBioSaxsToSASv1_1"
    __strControlledPluginISPyBAnalysis = "EDPluginHPLCPrimayDataISPyBv1_0"

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsHPLCv1_0)
        self.xsDataResult = XSDataResultBioSaxsHPLCv1_0()
        self.runId = None
        self.FrameId = None
        self.hplc_run = None
        self.curve = None
        self.subtracted = None
        self.lstExecutiveSummary = []
        self.modeling = True

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.rawImage, "No raw image")

    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.preProcess")
        sdi = self.dataInput
        if sdi.runId is not None:
            self.runId = sdi.runId.value
        else:
            path = sdi.rawImage.path.value
            if "_" in path:
                self.runId = path[::-1].split("_", 1)[1][::-1]
            else:
                self.runId = path

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.process")
        if self.runId in EDPluginBioSaxsHPLCv1_5.dictHPLC:
            self.processRun(EDPluginBioSaxsHPLCv1_5.dictHPLC[self.runId])
            try:
                edpluginIsPyB = self.loadPlugin(self.strControlledPluginISPyB)
                edpluginIsPyB.dataInput = XSDataInputBioSaxsISPyB_HPLCv1_0(sample=self.dataInput.sample,
                                                                           hdf5File=self.xsDataResult.hplcFile,
                                                                           jsonFile=XSDataFile(XSDataString(self.json)),
                                                                           hplcPlot=self.xsDataResult.hplcImage)
                edpluginIsPyB.executeSynchronous()
                self.dataOutputBioSaxsISPyB_HPLC = edpluginIsPyB.xsdResult
            except Exception as error:
                traceback.print_stack()
                self.ERROR("EDPluginBioSaxsFlushHPLCv1_5 calling to EDPluginBioSaxsISPyB_HPLCv1_0: %s" % error)

            self.processMerges(EDPluginBioSaxsHPLCv1_5.dictHPLC[self.runId])

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.postProcess")
        self.synchronizePlugins()

    def finallyProcess(self, _edObject=None):
        EDPluginControl.finallyProcess(self)
        executiveSummary = os.linesep.join(self.lstExecutiveSummary)
        self.xsDataResult.status = XSDataStatus(executiveSummary=XSDataString(executiveSummary))
        self.dataOutput = self.xsDataResult

    def processRun(self, run):
        # EDPluginBioSaxsHPLCv1_5.dictHPLC[self.runId].reset()
        for idx in run.frames:
            run.frames[idx].purge_memory()
        run.dump_json()

        hdf5 = run.save_hdf5(name=self.__class__.__name__)
        self.json = os.path.splitext(hdf5)[0] + ".json"
        self.xsDataResult.hplcFile = XSDataFile(XSDataString(hdf5))
        self.xsDataResult.hplcImage = XSDataFile(XSDataString(run.make_plot()))
        peaks = []
        try:
            peaks = run.analyse()
        except Exception as error:
            traceback.print_stack()
            self.ERROR("EDPluginBioSaxsFlushHPLCv1_5:  Error in run.analyse() %s" % error)
        for group in peaks:
                self.lstExecutiveSummary.append("Merging frames from %s to %s" % (group[0], group[-1]))
                xsdFrames = [XSDataFile(XSDataString(run.frames[i].subtracted)) for i in group]
                outname = os.path.splitext(run.frames[group[0]].subtracted)[0] + "_aver_%s.dat" % group[-1]
                edpugin = self.loadPlugin(self.strControlledPluginDatAver)
                edpugin.dataInput = XSDataInputDataver(outputCurve=XSDataFile(XSDataString(outname)), inputCurve=xsdFrames)
                edpugin.connectSUCCESS(self.doSuccessDatAver)
                edpugin.connectFAILURE(self.doFailureDatAver)
                edpugin.execute()
                run.merge_curves.append(outname)
                run.merge_analysis[outname] = None
                run.merge_Rg[outname] = None
                run.merge_framesDIC[outname] = [group[0], group[-1]]
        # Append to hdf5
        run.append_hdf5()

    def processMerges(self, run):
        # run analysis of merges

        for merge in run.merge_curves:
            if os.path.exists(merge):
                print merge
                xsdSubtractedCurve = XSDataFile(XSDataString(merge))
                self.__edPluginSaxsAnalysis = self.loadPlugin(self.__strControlledPluginSaxsAnalysis)

                self.__edPluginSaxsAnalysis.dataInput = XSDataInputSaxsAnalysis(

                                                                                    scatterCurve=xsdSubtractedCurve,
                                                                                    autoRg=run.merge_Rg[merge],
                                                                                    graphFormat=XSDataString("png")
                                                                                )
                self.__edPluginSaxsAnalysis.connectSUCCESS(self.doSuccessSaxsAnalysis)
                self.__edPluginSaxsAnalysis.connectFAILURE(self.doFailureSaxsAnalysis)
                self.__edPluginSaxsAnalysis.executeSynchronous()
                xsdBuffer = XSDataFile(XSDataString(run.buffer))
                xsdStartFrame = XSDataString(run.merge_framesDIC[merge][0])
                xsdEndFrame = XSDataString(run.merge_framesDIC[merge][-1])
                self.__edPluginISPyBAnalysis = self.loadPlugin(self.__strControlledPluginISPyBAnalysis)
                try:
                    inputBioSaxsISPyB = XSDataInputBioSaxsISPyBv1_0(
                                                                sample=self.dataOutputBioSaxsISPyB_HPLC.getSample(),
                                                                autoRg=run.merge_analysis[merge].autoRg,
                                                                gnom=run.merge_analysis[merge].gnom,
                                                                volume=run.merge_analysis[merge].volume,
                                                                bestBuffer=xsdBuffer,
                                                                scatterPlot=run.merge_analysis[merge].scatterPlot,
                                                                guinierPlot=run.merge_analysis[merge].guinierPlot,
                                                                kratkyPlot=run.merge_analysis[merge].kratkyPlot,
                                                                densityPlot=run.merge_analysis[merge].densityPlot
                                                )

                    xsdISPyBin = XSDataInputBioSaxsISPyBHPLCv1_0(
                                                                experimentId=self.dataOutputBioSaxsISPyB_HPLC.getExperimentId(),
                                                                startFrame=xsdStartFrame,
                                                                endFrame=xsdEndFrame,
                                                                dataInputBioSaxs=inputBioSaxsISPyB
                                                                )
                    self.__edPluginISPyBAnalysis.dataInput = xsdISPyBin
                    self.__edPluginISPyBAnalysis.connectSUCCESS(self.doSuccessISPyBAnalysis)
                    self.__edPluginISPyBAnalysis.connectFAILURE(self.doFailureISPyBAnalysis)
                    self.__edPluginISPyBAnalysis.executeSynchronous()
                except Exception as error:
                    traceback.print_stack()
                    self.ERROR("EDPluginBioSaxsFlushHPLCv1_5 calling to EDPluginHPLCPrimayDataISPyBv1_0: %s" % error)

        self.synchronizePlugins()
        # There were some recurring issues with dammin slowing down slavia, therefore I commented this out for the time being
        # Martha, 11.7.2014
        mergeNumber = 1
        # print run.merge_curves
        if self.modeling == True:
            for merge in run.merge_curves:
                if run.merge_analysis[merge] is not None and 15.0 >= run.merge_analysis[merge].autoRg.rg.value >= 1.0:
                    try:
                        xsdSubtractedCurve = XSDataFile(XSDataString(merge))
                        xsdGnomFile = XSDataFile(XSDataString(run.merge_analysis[merge].gnom.gnomFile.path.value))
                        destination = XSDataFile(XSDataString(os.path.join(os.path.dirname(os.path.dirname(merge)), "ednaSAS")))
                        self.__edPluginSaxsToSAS = self.loadPlugin(self.__strControlledPluginSaxsModeling)
                        # print "Changing measurentID by runMerge"
                        # In order to keep dammin models in different folder a measurementId should be given
                        self.__edPluginISPyBAnalysis.xsDataResult.dataInputBioSaxs.sample.measurementID.value = mergeNumber
                        # print "------------>  MeasurementId changed " + str(self.__edPluginISPyBAnalysis.xsDataResult.dataInputBioSaxs.sample.measurementID.value)
                        self.__edPluginSaxsToSAS.dataInput = XSDataInputBioSaxsToSASv1_0(
                                                                                             sample=self.__edPluginISPyBAnalysis.xsDataResult.dataInputBioSaxs.sample,
                                                                                             subtractedCurve=xsdSubtractedCurve,
                                                                                             gnomFile=xsdGnomFile,
                                                                                             destinationDirectory=destination)
                        self.__edPluginSaxsToSAS.connectSUCCESS(self.doSuccessSaxsToSAS)
                        self.__edPluginSaxsToSAS.connectFAILURE(self.doFailureSaxsToSAS)
                        mergeNumber = mergeNumber + 1;
                        # self.__edPluginSaxsToSAS.executeSynchronous()
                        self.__edPluginSaxsToSAS.execute()
                    except Exception as error:
                        traceback.print_stack()
                        self.ERROR("EDPluginBioSaxsFlushHPLCv1_5 calling to EDPluginBioSaxsToSASv1_1: %s" % error)
        self.synchronizePlugins()

    def doSuccessDatAver(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doSuccessDatAver")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsFlushHPLCv1_5.doSuccessDatAver")
        if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.outputCurve:
            outCurve = _edPlugin.dataOutput.outputCurve.path.value
            if os.path.exists(outCurve):
                self.xsDataResult.mergedCurves += [_edPlugin.dataOutput.outputCurve]
            else:
                strErr = "DatAver claimed merged curve is in %s but no such file !!!" % outCurve
                self.ERROR(strErr)
                self.lstExecutiveSummary.append(strErr)
                self.setFailure()

    def doFailureDatAver(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doFailureDatAver")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsFlushHPLCv1_5.doFailureDatAver")
        if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
            self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
        else:
            self.lstExecutiveSummary.append("Edna plugin DatAver failed.")
        self.setFailure()

    def doSuccessSaxsAnalysis(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doSuccessSaxsAnalysis")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsFlushHPLCv1_5.doSuccessSaxsAnalysis")
        self.retrieveMessages(_edPlugin)
        run = EDPluginBioSaxsHPLCv1_5.dictHPLC[self.runId]
        curvename = _edPlugin.dataOutput.autoRg.filename.path.value
        run.merge_analysis[curvename] = _edPlugin.dataOutput
#         run.merge_xsScatterPlot[curvename] = _edPlugin.dataOutput.scatterPlot
#         run.merge_xsGuinierPlot[curvename] = _edPlugin.dataOutput.guinierPlot
#         run.merge_xsKratkyPlot[curvename] = _edPlugin.dataOutput.kratkyPlot
#         run.merge_xsDensityPlot[curvename] = _edPlugin.dataOutput.densityPlot
        self.addExecutiveSummaryLine(_edPlugin.dataOutput.status.executiveSummary.value)

    def doFailureSaxsAnalysis(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doFailureSaxsAnalysis")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsFlushHPLCv1_5.doFailureSaxsAnalysis")
        self.retrieveMessages(_edPlugin)
        strErr = "Error in Processing of EDNA SaxsAnalysis = AutoRg => datGnom => datPorod"
        self.ERROR(strErr)
        self.addExecutiveSummaryLine(strErr)
        self.setFailure()

    def doSuccessSaxsToSAS(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doSuccessSaxsToSAS")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsFlushHPLCv1_5.doSuccessSaxsToSAS")

    def doFailureSaxsToSAS(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doFailureSaxsToSAS")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsFlushHPLCv1_5.doFailureSaxsToSAS")
        self.retrieveMessages(_edPlugin)
        strErr = "Error in Modeling of merged data"
        self.ERROR(strErr)
        self.addExecutiveSummaryLine(strErr)
        self.setFailure()

    def doSuccessISPyBAnalysis(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doSuccessISPyBAnalysis")
        self.addExecutiveSummaryLine("Registered analysis in ISPyB")
        self.retrieveMessages(_edPlugin)

    def doFailureISPyBAnalysis(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsFlushHPLCv1_5.doFailureISPyBAnalysis")
        self.addExecutiveSummaryLine("Failed to register analysis in ISPyB")
        self.retrieveMessages(_edPlugin)
