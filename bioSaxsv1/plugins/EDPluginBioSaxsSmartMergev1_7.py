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
from __future__ import with_statement

__author__ = "Jérôme Kieffer, Alejandro Antolinos, Martha Brennich"
__contact__ = "Jerome.Kieffer@esrf.fr"
__license__ = "GPLv3+"
__copyright__ = "2015, ESRF Grenoble, EMBL Grenoble"
__date__ = "20150117"
__status__ = "production"

import os, shutil, time
from math import log
from EDPluginControl import EDPluginControl
from EDFactoryPluginStatic import EDFactoryPluginStatic
from EDFactoryPlugin        import edFactoryPlugin
from EDThreading import Semaphore
edFactoryPlugin.loadModule("XSDataExecCommandLine")
EDFactoryPluginStatic.loadModule("XSDataEdnaSaxs")
EDFactoryPluginStatic.loadModule("XSDataBioSaxsv1_0")
EDFactoryPluginStatic.loadModule("XSDataWaitFilev1_0")
from XSDataCommon           import XSDataString, XSDataStatus, XSDataFile, XSDataTime, XSDataInteger
from XSDataBioSaxsv1_0      import XSDataInputBioSaxsSmartMergev1_0, XSDataResultBioSaxsSmartMergev1_0, \
                            XSDataInputBioSaxsISPyBv1_0
from XSDataEdnaSaxs         import XSDataInputDataver, XSDataInputAutoSub, XSDataInputDatop, XSDataInputSaxsAnalysis
from XSDataWaitFilev1_0     import XSDataInputWaitMultiFile
from XSDataExecCommandLine  import XSDataInputRsync
from freesas.cormap import gof
from numpy import loadtxt

def cmp(a, b):
    """
    Helper function to sort XSDataFile object
    """
    strA = a.path.value
    strB = a.path.value
    if strA > strB:
        return 1
    elif strA < strB:
        return -1
    else:
        return 0


class EDPluginBioSaxsSmartMergev1_7(EDPluginControl):
    """
    This plugin takes a set of input data files (1D SAXS) measure
    their differences (versus previous and versus first) and merge those which are equivalent
    v1.1 adds information from AutoSub
    
    Controlled plugins:
     - Execplugins/WaitMultifile
    
     - EdnaSaxs/Atsas/DatAver
     - EdnaSaxs/Atsas/AutoSub
     - EdnaSaxs/SaxsAnalysis
    """
    # dictAve = {} #key=?; value=path to average file
    lastBuffer = None
    lastSample = None
    dictFrames = {} #dictionary key = filename of average, value = list of files used for average
    __strControlledPluginDataver = "EDPluginExecDataverv2_0"
    #__strControlledPluginDatcmp = "EDPluginExecDatcmpv2_0"
    __strControlledPluginWaitFile = "EDPluginWaitMultiFile"
    __strControlledPluginAutoSub = "EDPluginAutoSubv1_1"
    __strControlledPluginSaxsAnalysis = "EDPluginControlSaxsAnalysisv1_0"
    __strControlledPluginSaxsISPyB = "EDPluginBioSaxsISPyBv1_0"
    __configured = False
    CONF_MINIMUM_CURVE_FILE_SIZE = "MinCurveFileSize"
    minimumCurveFileSize = 10000
    semaphore = Semaphore()
    cpRsync = "EDPluginExecRsync"

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)

        #self.__edPluginExecDatcmp = None
        self.__edPluginExecDataver = None
        self.__edPluginExecWaitFile = None
        self.__edPluginExecAutoSub = None
        self.__edPluginSaxsAnalysis = None
        self.__edPluginSaxsISPyB = None
        self.setXSDataInputClass(XSDataInputBioSaxsSmartMergev1_0)
        #self.__edPluginExecDatCmp = None
        self.lstInput = []
        self.curves = []
        self.forgetLastSample = False
        self.lstMerged = []
        self.lstDiscarded = []
        self.lstXsdInput = []
        self.absoluteFidelity = 0.01#None
        self.relativeFidelity = 0.01#None
        self.dictSimilarities = {}  # key: 2-tuple of images, similarities
        self.lstStrInput = []
        self.autoRg = None
        self.gnom = None
        self.volume = None
        self.rti = None
        self.strRadiationDamage = None
        self.strMergedFile = None
        self.lstSub = []
        self.strSubFile = None
        self.fConcentration = None
        self.xsDataResult = XSDataResultBioSaxsSmartMergev1_0()
        self.xsBestBuffer = None
        self.bestBufferType = ""
        self.bufferFrames = []
        self.xsScatterPlot = None
        self.xsGuinierPlot = None
        self.xsKratkyPlot = None
        self.xsDensityPlot = None
        self.xsdSubtractedCurve = None
        self.outdir = None #directory on rnice for analysis results to go to

    def configure(self):
        """
        Configures the plugin from the configuration file with the following parameters:
         - curve_file_size: minimum size of the file.
        """
        EDPluginControl.configure(self)
        if not self.__configured:
            with self.semaphore:
                if not self.__configured:
                    self.DEBUG("EDPluginBioSaxsSmartMergev1_7.configure")
                    min_size = self.config.get(self.CONF_MINIMUM_CURVE_FILE_SIZE)
                    if min_size is None:
                        strMessage = 'EDPluginBioSaxsSmartMergev1_7.configure: %s Configuration parameter missing: \
            %s, defaulting to "%s"' % (self.getBaseName(), self.CONF_MINIMUM_CURVE_FILE_SIZE, self.minimumCurveFileSize)
                        self.WARNING(strMessage)
                        self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                    else:
                        self.__class__.minimumCurveFileSize = float(min_size)
                    self.__class__.__configured = True

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.inputCurves, "Input curve list is empty")
        self.checkMandatoryParameters(self.dataInput.mergedCurve, "Output curve filename  is empty")


    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.preProcess")
        # Load the execution plugin
        if self.dataInput.absoluteFidelity is not None:
            self.absoluteFidelity = self.dataInput.absoluteFidelity.value
            if self.absoluteFidelity == 0.0:
                self.absoluteFidelity = None
        if self.dataInput.relativeFidelity is not None:
            self.relativeFidelity = self.dataInput.relativeFidelity.value
            if self.relativeFidelity == 0.0:
                self.relativeFidelity = None
        self.lstInput = self.dataInput.inputCurves
        self.lstStrInput = [i.path.value for i in self.lstInput]
        self.__edPluginExecWaitFile = self.loadPlugin(self.__strControlledPluginWaitFile)

        if self.dataInput.mergedCurve is not None:
            self.strMergedFile = self.dataInput.mergedCurve.path.value

        if self.dataInput.subtractedCurve is not None:
            self.strSubFile = self.dataInput.subtractedCurve.path.value
            dirname = os.path.dirname(self.strSubFile)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)

    def getFramesByFilename(self, filename):
        if filename is not None:
            if filename in self.__class__.dictFrames.keys():
                return self.__class__.dictFrames[filename]
        return {"averaged": [], "discarded": []}

    def getAveragedFrameByFilename(self, filename):
        frames = self.getFramesByFilename(filename)
        if frames is not None:
            return frames["averaged"]
        return []


    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.process")

        xsdwf = XSDataInputWaitMultiFile(timeOut=XSDataTime(30),
                                         expectedSize=XSDataInteger(self.minimumCurveFileSize),
                                         expectedFile=[XSDataFile(i.path) for i in self.lstInput])
        self.__edPluginExecWaitFile.setDataInput(xsdwf)
        self.__edPluginExecWaitFile.connectFAILURE(self.doFailureExecWait)
        self.__edPluginExecWaitFile.connectSUCCESS(self.doSuccessExecWait)
        self.__edPluginExecWaitFile.executeSynchronous()
        if self.isFailure():
            return
        if len(self.lstInput) == 1:
            inp = self.lstInput[0].path.value
            dst = self.dataInput.mergedCurve.path.value
            if not os.path.isdir(os.path.dirname(dst)):
                self.error("Output directory for %s does not exist" % dst)
                os.makedirs(os.path.dirname(dst))
            if not os.path.exists(inp):
                self.warning("Input %s does not (yet?) exist" % inp)
                time.sleep(1.0)
            shutil.copyfile(inp, dst)
            self.addExecutiveSummaryLine("Got only one frame ... nothing to merge !!!")
        else:
            self.lstMerged = []
            
            if (self.absoluteFidelity is not None) or (self.relativeFidelity is not None):
                if self.absoluteFidelity is not None :
                    for idx, oneFile in enumerate(self.lstInput[1:]):
                        self.DEBUG("Calculating similarity of 0 and %s" % idx)
                        self.datcmp(self.lstInput[0], oneFile)
#                         edPluginExecAbsoluteFidelity = self.loadPlugin(self.__strControlledPluginDatcmp)
#                         xsd = XSDataInputDatcmp(inputCurve=[self.lstInput[0], oneFile])
#                         edPluginExecAbsoluteFidelity.setDataInput(xsd)
#                         edPluginExecAbsoluteFidelity.connectFAILURE(self.doFailureExecDatcmp)
#                         edPluginExecAbsoluteFidelity.connectSUCCESS(self.doSuccessExecDatcmp)
#                         edPluginExecAbsoluteFidelity.execute()
                if (self.relativeFidelity is not None):
                    if (self.absoluteFidelity is  None):
                        self.DEBUG("Calculating similarity of 0 and 1")
                        self.datcmp(self.lstInput[0],self.lstInput[1])
#                         edPluginExecAbsoluteFidelity = self.loadPlugin(self.__strControlledPluginDatcmp)
#                         xsd = XSDataInputDatcmp(inputCurve=[self.lstInput[0], self.lstInput[1] ])
#                         edPluginExecAbsoluteFidelity.setDataInput(xsd)
#                         edPluginExecAbsoluteFidelity.connectFAILURE(self.doFailureExecDatcmp)
#                         edPluginExecAbsoluteFidelity.connectSUCCESS(self.doSuccessExecDatcmp)
#                         edPluginExecAbsoluteFidelity.execute()
                    if (len(self.lstInput) > 2):
                        for idx, oneFile in enumerate(self.lstInput[2:]):
                            self.DEBUG("Calculating similarity of %s and %s" % (idx, idx + 1))
                            self.datcmp(self.lstInput[idx + 1], oneFile)
#                             edPluginExecRelativeFidelity = self.loadPlugin(self.__strControlledPluginDatcmp)
#                             xsd = XSDataInputDatcmp(inputCurve=[self.lstInput[idx + 1], oneFile])
#                             edPluginExecRelativeFidelity.setDataInput(xsd)
#                             edPluginExecRelativeFidelity.connectFAILURE(self.doFailureExecDatcmp)
#                             edPluginExecRelativeFidelity.connectSUCCESS(self.doSuccessExecDatcmp)
#                             edPluginExecRelativeFidelity.execute()
            self.synchronizePlugins()

            for idx, oneFile in enumerate(self.lstInput):
                if idx == 0:
                    self.lstMerged.append(oneFile)
                elif (self.absoluteFidelity is not None) and (self.relativeFidelity is not None):
                    if (idx - 1, idx) not in self.dictSimilarities:
                        self.ERROR("dict missing %i,%i: \n" % (idx - 1, idx) + "\n".join([ "%s: %s" % (key, self.dictSimilarities[key]) for key in self.dictSimilarities]))
                        self.resynchronize()

                    if (0, idx) not in self.dictSimilarities:
                        self.ERROR("dict missing %i,%i: \n" % (0, idx) + "\n".join([ "%s: %s" % (key, self.dictSimilarities[key]) for key in self.dictSimilarities]))
                        self.resynchronize()

                    if (self.dictSimilarities[(0, idx)] >= self.absoluteFidelity) and (self.dictSimilarities[(idx - 1, idx)] >= self.relativeFidelity):
                        self.lstMerged.append(oneFile)

                elif (self.absoluteFidelity is not None) :
                    if (0, idx) not in self.dictSimilarities:
                        self.ERROR("dict missing %i,%i: \n" % (0, idx) + "\n".join([ "%s: %s" % (key, self.dictSimilarities[key]) for key in self.dictSimilarities]))
                        self.resynchronize()

                    if (self.dictSimilarities[(0, idx)] >= self.absoluteFidelity):
                        self.lstMerged.append(oneFile)
                   
                elif (self.relativeFidelity is not None) :
                    if (idx - 1, idx) not in self.dictSimilarities:
                        self.ERROR("dict missing %i,%i: \n" % (idx - 1, idx) + "\n".join([ "%s: %s" % (key, self.dictSimilarities[key]) for key in self.dictSimilarities]))
                        self.resynchronize()

                    if (self.dictSimilarities[(idx - 1, idx)] >= self.relativeFidelity):
                        self.lstMerged.append(oneFile)
                    
                else:
                    self.lstMerged.append(oneFile)
            self.lstMerged.sort(cmp)
            if len(self.lstMerged) != len(self.lstInput):
                self.strRadiationDamage = "Radiation damage detected, merged %i curves" % len(self.lstMerged)
                self.WARNING(self.strRadiationDamage)
                self.addExecutiveSummaryLine("WARNING: " + self.strRadiationDamage)
            self.addExecutiveSummaryLine("Merging files: " + " ".join([os.path.basename(i.path.value) for i in self.lstMerged]))
            if len(self.lstMerged) == 1:
                self.rewriteHeader(self.lstMerged[0].path.value, self.strMergedFile)
            else:
                self.__edPluginExecDataver = self.loadPlugin(self.__strControlledPluginDataver)
                xsd = XSDataInputDataver(inputCurve=self.lstMerged)
                # outputCurve=self.dataInput.mergedCurve,
                self.__edPluginExecDataver.setDataInput(xsd)
                self.__edPluginExecDataver.connectSUCCESS(self.doSuccessExecDataver)
                self.__edPluginExecDataver.connectFAILURE(self.doFailureExecDataver)
                self.__edPluginExecDataver.executeSynchronous()

        if (self.fConcentration == 0) and (self.strSubFile is not None):
            if (self.__class__.lastBuffer is not None) and (self.__class__.lastSample is not None):
                self.__edPluginExecAutoSub = self.loadPlugin(self.__strControlledPluginAutoSub)
                base = "_".join(os.path.basename(self.__class__.lastSample.path.value).split("_")[:-1])
                suff = os.path.basename(self.strSubFile).split("_")[-1]
                sub = os.path.join(os.path.dirname(self.strSubFile), base + "_" + suff)
                self.xsdSubtractedCurve = XSDataFile(XSDataString(sub))
                #self.curves.append(xsdSubtractedCurve)
                self.__edPluginExecAutoSub.dataInput = XSDataInputAutoSub(sampleCurve=self.__class__.lastSample,
                                         buffers=[self.__class__.lastBuffer, self.dataInput.mergedCurve],
                                         subtractedCurve=self.xsdSubtractedCurve)
                self.__edPluginExecAutoSub.connectSUCCESS(self.doSuccessExecAutoSub)
                self.__edPluginExecAutoSub.connectFAILURE(self.doFailureExecAutoSub)
                self.__edPluginExecAutoSub.executeSynchronous()

                if self.isFailure():
                    return

                self.__edPluginSaxsAnalysis = self.loadPlugin(self.__strControlledPluginSaxsAnalysis)
                self.__edPluginSaxsAnalysis.dataInput = XSDataInputSaxsAnalysis(scatterCurve=self.xsdSubtractedCurve,
                                                                                autoRg=self.autoRg,
                                                                                graphFormat=XSDataString("png"))
                self.__edPluginSaxsAnalysis.connectSUCCESS(self.doSuccessSaxsAnalysis)
                self.__edPluginSaxsAnalysis.connectFAILURE(self.doFailureSaxsAnalysis)
                self.__edPluginSaxsAnalysis.executeSynchronous()

            self.__class__.lastBuffer = self.dataInput.mergedCurve
            #self.__class__.lastSample = None #Information neededfor transfer to ISPyB
            self.forgetLastSample = True 
        else:
            self.__class__.lastSample = self.dataInput.mergedCurve

        if self.dataInput.sample and self.dataInput.sample.login and self.dataInput.sample.passwd and self.dataInput.sample.measurementID:
            self.addExecutiveSummaryLine("Registering to ISPyB")
            self.lstDiscarded = list(set(self.lstInput) - set(self.lstMerged))
            self.__class__.dictFrames[self.dataInput.mergedCurve] = {'averaged':self.lstMerged, 'discarded':self.lstDiscarded}
            self.__edPluginSaxsISPyB = self.loadPlugin(self.__strControlledPluginSaxsISPyB)
            if len(self.lstInput) > 1:
                frameAverage = XSDataInteger(len(self.lstInput))
                frameMerged = XSDataInteger(len(self.lstMerged))
            else:
                frameMerged = frameAverage = XSDataInteger(1)
            self.curves = [XSDataFile(i.path) for i in self.lstInput]
            self.discardedCurves = [XSDataFile(i.path) for i in self.lstDiscarded]
            self.mergedCurves = [XSDataFile(i.path) for i in self.lstMerged]

            averageFilePath = None
            if self.strMergedFile is not None:
                averageFilePath = XSDataFile(XSDataString(self.strMergedFile))

            

            self.sampleFrames = self.getAveragedFrameByFilename(self.__class__.lastSample)
            lastSample = None
            if self.__class__.lastSample is not None:
                lastSample = self.__class__.lastSample

            subtractedCurve = None	
            if self.xsdSubtractedCurve is not None:
                subtractedCurve = self.xsdSubtractedCurve
                xsdin = XSDataInputBioSaxsISPyBv1_0(sample=self.dataInput.sample,
                                                         autoRg=self.autoRg,
                                                         gnom=self.gnom,
                                                         volume=self.volume,
                                                         frameAverage=frameAverage,
                                                         frameMerged=frameMerged,
                                                         curves=self.curves,
            					                         discardedFrames=self.discardedCurves,
            					                         averagedFrames=self.mergedCurves,
                 					                     averageFilePath=averageFilePath,
            					                         bufferFrames = self.bufferFrames,
            					                         sampleFrames = self.sampleFrames,
                                                         bestBuffer=self.xsBestBuffer,
            					                         averageSample=lastSample,
                                                         scatterPlot=self.xsScatterPlot,
                                                         guinierPlot=self.xsGuinierPlot,
                                                         kratkyPlot=self.xsKratkyPlot ,
                                                         densityPlot=self.xsDensityPlot,
    						                             subtractedFilePath=subtractedCurve
    #                                                     destination=self.dataInput.sample.ispybDestination #duplicate, already in sample
                                                   )
                self.__edPluginSaxsISPyB.dataInput = xsdin
                self.__edPluginSaxsISPyB.connectSUCCESS(self.doSuccessISPyB)
                self.__edPluginSaxsISPyB.connectFAILURE(self.doFailureISPyB)
                self.__edPluginSaxsISPyB.execute()
            
        
       # transfer analysis data to correct location on nice
        if self.gnom is not None:
            self.outdir = os.path.join(os.path.dirname(os.path.dirname(self.lstStrInput[0])), "ednaAnalysis")
            basename = os.path.basename(os.path.splitext(self.gnom.gnomFile.path.value)[0])
            self.outdir = os.path.join(self.outdir, basename)
            if not os.path.isdir(self.outdir):
                os.makedirs(self.outdir)
             #self.outFile = os.path.join(outdir, "NoResults.html")
            workingdir = os.path.dirname(self.gnom.gnomFile.path.value)
     
            self.pluginRsync = self.loadPlugin(self.cpRsync)
            self.pluginRsync.dataInput = XSDataInputRsync(source=XSDataFile(XSDataString(workingdir)) ,
                                                                       destination=XSDataFile(XSDataString(self.outdir)),
                                                                       options=XSDataString("-avx"))
     
            self.pluginRsync.connectSUCCESS(self.doSuccessExecRsync)
            self.pluginRsync.connectFAILURE(self.doFailureExecRsync)
            self.pluginRsync.executeSynchronous()
    
            

        if self.forgetLastSample:
	    #Also redefine dictionary to contain the buffer just processed?
            self.__class__.lastSample = None

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.postProcess")
        # Create some output data
        self.xsDataResult.mergedCurve = self.dataInput.mergedCurve
        if self.strSubFile is not None and os.path.isfile(self.strSubFile):
            self.xsDataResult.subtractedCurve = XSDataFile(XSDataString(self.strSubFile))
        self.xsDataResult.autoRg = self.autoRg
        self.xsDataResult.gnom = self.gnom
        self.xsDataResult.volume = self.volume
        self.xsDataResult.rti = self.rti

    def finallyProcess(self, _edObject=None):
        EDPluginControl.finallyProcess(self)

        self.xsDataResult.status = XSDataStatus(message=self.getXSDataMessage(),
                                          executiveSummary=XSDataString(os.linesep.join(self.getListExecutiveSummaryLines())))
        self.dataOutput = self.xsDataResult

    def resynchronize(self):
        """
        Sometimes plugins are not started or not yet finished... 
        """
        time.sleep(1)  # this is just a way to give the focus to other threads
        self.synchronizePlugins()
        self.ERROR("I slept a second, waiting for sub-plugins to finish")


    def rewriteHeader(self, infile=None, output=None, hdr="#", linesep=os.linesep):
        """
        Write the output file with merged data with the good header.

# BSA sample
#Sample c= 0.0 mg/ml
#
# Merged from: file1
# Merged from: file2
# Merged from: file4
#
# Sample Information:
# Concentration: 0.0
# Code: BSA

        @param infile: path of the original data file where original data are taken from
        @param outfile: path of the destination file (by default, from XSD)
        @param hdr: header marker, often "#"
        @param linesep: line separator, usually "\n" or "\r\n" depending on the Operating System 
        @return: None
        """
        Code = Concentration = None
        frameMax = exposureTime = measurementTemperature = storageTemperature = None
        originalFile = self.lstMerged[0].path.value
        headers = [line.strip() for line in open(originalFile) if line.startswith("#")]
        Comments = headers[0][1:].strip()
        for line in headers:
            if "title =" in  line:
                Comments = line.split("=", 1)[1].strip()
            elif "Comments =" in line:
                Comments = line.split("=", 1)[1].strip()
            elif "Concentration:" in line:
                Concentration = line.split(":", 1)[1].strip()
            elif "Concentration =" in line:
                Concentration = line.split("=", 1)[1].strip()
            elif "Code =" in line:
                Code = line.split("=", 1)[1].strip()
            elif "Code:" in line:
                Code = line.split(":", 1)[1].strip()
            elif "Storage Temperature" in line:
                storageTemperature = line.split(":", 1)[1].strip()
            elif "Measurement Temperature" in line:
                measurementTemperature = line.split(":", 1)[1].strip()
            elif "Time per frame" in line:
                exposureTime = line.split("=", 1)[1].strip()
            elif "Frame" in line:
                frameMax = line.split()[-1]
        try:
            c = float(Concentration)
        except Exception:
            c = -1.0
        self.fConcentration = c
        lstHeader = ["%s %s" % (hdr, Comments), "%s Sample c= %s mg/ml" % (hdr, Concentration), hdr]
        if frameMax is not None:
            lstHeader.append(hdr + " Number of frames collected: %s" % frameMax)
        if exposureTime is not None:
            lstHeader.append(hdr + " Exposure time per frame: %s" % exposureTime)

        if self.strRadiationDamage is None:
            lstHeader.append("%s No significant radiation damage detected, merging %i files" % (hdr, len(self.lstMerged)))
        else:
            lstHeader.append("%s WARNING: %s" % (hdr, self.strRadiationDamage))
        lstHeader += [hdr + " Merged from: %s" % i.path.value for i in self.lstMerged]
        if self.lstSub:
            lstHeader.append(hdr)
            lstHeader += ["%s %s" % (hdr, i) for i in self.lstSub]
        lstHeader += [hdr, hdr + " Sample Information:"]
        if storageTemperature is not None:
            lstHeader.append(hdr + " Storage Temperature (degrees C): %s" % storageTemperature)
        if measurementTemperature is not None:
            lstHeader.append(hdr + " Measurement Temperature (degrees C): %s" % measurementTemperature)

        lstHeader += [hdr + " Concentration: %s" % Concentration, "# Code: %s" % Code]
        if infile is None:
            infile = self.strMergedFile
        if output is None:
            output = self.strMergedFile
        data = linesep.join([ i.strip() for i in open(infile) if not i.startswith(hdr)])
        with open(output, "w") as outfile:
            outfile.write(linesep.join(lstHeader))
            if not data.startswith(linesep):
                outfile.write(linesep)
            outfile.write(data)


    def doSuccessExecWait(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessExecWait")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_6.doSuccessExecWait")
        self.retrieveMessages(_edPlugin)
        xsdo = _edPlugin.dataOutput
        self.DEBUG("ExecWait Output:%s" % xsdo.marshal())
        if (xsdo.timedOut is not None) and  bool(xsdo.timedOut.value):
            strErr = "Error in waiting for all input files to arrive"
            self.addExecutiveSummaryLine("EDPluginBioSaxsSmartMergev1_6.doSuccessExecWait :" + strErr)
            self.ERROR(strErr)
            self.setFailure()


    def doFailureExecWait(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureExecWait")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_6.doFailureExecWait")
        self.retrieveMessages(_edPlugin)
        strErr = "Error in waiting for all input files to arrive"
        self.ERROR(strErr)
        self.addExecutiveSummaryLine(strErr)
        self.setFailure()


    def doSuccessExecDataver(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessExecDataver")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doSuccessExecDataver")
        self.retrieveMessages(_edPlugin)
        self.rewriteHeader(_edPlugin.dataOutput.outputCurve.path.value, output=self.strMergedFile)


    def doFailureExecDataver(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureExecDataver")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_6.doFailureExecDataver")
        self.retrieveMessages(_edPlugin)
        strErr = "Error in Processing of Atsas 'dataver'"
        self.addExecutiveSummaryLine(strErr)
        self.ERROR(strErr)
        self.setFailure()

    def datcmp(self, file1, file2):
    #def doSuccessExecDatcmp(self, _edPlugin=None):
#         self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessExecDatcmp")
#         self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doSuccessExecDatcmp")
#         self.retrieveMessages(_edPlugin)
        data1 = loadtxt(file1.path.value)
        data2 = loadtxt(file2.path.value)
        with self.locked():
            try: 
                cor = gof(data1,data2)
                fidelity = cor.P
            except Exception as err:
                strErr = "No fidelity in output of datcmp"
                self.error(strErr)
                self.addExecutiveSummaryLine(strErr)
                self.addExecutiveSummaryLine(str(err))
                #self.lstExecutiveSummary.append(strErr)
                #self.lstExecutiveSummary.append(err)
                if not "maximum recursion depth exceeded" in str(err): 
                    self.setFailure()
                fidelity = 0
            lstIdx = [self.lstStrInput.index(file1.path.value), self.lstStrInput.index(file2.path.value)]
            lstIdx.sort()
            self.dictSimilarities[tuple(lstIdx)] = fidelity
            lstIdx.reverse()
            self.dictSimilarities[tuple(lstIdx)] = fidelity
            if fidelity <= 0:
                logFid = "infinity"
            else:
                logFid = "%.2f" % (-log(fidelity, 10))
            self.addExecutiveSummaryLine("-log(Fidelity) between %s and %s is %s" % (os.path.basename(file1.path.value), os.path.basename(file2.path.value), logFid))


#     def doFailureExecDatcmp(self, _edPlugin=None):
#         self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureExecDatcmp")
#         self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_6.doFailureExecDatcmp")
#         self.retrieveMessages(_edPlugin)
#         self.addExecutiveSummaryLine("Failure in Processing of Atsas 'datcmp'")
#         self.setFailure()


    def doSuccessExecAutoSub(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessExecAutoSub")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doSuccessExecAutoSub")
        self.retrieveMessages(_edPlugin)
        self.autoRg = _edPlugin.dataOutput.autoRg
        if _edPlugin.dataOutput.subtractedCurve is not None:
            subcurve = _edPlugin.dataOutput.subtractedCurve
            if os.path.exists(subcurve.path.value):
                self.strSubFile = subcurve.path.value
        self.xsBestBuffer = _edPlugin.dataOutput.bestBuffer  
        self.bestBufferType = _edPlugin.dataOutput.bestBufferType.value
        if self.bestBufferType == 'average':
            self.bufferFrames = self.lstMerged
	    #if self.__class__.dictFrames[self.__class__.lastBuffer] is not None:
	    #	if self.__class__.dictFrames[self.__class__.lastBuffer]['averaged'] is not None:
	    #	    self.bufferFrames = self.bufferFrames + self.__class__.dictFrames[self.__class__.lastBuffer]['averaged'] 
            self.bufferFrames = self.bufferFrames + self.getAveragedFrameByFilename(self.__class__.lastBuffer)

        elif self.bestBufferType == self.__class__.lastBuffer:
            #if self.__class__.dictFrames[self.__class__.lastBuffer] is not None:
	    #    if self.__class__.dictFrames[self.__class__.lastBuffer]['averaged'] is not None:
	    #        self.bufferFrames = self.__class__.dictFrames[self.__class__.lastBuffer]['averaged']
            self.bufferFrames = self.getAveragedFrameByFilename(self.__class__.lastBuffer)
        else:
            self.bufferFrames = self.lstMerged
        self.addExecutiveSummaryLine(_edPlugin.dataOutput.status.executiveSummary.value)


    def doFailureExecAutoSub(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureExecAutoSub")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doFailureExecAutoSub")
        self.retrieveMessages(_edPlugin)
        strErr = "Error in Processing of EDNA 'AutoSub'"
        # clean up "EDNA memory"
        self.__class__.lastBuffer = self.dataInput.mergedCurve
        self.__class__.lastSample = None
        self.ERROR(strErr)
        self.addExecutiveSummaryLine(strErr)
        self.setFailure()

    def doSuccessSaxsAnalysis(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessSaxsAnalysis")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doSuccessSaxsAnalysis")
        self.retrieveMessages(_edPlugin)
        self.gnom = _edPlugin.dataOutput.gnom
        self.volume = _edPlugin.dataOutput.volume
        self.rti = _edPlugin.dataOutput.rti
        self.xsScatterPlot = _edPlugin.dataOutput.scatterPlot
        self.xsGuinierPlot = _edPlugin.dataOutput.guinierPlot
        self.xsKratkyPlot = _edPlugin.dataOutput.kratkyPlot
        if hasattr(_edPlugin.dataOutput, 'kratkyRgPlot'):
            self.xsKratkyRgPlot = _edPlugin.dataOutput.kratkyRgPlot
        if hasattr(_edPlugin.dataOutput, 'kratkyVcPlot'):
            self.xsKratkyVcPlot = _edPlugin.dataOutput.kratkyVcPlot

        self.xsDensityPlot = _edPlugin.dataOutput.densityPlot
        self.addExecutiveSummaryLine(_edPlugin.dataOutput.status.executiveSummary.value)

    def doFailureSaxsAnalysis(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureSaxsAnalysis")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_6.doFailureSaxsAnalysis")
        self.retrieveMessages(_edPlugin)
        strErr = "Error in Processing of EDNA SaxsAnalysis = AutoRg => datGnom => datPorod"
        self.ERROR(strErr)
        self.addExecutiveSummaryLine(strErr)
        self.setFailure()

    def doSuccessISPyB(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessISPyB")
        self.addExecutiveSummaryLine("Registered in ISPyB")
        self.retrieveMessages(_edPlugin)

    def doFailureISPyB(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureISPyB")
        self.addExecutiveSummaryLine("Failed to registered in ISPyB")
        self.retrieveMessages(_edPlugin)


    def doSuccessExecRsync(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doSuccessExecRsync")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doSuccessExecRsync")
        self.retrieveMessages(_edPlugin)
        self.gnom.gnomFile = XSDataFile(XSDataString(os.path.join(self.outdir,os.path.split(self.gnom.gnomFile.path.value)[1])))


    def doFailureExecRsync(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsSmartMergev1_7.doFailureExecRsync")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsSmartMergev1_7.doFailureExecRsync")
        self.retrieveMessages(_edPlugin)
        self.setFailure()
