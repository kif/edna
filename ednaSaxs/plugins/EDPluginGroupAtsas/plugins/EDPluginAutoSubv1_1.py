# coding: utf-8
#
#    Project: EdnaSaxs/Atsas
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
__copyright__ = "2011 ESRF"
__status__ = "production"
__date__ = "28/08/2015"

import os
import shutil
from EDPluginControl import EDPluginControl
from XSDataCommon import XSDataFile, XSDataString, XSDataStatus
from XSDataEdnaSaxs import XSDataInputAutoSub, XSDataInputDataver, \
    XSDataInputAutoRg, XSDataInputDatop, XSDataAutoRg, \
    XSDataLength, XSDataBoolean, XSDataInteger, XSDataDouble
from XSDataEdnaSaxs import XSDataResultAutoSub
from EDThreading import Semaphore
_sem = Semaphore()
import numpy
from freesas.cormap import gof


def copy(src, dst):
    """
    Operating system independent copy that does symbolic links on posix system.
    
    This implements a lock to be extremely thread-safe as BsxCube call many times the same plugin ... at the same time 
    """
    with _sem:
        if os.path.exists(src):
            if os.name == "posix":
                if os.path.islink(dst):
                    os.unlink(dst)
                os.symlink(src, dst)
            else:
                if os.path.exists(dst):
                    os.unlink(dst)
                shutil.copyfile(src, dst)


class EDPluginAutoSubv1_1(EDPluginControl):
    """
    This is an EDNA plugin that tries to mimic the behavour of autosub:
    - find the best buffer (usually the merged one)
    - subtract data with the best buffer  
    """
    BUFFER_SIMILARITY = 1e-300

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputAutoSub)
        self.__strPluginDatop = "EDPluginExecDatopv2_0"
        #self.__strPluginDatcmp = "EDPluginExecDatcmpv3_0"
        self.__strPluginDataver = "EDPluginExecDataverv2_0"
        self.__strPluginAutoRg = "EDPluginExecAutoRgv1_1"
        self.__edPluginDatcmp = None
        self.__edPluginDatop = None
        self.__edPluginDataver = None
        self.__edPluginAutoRg = None
        self.xsDataResult = XSDataResultAutoSub()

        self.buffers = []
        self.sampleCurve = None
        self.subtractedCurve = None
        self.outdir = None
        self.bestBuffer = None
        self.averBuffer = None
        self.tempSub = [] #Subtraction for finding best buffer
        self.actualBestBuffer = None
        self.bestBufferType = "" #"average" or the file name of the best buffer
        self.fidelity = None
        self.lstProcessLog = []
        self.dictRg = {} #key: filename, value = (Rg,I0)
        self.fConcentration = -1
        self.headers = {}
        self.bufferData = {}

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginAutoSubv1_1.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.sampleCurve, "No sample curve")
        self.checkMandatoryParameters(self.dataInput.buffers, "No buffer curves")

    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginAutoSubv1_1.preProcess")
        # Load the execution plugin
        self.__edPluginDataver = self.loadPlugin(self.__strPluginDataver)

        if self.dataInput.subtractedCurve is not None:
            self.subtractedCurve = self.dataInput.subtractedCurve.path.value
        self.buffers = [i.path.value for i in self.dataInput.buffers]
        self.sampleCurve = self.dataInput.sampleCurve.path.value
        if self.subtractedCurve is None:
            self.outdir = os.path.dirname(self.sampleCurve)
        else:
            self.outdir = os.path.dirname(self.subtractedCurve)

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginAutoSubv1_1.process")

        basename = os.path.splitext(os.path.basename(self.sampleCurve))[0]
        self.bestBuffer = os.path.join(self.outdir, "%s_bestbuffer.dat" % basename)
        self.averBuffer = os.path.join(self.outdir, "%s_averbuffer.dat" % basename)
        self.data = numpy.loadtxt(self.sampleCurve)
        for idx, name in enumerate(self.buffers):
            copy(name, os.path.join(self.outdir, "%s_buf%i.dat" % (basename, idx)))
        if len(self.buffers) == 1:
            self.actualBestBuffer = self.buffers[0]
            self.bestBufferType = str(self.actualBestBuffer)
            self.bufferIntensity = numpy.loadtxt(self.buffers[0])
        else:
#             self.__edPluginDataver.dataInput = XSDataInputDataver(inputCurve=self.dataInput.buffers,
#                                                                   outputCurve=XSDataFile(XSDataString(self.averBuffer)))
#             self.__edPluginDataver.connectSUCCESS(self.doSuccessExecDataver)
#             self.__edPluginDataver.connectFAILURE(self.doFailureExecDataver)
#             self.__edPluginDataver.execute()

            if len(self.buffers) == 2:
                self.bufferData[self.buffers[0]] = numpy.loadtxt(self.buffers[0])
                self.bufferData[self.buffers[1]] = numpy.loadtxt(self.buffers[1])
                
                #Average
                AverageIntensity = (self.bufferData[self.buffers[0]][:,1] + self.bufferData[self.buffers[1]][:,1])/2
                StdErr = numpy.sqrt((self.bufferData[self.buffers[0]][:,2]  * self.bufferData[self.buffers[0]][:,2]+ \
                        self.bufferData[self.buffers[1]][:,2]  * self.bufferData[self.buffers[1]][:,2])/2)
                self.averageBuffer = numpy.vstack((self.bufferData[self.buffers[0]][:,0], AverageIntensity, StdErr)).T
                with open(self.averBuffer, "w") as outfile:
                    numpy.savetxt(outfile,self.averageBuffer)
                
                self.datcmp(self.bufferData[self.buffers[0]],self.bufferData[self.buffers[1]])
#                 edPluginDatcmp = self.loadPlugin(self.__strPluginDatcmp)
#                 edPluginDatcmp.dataInput = XSDataInputDatcmp(inputCurve=self.dataInput.buffers)
#                 edPluginDatcmp.connectSUCCESS(self.doSuccessExecDatcmp)
#                 edPluginDatcmp.connectFAILURE(self.doFailureExecDatcmp)
#                 edPluginDatcmp.executeSynchronous()
                if self.isFailure() or (self.fidelity is None):
                    return
                if self.fidelity < self.BUFFER_SIMILARITY:
                    print "Different buffers"
                    # buffer are not the same: keeping the one with lowest Rg/I0 for the subtraction
                    #data = numpy.loadtxt(self.sampleCurve.path.value)
                    for buff in self.buffers: 
                        tempName = os.path.splitext(buff)[0] + '.sub'
                        self.doSubtract(self.data,self.bufferData[buff],tempName)
#                         Isub = data [:,1] - self.bufferData[buff][:,1]
#                         StdErr = numpy.sqrt(data[:,2] * data[:,2] + \
#                         self.bufferData[buff][:,2] * self.bufferData[buff][:,2])
# 
#                         with open(tempName, "w") as outfile:
#                             numpy.savetxt(outfile, numpy.vstack((data[:,0], Isub, StdErr)).T)

                       
#                         edPluginDatop = self.loadPlugin(self.__strPluginDatop)
#                         edPluginDatop.dataInput = XSDataInputDatop(operation=XSDataString("SUB"),
#                                    outputCurve=XSDataFile(XSDataString(os.path.splitext(buff)[0] + '.sub')),
#                                    inputCurve=[XSDataFile(XSDataString(self.sampleCurve)), XSDataFile(XSDataString(buff))])
#                         edPluginDatop.connectSUCCESS(self.doSuccessExecDatop)
#                         edPluginDatop.connectFAILURE(self.doFailureExecDatop)
#                         edPluginDatop.executeSynchronous()
                        self.tempSub += [XSDataFile(XSDataString(tempName))]
                    
                    edpluginRg = self.loadPlugin(self.__strPluginAutoRg)
                    edpluginRg.dataInput = XSDataInputAutoRg(inputCurve=self.tempSub)
                    edpluginRg.connectSUCCESS(self.doSuccessExecAutoRg)
                    edpluginRg.connectFAILURE(self.doFailureExecAutoRg)
                    edpluginRg.executeSynchronous()
                    actualBestTmp = self.dictRg.keys()[self.dictRg.values().index(min(self.dictRg.values()))]
                    self.actualBestBuffer = os.path.splitext(actualBestTmp)[0] + '.dat'
                    self.bestBufferType = str(self.actualBestBuffer)
                    self.bufferIntensity = self.bufferData[self.actualBestBuffer]
                else:
                    print "similar buffers"
                    self.actualBestBuffer = self.averBuffer
                    self.bestBufferType = "average"
                    self.bufferIntensity = self.averageBuffer
                print self.actualBestBuffer
            else:
                self.synchronizePlugins()
                strError = "You should specify exactly 2 buffers for guessing best buffer, I got: " + ", ".join(self.buffers)
                self.WARNING(strError)
                self.lstProcessLog.append(strError)
                self.actualBestBuffer = self.averBuffer
                self.bufferIntensity = self.averageBuffer

        copy(self.actualBestBuffer, self.bestBuffer)
        print "copied to ", self.bestBuffer
        self.lstProcessLog.append("Best buffer is %s" % self.actualBestBuffer)
        if self.isFailure() or not os.path.exists(self.bestBuffer) or not os.path.exists(self.sampleCurve):
            return
#         edPluginDatop = self.loadPlugin(self.__strPluginDatop)
#         edPluginDatop.dataInput = XSDataInputDatop(operation=XSDataString("SUB"),
#                                outputCurve=XSDataFile(XSDataString(self.subtractedCurve)),
#                                inputCurve=[XSDataFile(XSDataString(self.sampleCurve)), XSDataFile(XSDataString(self.bestBuffer))])
#         edPluginDatop.connectSUCCESS(self.doSuccessExecDatop)
#         edPluginDatop.connectFAILURE(self.doFailureExecDatop)
#         edPluginDatop.executeSynchronous()
        self.doSubtract(self.data, self.bufferIntensity, self.subtractedCurve)

        if self.isFailure() or not os.path.exists(self.subtractedCurve):
            return
        self.headers = self.parseHeaders(self.sampleCurve)
        edpluginRg = self.loadPlugin(self.__strPluginAutoRg)
        edpluginRg.dataInput = XSDataInputAutoRg(inputCurve=[XSDataFile(XSDataString(self.subtractedCurve))])
        edpluginRg.connectSUCCESS(self.doSuccessExecAutoRg)
        edpluginRg.connectFAILURE(self.doFailureExecAutoRg)
        edpluginRg.executeSynchronous()

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginAutoSubv1_1.postProcess")
        # Create some output data

        self.xsDataResult.status = XSDataStatus(executiveSummary=XSDataString(os.linesep.join(self.lstProcessLog)))
        self.xsDataResult.subtractedCurve = XSDataFile(XSDataString(self.subtractedCurve))
        self.xsDataResult.bestBufferType = XSDataString(self.bestBufferType)
        if self.actualBestBuffer:
            self.xsDataResult.bestBuffer = XSDataFile(XSDataString(self.actualBestBuffer))
        self.dataOutput = self.xsDataResult

#     def doSuccessExecDatop(self, _edPlugin=None):
#         self.DEBUG("EDPluginAutoSubv1_1.doSuccessExecDatop")
#         self.retrieveSuccessMessages(_edPlugin, "EDPluginAutoSubv1_1.doSuccessExecDatop")
# 
#     def doFailureExecDatop(self, _edPlugin=None):
#         self.DEBUG("EDPluginAutoSubv1_1.doFailureExecDatop")
#         self.retrieveFailureMessages(_edPlugin, "EDPluginAutoSubv1_1.doFailureExecDatop")
#         self.lstProcessLog.append("Failure in datop (subtraction of the buffer")
#         self.setFailure()

    def doSuccessExecDataver(self, _edPlugin=None):
        self.DEBUG("EDPluginAutoSubv1_1.doSuccessExecDataver")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginAutoSubv1_1.doSuccessExecDataver")

    def doFailureExecDataver(self, _edPlugin=None):
        self.DEBUG("EDPluginAutoSubv1_1.doFailureExecDataver")
        self.retrieveFailureMessages(_edPlugin, "EDPluginAutoSubv1_1.doFailureExecDataver")
        self.lstProcessLog.append("Failure in dataver (averaging of buffer")
        self.setFailure()

    def doSuccessExecAutoRg(self, _edPlugin=None):
        self.DEBUG("EDPluginAutoSubv1_1.doSuccessExecAutoRg")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginAutoSubv1_1.doSuccessExecAutoRg")
        for data in _edPlugin.dataOutput.autoRgOut:
            self.dictRg[data.filename.path.value] = (data.rg.value, data.i0.value)

        if (len(_edPlugin.dataOutput.autoRgOut) == 1) and os.path.exists(self.subtractedCurve):  # rewrite Headers
            res = _edPlugin.dataOutput.autoRgOut[0]
            #Scale I0 by concentration
            if self.fConcentration > 1e-6:
                res.i0.value = res.i0.value / self.fConcentration
                res.i0Stdev.value = res.i0Stdev.value / self.fConcentration
            lstRg = []
            if res.rg.value < 1e-6:
                lstRg.append("AutoRg: Rg   =   %.4f +/- %.2f" % (res.rg.value, res.rgStdev.value))
            else:
                lstRg.append("AutoRg: Rg   =   %.4f +/- %.2f ( %.1f%%)" % (res.rg.value, res.rgStdev.value, 100. * res.rgStdev.value / res.rg.value))
            lstRg.append("AutoRg: I(0) =   %.2f +/- %.4f" % (res.i0.value, res.i0Stdev.value))
            lstRg.append("AutoRg: Points   %i to %i ( %i total)" % (res.firstPointUsed.value, res.lastPointUsed.value , 1 + res.lastPointUsed.value - res.firstPointUsed.value))
            lstRg.append("AutoRg: Quality: %.1f%%" % (res.quality.value * 100.0))

            self.rewriteHeader(output=self.subtractedCurve,
                               infile=self.sampleCurve,
                               extraHeaders=["Buffer subtracted data file:",
                                             "Protein DatFile = %s" % self.sampleCurve,
                                             "Buffer DatFile = %s" % self.actualBestBuffer,
                                             ""] + lstRg,
                               scale=True)
            self.lstProcessLog += lstRg
            self.xsDataResult.autoRg = res

    def doFailureExecAutoRg(self, _edPlugin=None):
        self.DEBUG("EDPluginAutoSubv1_1.doFailureExecAutoRg")
        self.retrieveFailureMessages(_edPlugin, "EDPluginAutoSubv1_1.doFailureExecAutoRg")
        self.lstProcessLog.append("Failure in AutoRg")
        if len(_edPlugin.dataInput.inputCurve) == 2:
            # we were comparing 2 buffers but were analyzed by AutoRg
            for fn in self.buffers:
                self.dictRg[fn] = (0, numpy.loadtxt(fn, unpack=True)[1].sum())
        elif os.path.exists(self.subtractedCurve):
            #we don't want to fail the subtraction plugin because the result has no Rg - default to 0
            res = XSDataAutoRg()
            res.rg = XSDataLength(0.0)
            res.rgStdev = XSDataLength(0.0)
            res.i0 = XSDataDouble(0.0)
            res.i0Stdev = XSDataDouble(0.0)
            res.firstPointUsed = XSDataInteger(0)
            res.lastPointUsed = XSDataInteger(0)
            res.quality = XSDataDouble(0.0)
            res.isagregated = XSDataDouble(bool(0))
            res.filename = XSDataFile(XSDataString(self.subtractedCurve))
            self.xsDataResult.autoRg = res
        else:
            self.setFailure()

    def datcmp(self,data1,data2):
        cor = gof(data1,data2)
        self.fidelity = cor.P
#     def doSuccessExecDatcmp(self, _edPlugin=None):
#         self.DEBUG("EDPluginAutoSubv1_1.doSuccessExecDatcmp")
#         self.retrieveSuccessMessages(_edPlugin, "EDPluginAutoSubv1_1.doSuccessExecDatcmp")
#         self.fidelity = _edPlugin.dataOutput.fidelity.value
# 
#     def doFailureExecDatcmp(self, _edPlugin=None):
#         self.DEBUG("EDPluginAutoSubv1_1.doFailureExecDatcmp")
#         self.retrieveFailureMessages(_edPlugin, "EDPluginAutoSubv1_1.doFailureExecDatcmp")
#         self.lstProcessLog.append("Failure in datcmp (comparison of 2 buffers)")
#         self.setFailure()

    def doSubtract(self, data, aBuffer, filename):
        tempName = filename
        Isub = data [:,1] - aBuffer[:,1]
        StdErr = numpy.sqrt(data[:,2] * data[:,2] + \
                aBuffer[:,2] * aBuffer[:,2])

        with open(tempName, "w") as outfile:
            numpy.savetxt(outfile, numpy.vstack((data[:,0], Isub, StdErr)).T)
        

    def parseHeaders(self, infile=None, hdr="#"):
        """
        Read the headers 
        @param infile: filename of the data file 
        @param hdr: header marker: # 
        @return: headers as a dictionary 
        """
        headers = {}
#         Code = Concentration = None
#         frameMax = exposureTime = measurementTemperature = storageTemperature = None

        headLines = [line.strip() for line in open(infile) if line.startswith(hdr)]
        headers["Comments"] = headLines[0][1:].strip()
        for line in headLines:
            if "title =" in line:
                headers["Comments"] = line.split("=", 1)[1].strip()
            elif "Comments =" in line:
                headers["Comments"] = line.split("=", 1)[1].strip()
            elif "Concentration:" in line:
                headers["Concentration"] = line.split(":", 1)[1].strip()
            elif "Concentration =" in line:
                headers["Concentration"] = line.split("=", 1)[1].strip()
            elif "Code =" in line:
                headers["Code"] = line.split("=", 1)[1].strip()
            elif "Code:" in line:
                headers["Code"] = line.split(":", 1)[1].strip()
            elif "Storage Temperature" in line:
                headers["storageTemperature"] = line.split(":", 1)[1].strip()
            elif "Measurement Temperature" in line:
                headers["measurementTemperature"] = line.split(":", 1)[1].strip()
            elif "Time per frame" in line:
                headers["exposureTime"] = line.split("=", 1)[1].strip()
#            elif "Frame" in line:
#                frameMax = line.split()[-1]
        try:
            self.fConcentration = float(headers["Concentration"])
        except Exception:
            self.fConcentration = -1
        return headers

    def rewriteHeader(self, infile=None, output=None, hdr="#", linesep=os.linesep, extraHeaders=None, scale=False):
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
        @param outfile: path of the destination file
        @param hdr: header marker, often "#"
        @param linesep: line separator, usually "\n" or "\r\n" depending on the Operating System
        @param scale: if True, divide by concentration   
        @return: None
        """
        if extraHeaders is None:
            extraHeaders = []

        lstHeader = [self.headers.get("Comments", "No Comment!"), "Sample c= %s mg/ml" % self.headers.get("Concentration", -1)]
        if scale and self.fConcentration != 0:
            lstHeader.append("Intensity scaled by concentration")
            df = [i.split() for i in open(output) if not i.startswith(hdr)]
            data = linesep.join(["%14s %14.6e %14.6e" % (w[0], float(w[1]) / self.fConcentration, float(w[2]) / self.fConcentration)
                                 for w in df if len(w) >= 3])
        else:
            data = linesep.join([i.strip() for i in open(output) if not i.startswith(hdr)])

        lstHeader += extraHeaders
        lstHeader.append("")
        if "exposureTime" in self.headers:
            lstHeader.append("Exposure time per frame: %s" % self.headers["exposureTime"])
        lstHeader += ["", "Sample Information:"]
        if "storageTemperature" in self.headers:
            lstHeader.append("Storage Temperature (degrees C): %s" % self.headers["storageTemperature"])
        if "measurementTemperature" in self.headers:
            lstHeader.append("Measurement Temperature (degrees C): %s" % self.headers["measurementTemperature"])
        lstHeader += ["Concentration: %s" % self.headers.get("Concentration", "-1"), "Code: %s" % self.headers.get("Code", "None")]
        with open(output, "w") as outfile:
            outfile.write(linesep.join([hdr + " " + i for i in lstHeader]))
            outfile.write(linesep)
            outfile.write(data)
