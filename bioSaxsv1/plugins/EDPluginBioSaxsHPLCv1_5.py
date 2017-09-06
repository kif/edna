# coding: utf-8
#
#    Project:BioSaxs
#             http://www.edna-site.org
#
#    Copyright (C) 2014 ESRF
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
__copyright__ = "2014 ESRF"
__date__ = "16/09/2016"
__status__ = "development"

import os
import fabio
import numpy
from EDPluginControl import EDPluginControl
from EDThreading import Semaphore
from EDUtilsArray import EDUtilsArray

from EDFactoryPlugin import edFactoryPlugin
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
edFactoryPlugin.loadModule("XSDataEdnaSaxs")
from XSDataBioSaxsv1_0 import XSDataInputBioSaxsHPLCv1_0, XSDataResultBioSaxsHPLCv1_0, \
                              XSDataInputBioSaxsProcessOneFilev1_0, XSDataRamboTainer
from XSDataEdnaSaxs import XSDataInputDatcmp, XSDataInputAutoRg
from XSDataCommon import XSDataFile, XSDataString, XSDataStatus, XSDataTime, \
    XSDataDouble

from EDUtilsBioSaxs import HPLCframe, HPLCrun, RamboTainerInvariant

from freesas.cormap import gof
from freesas.autorg import autoRg, InsufficientDataError


class EDPluginBioSaxsHPLCv1_5(EDPluginControl):
    """
    plugin for processing Saxs data coming from HPLC

    runs subsequently:
    *ProcessOneFile,
    *subtraction of buffer
    *AutoRg

    todo:
    only store references: Wait for flush to construct HDF5 file and (possibly) web pages with PNG graphs
    
    Changes since v1.0:
    * do not perfom Gnom not Porod analysis
    * return Summed Intensity and time stamp
    * fix issue with timestamp.
    Changes since v1.1:
    * Read configuration for a buffer threshold 
    Changes since v1.2:
    lock during buffer curve creation.
    """

    strControlledPluginProcessOneFile = "EDPluginBioSaxsProcessOneFilev1_5"
    strControlledPluginDatop = "EDPluginExecDatopv2_0"
    strControlledPluginAutoRg = "EDPluginExecAutoRgv1_1"
    #strControlledPluginDatCmp = "EDPluginExecDatcmpv2_0"
    strControlledPluginDatAver = "EDPluginExecDataverv2_0"
    dictHPLC = {}  # key=runId, value= HPLCrun instance
    _sem = Semaphore()
    buffer_sem = Semaphore()
    SIMILARITY_THRESHOLD_SAMPLE_KEY = "similaritySample"
    SIMILARITY_THRESHOLD_BUFFER_KEY = "similarityBuffer"
    SIMILARITY_THRESHOLD_SAMPLE_DEFAULT = 0.01
    SIMILARITY_THRESHOLD_BUFFER_DEFAULT = 0.1
    SIMILARITY_THRESHOLD_SAMPLE = None
    SIMILARITY_THRESHOLD_BUFFER = None

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsHPLCv1_0)
        self.edPluginProcessOneFile = None
        self.edPluginSubtract = None
        self.edPluginAutoRg = None
        #self.edPluginDatCmp = None
        self.xsDataResult = XSDataResultBioSaxsHPLCv1_0()
        self.runId = None
        self.frameId = None
        self.frame = None
        self.hplc_run = None
        self.curve = None
        self.subtracted = None
        self.lstExecutiveSummary = []
        self.isBuffer = False
        self.intensity = None
        self.stdError = None
        self.ATSASRg = False

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.rawImage, "No raw image")
        self.checkMandatoryParameters(self.dataInput.sample, "no Sample parameter")
        self.checkMandatoryParameters(self.dataInput.experimentSetup, "No experimental setup parameter")

    def configure(self):
        """
        Configures the HPLC plugin by reading from the configuration file 
         - The threshold for similarity with sample & buffer
        """
        EDPluginControl.configure(self)
        if self.__class__.SIMILARITY_THRESHOLD_SAMPLE is None:
            with self._sem:
                if self.__class__.SIMILARITY_THRESHOLD_SAMPLE is None:
                    self.DEBUG("EDPluginBioSaxsHPLCv1_5.configure")
                    self.__class__.SIMILARITY_THRESHOLD_BUFFER = float(self.config.get(self.SIMILARITY_THRESHOLD_BUFFER_KEY, self.SIMILARITY_THRESHOLD_BUFFER_DEFAULT))
                    self.__class__.SIMILARITY_THRESHOLD_SAMPLE = float(self.config.get(self.SIMILARITY_THRESHOLD_SAMPLE_KEY, self.SIMILARITY_THRESHOLD_SAMPLE_DEFAULT))

    def checkRun(self):
        if self.hplc_run.deleted:
            str_err = "Processing of HPLC run already finished. Aborted processing of frame"
            self.ERROR(str_err)
            self.setFailure()
            raise RuntimeError(str_err)

    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.preProcess")
        sdi = self.dataInput
        if sdi.runId is not None:
            self.runId = sdi.runId.value
        else:
            path = sdi.rawImage.path.value
            if "_" in path:
                self.runId = path[::-1].split("_", 1)[1][::-1]
            else:
                self.runId = path
        with self._sem:
            if self.runId not in self.dictHPLC:
                self.dictHPLC[self.runId] = HPLCrun(self.runId)
        self.hplc_run = self.dictHPLC[self.runId]
        self.checkRun()  # Only continue if run is not deleted!

        if sdi.frameId is not None:
            self.frameId = sdi.frameId.value
        else:
            path = sdi.rawImage.path.value
            if "_" in path:
                digits = os.path.splitext(os.path.basename(path))[0].split("_")[-1]
                try:
                    self.frameId = int(digits)
                except ValueError:
                    self.WARNING("frameId is supposed to be an integer, I got %s" % digits)
                    self.frameId = digits
            else:
                self.warning("using frameID=0 in tests, only")
                self.frameId = 0
        with self._sem:
            self.frame = HPLCframe(self.runId, self.frameId)
            self.hplc_run.frames[self.frameId] = self.frame

        if sdi.bufferCurve and os.path.exists(sdi.bufferCurve.path.value):
            with self._sem:
                self.hplc_run.buffer = sdi.bufferCurve.path.value

        if self.hplc_run.hdf5_filename:
            hplc = self.hplc_run.hdf5_filename
        elif sdi.hplcFile:
            hplc = sdi.hplcFile.path.value
        else:
            path = sdi.rawImage.path.value
            dirname = os.path.dirname(os.path.dirname(path))
            hplc_dir = os.path.join(dirname, "HPLC")
            if not os.path.exists(hplc_dir):
                try:
                    os.mkdir(hplc_dir)
                except:
                    pass  # I don't care

            if "_" in path:
                hplc = "_".join(os.path.splitext(os.path.basename(path))[0].split("_")[:-1]) + ".h5"
            else:
                hplc = os.path.splitext(os.path.basename(path))[0] + ".h5"
            hplc = os.path.join(hplc_dir, hplc)
        if not self.hplc_run.hdf5_filename:
            with self._sem:
                self.hplc_run.init_hdf5(hplc)
#        self.xsDataResult.hplcFile = XSDataFile(XSDataString(hplc))

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.process")

        xsdIn = XSDataInputBioSaxsProcessOneFilev1_0(rawImage=self.dataInput.rawImage,
                                                     sample=self.dataInput.sample,
                                                     experimentSetup=self.dataInput.experimentSetup,
                                                     rawImageSize=self.dataInput.rawImageSize,
                                                     normalizedImage=self.dataInput.normalizedImage,
                                                     integratedCurve=self.dataInput.integratedCurve,
                                                     runId=self.dataInput.runId,
                                                     frameId=self.dataInput.frameId)
        self.edPluginProcessOneFile = self.loadPlugin(self.strControlledPluginProcessOneFile)
        self.edPluginProcessOneFile.dataInput = xsdIn
        self.edPluginProcessOneFile.connectSUCCESS(self.doSuccessProcessOneFile)
        self.edPluginProcessOneFile.connectFAILURE(self.doFailureProcessOneFile)
        self.edPluginProcessOneFile.executeSynchronous()

        if self.isFailure():
            return
        if self.hplc_run.first_curve == self.curve:
            return
#         xsdIn = XSDataInputDatcmp(inputCurve=[XSDataFile(XSDataString(self.hplc_run.first_curve)),
#                                               XSDataFile(XSDataString(self.curve))])
#         self.edPluginDatCmp = self.loadPlugin(self.strControlledPluginDatCmp)
#         self.edPluginDatCmp.dataInput = xsdIn
#         self.edPluginDatCmp.connectSUCCESS(self.doSuccessDatCmp)
#         self.edPluginDatCmp.connectFAILURE(self.doFailureDatCmp)
#         self.edPluginDatCmp.executeSynchronous()

        self.datcmp(self.hplc_run.firstCurveIntensity, self.intensity)

        if self.isFailure() or self.isBuffer:

            return


        #######################
        # # DatOp subtraction ##
        #######################

        if self.dataInput.subtractedCurve is not None:
            self.subtracted = self.dataInput.subtractedCurve.path.value
        else:
            self.subtracted = os.path.splitext(self.curve)[0] + "_sub.dat"

        Isub = self.intensity - self.hplc_run.buffer_I
        StdErr = numpy.sqrt(self.stdError * self.stdError + \
                            self.hplc_run.buffer_Stdev * self.hplc_run.buffer_Stdev)
        self.SAXSdata= numpy.vstack((self.hplc_run.q, Isub, StdErr)).T
        with open(self.subtracted, "w") as outfile:
            #print "Writing file"
            numpy.savetxt(outfile, self.SAXSdata)
        self.xsDataResult.subtractedCurve = XSDataFile(XSDataString(self.subtracted))
        self.frame.subtracted = self.subtracted
        
        if self.subtracted and os.path.exists(self.subtracted) and self.ATSASRg:
            self.edPluginAutoRg = self.loadPlugin(self.strControlledPluginAutoRg)
            self.edPluginAutoRg.dataInput = XSDataInputAutoRg(inputCurve=[XSDataFile(XSDataString(self.subtracted))])
            self.edPluginAutoRg.connectSUCCESS(self.doSuccessAutoRg)
            self.edPluginAutoRg.connectFAILURE(self.doFailureAutoRg)
            self.edPluginAutoRg.executeSynchronous()
        try:
            #print self.subtracted, autoRg(self.SAXSdata)
            self.frame.RgF, self.frame.Rg_StdevF, self.frame.I0F, self.frame.I0_StdevF, self.frame.Rg_imin, self.frame.Rg_imax  =   autoRg(self.SAXSdata)
        except InsufficientDataError:
            print "Not enough usable data to run autorg"
        except SystemError as SE:
            print "System Error"
            print SE
        except Exception as Err:
            print Err
            raise 
        else: 
            if self.ATSASRg == False:
                self.frame.Rg = max(0,self.frame.RgF)
                self.frame.Rg_Stdev =  max(0,self.frame.Rg_StdevF)
                self.frame.I0 =  max(0,self.frame.I0F)
                self.frame.I0_Stdev =  max(0,self.frame.I0_StdevF)
                self.frame.quality = 0

                """
                Calculate the invariants Vc and Qr from the Rambo&Tainer 2013 Paper,
                also the the mass estimate based on Qr for proteins
                """
            if self.subtracted and os.path.exists(self.subtracted):
                self.subtracted_data = numpy.loadtxt(self.subtracted)
                if self.subtracted_data is not None and\
                    self.frame.RgF and self.frame.Rg_StdevF and self.frame.I0F and self.frame.I0_StdevF:
                    dictRTI = RamboTainerInvariant(self.subtracted_data, self.frame.RgF,
                                                   self.frame.Rg_StdevF, self.frame.I0F,
                                                   self.frame.I0_StdevF, self.frame.Rg_imin)
    #             {'Vc': vc[0], 'dVc': vc[1], 'Qr': qr, 'dQr': dqr, 'mass': mass, 'dmass': dmass}
                    self.frame.Vc = dictRTI.get("Vc")
                    self.frame.Vc_Stdev = dictRTI.get("dVc")
                    self.frame.Qr = dictRTI.get("Qr")
                    self.frame.Qr_Stdev = dictRTI.get("dQ")
                    self.frame.mass = dictRTI.get("mass")
                    self.frame.mass_Stdev = dictRTI.get("dmass")
                    xsdRTI = XSDataRamboTainer(vc=XSDataDouble(self.frame.Vc),
                                               qr=XSDataDouble(self.frame.Qr),
                                               mass=XSDataDouble(self.frame.mass),
                                               dvc=XSDataDouble(self.frame.Vc_Stdev),
                                               dqr=XSDataDouble(self.frame.Qr_Stdev),
                                               dmass=XSDataDouble(self.frame.mass_Stdev))
                    self.xsDataResult.rti = xsdRTI



    def postProcess(self, _edObject=None):
        """
        after process
        """
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.postProcess")
        if self.hplc_run.buffer:
            self.xsDataResult.bufferCurve = XSDataFile(XSDataString(self.hplc_run.buffer))
        if self.curve:
            self.xsDataResult.integratedCurve = XSDataFile(XSDataString(self.curve))
        if self.subtracted:
            self.xsDataResult.subtractedCurve = XSDataFile(XSDataString(self.subtracted))

    def finallyProcess(self, _edObject=None):
        EDPluginControl.finallyProcess(self)
        executiveSummary = os.linesep.join(self.lstExecutiveSummary)
        self.xsDataResult.status = XSDataStatus(executiveSummary=XSDataString(executiveSummary))
        self.dataOutput = self.xsDataResult
        if self.frame:
            self.frame.processing = False


    def average_buffers(self):
        """
        Average out all buffers
        """

        self.DEBUG("Averaging")
        print "Avergaring"
        nb_frames = len(self.hplc_run.for_buffer)
        print nb_frames
        self.hplc_run.for_buffer.sort()
        filename = self.hplc_run.first_curve[::-1].split("_", 1)[1][::-1] + \
                        "_buffer_aver_%04i.dat" % nb_frames
        self.lstExecutiveSummary.append("Averaging out buffer to %s: %s" % (filename,
                                        ", ".join([str(i) for i in self.hplc_run.for_buffer])))

        q = self.hplc_run.q
        I = self.hplc_run.for_buffer_sum_I / nb_frames
        err = numpy.sqrt(self.hplc_run.for_buffer_sum_sigma2) / nb_frames
        m = numpy.vstack((q, I, err))

        with open(filename, "w") as outfile:
            # print "Saving buffer"
            numpy.savetxt(outfile, m.T)
        with self._sem:
            self.hplc_run.buffer = filename
            self.hplc_run.buffer_I = I
            self.hplc_run.buffer_Stdev = err

    def doSuccessProcessOneFile(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.doSuccessProcessOneFile")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsHPLCv1_5.doSuccessProcessOneFile")
        if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
            self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
        output = _edPlugin.dataOutput
        if not output.integratedCurve:
            strErr = "Edna plugin ProcessOneFile did not produce integrated curve"
            self.ERROR(strErr)
            self.lstExecutiveSummary.append(strErr)
            self.setFailure()
            return
        self.curve = output.integratedCurve.path.value
        if not os.path.exists(self.curve):
            strErr = "Edna plugin ProcessOneFile: integrated curve not on disk !!"
            self.ERROR(strErr)
            self.lstExecutiveSummary.append(strErr)
            self.setFailure()
            return
        self.xsDataResult.integratedCurve = output.integratedCurve
        self.xsDataResult.normalizedImage = output.normalizedImage
        self.xsDataResult.dataI = output.dataI
        self.xsDataResult.dataQ = output.dataQ
        self.xsDataResult.dataStdErr = output.dataStdErr
        self.intensity = EDUtilsArray.xsDataToArray(output.dataI)
        self.stdError = EDUtilsArray.xsDataToArray(output.dataStdErr)
        if output.experimentSetup and output.experimentSetup.timeOfFrame:
            startTime = output.experimentSetup.timeOfFrame.value
        else:
            try:
                startTime = float(fabio.openheader(self.dataInput.rawImage.path.value).header["time_of_day"])
            except Exception:
                self.ERROR("Unable to read time_of_day in header of %s" % self.dataInput.rawImage.path.value)
                startTime = 0

        if not self.hplc_run.first_curve:
            with self._sem:
                if True: #not self.hplc_run.first_curve:
                    # Populate the buffer with the first curve if needed
                    self.hplc_run.first_curve = self.curve
                    self.hplc_run.start_time = startTime
                    self.hplc_run.q = EDUtilsArray.xsDataToArray(output.dataQ)
                    self.hplc_run.size = self.hplc_run.q.size
                    self.hplc_run.buffer_I = self.intensity
                    self.hplc_run.buffer_Stdev = self.stdError
                    self.hplc_run.firstCurveIntensity = self.intensity
                    self.hplc_run.for_buffer_sum_I = self.intensity
                    self.hplc_run.for_buffer_sum_sigma2 = self.stdError ** 2
                    self.hplc_run.for_buffer.append(self.frameId)
        self.frame.curve = self.curve
        self.frame.time = startTime
        self.xsDataResult.timeStamp = XSDataTime(value=(startTime - self.hplc_run.start_time))
#         self.calcIntensity()
        self.frame.sum_I = float(self.intensity.sum())
        self.xsDataResult.summedIntensity = XSDataDouble(value=self.frame.sum_I)

    def doFailureProcessOneFile(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.doFailureProcessOneFile")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsHPLCv1_5.doFailureProcessOneFile")
        if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
            self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
        else:
            self.lstExecutiveSummary.append("Edna plugin ProcessOneFile failed.")
        self.setFailure()

    def doSuccessAutoRg(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.doSuccessAutoRg")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsHPLCv1_5.doSuccessAutoRg")
        if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
            self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
        try:
            rg = _edPlugin.dataOutput.autoRgOut[0]
        except:
            self.ERROR("AutoRg returned nothing !")
        else:
            if rg.rg:
                self.frame.Rg = rg.rg.value
            if rg.rgStdev:
                self.frame.Rg_Stdev = rg.rgStdev.value
            if rg.i0:
                self.frame.I0 = rg.i0.value
            if rg.i0Stdev:
                self.frame.I0_Stdev = rg.i0Stdev.value
            if rg.quality:
                self.frame.quality = rg.quality.value
            self.xsDataResult.autoRg = rg
      
     
        """
        Calculate the invariants Vc and Qr from the Rambo&Tainer 2013 Paper,
        also the the mass estimate based on Qr for proteins
        """
        if self.subtracted and os.path.exists(self.subtracted):
            self.subtracted_data = numpy.loadtxt(self.subtracted)
            if self.subtracted_data is not None and\
                self.frame.Rg and self.frame.Rg_Stdev and self.frame.I0 and self.frame.I0_Stdev:
                dictRTI = RamboTainerInvariant(self.subtracted_data, self.frame.Rg,
                                               self.frame.Rg_Stdev, self.frame.I0,
                                               self.frame.I0_Stdev, rg.firstPointUsed.value)
#             {'Vc': vc[0], 'dVc': vc[1], 'Qr': qr, 'dQr': dqr, 'mass': mass, 'dmass': dmass}
                self.frame.Vc = dictRTI.get("Vc")
                self.frame.Vc_Stdev = dictRTI.get("dVc")
                self.frame.Qr = dictRTI.get("Qr")
                self.frame.Qr_Stdev = dictRTI.get("dQ")
                self.frame.mass = dictRTI.get("mass")
                self.frame.mass_Stdev = dictRTI.get("dmass")
                xsdRTI = XSDataRamboTainer(vc=XSDataDouble(self.frame.Vc),
                                           qr=XSDataDouble(self.frame.Qr),
                                           mass=XSDataDouble(self.frame.mass),
                                           dvc=XSDataDouble(self.frame.Vc_Stdev),
                                           dqr=XSDataDouble(self.frame.Qr_Stdev),
                                           dmass=XSDataDouble(self.frame.mass_Stdev))
                self.xsDataResult.rti = xsdRTI

    def doFailureAutoRg(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsHPLCv1_5.doFailureAutoRg")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsHPLCv1_5.doFailureAutoRg")
        strErr = "Edna plugin AutoRg failed."
        if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
            self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
            self.lstExecutiveSummary.append(strErr)
        else:
            self.lstExecutiveSummary.append(strErr)
#         self.setFailure()

    #def doSuccessDatCmp(self, _edPlugin=None):
    def datcmp(self,intensity1,intensity2):
        self.DEBUG("DatCmp")
        #self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsHPLCv1_5.doSuccessDatCmp")
        #if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
        #   self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
        #if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.fidelity:
        #    fidelity = _edPlugin.dataOutput.fidelity.value
        try: 
            cor = gof(intensity1,intensity2)
            fidelity = cor.P
        except Exception as err:
            strErr = "No fidelity in output of datcmp"
            self.error(strErr)
            self.lstExecutiveSummary.append(strErr)
            self.lstExecutiveSummary.append(err)
            self.setFailure()
            fidelity = 0

        if not self.hplc_run.deleted and self.hplc_run.buffer is None:
            if fidelity > self.SIMILARITY_THRESHOLD_SAMPLE:
            
                self.isBuffer = True
                if fidelity > self.SIMILARITY_THRESHOLD_BUFFER:
                    # self.frame.I = EDUtilsArray.xsDataToArray(self.xsDataResult.dataI)
                    # self.frame.q = EDUtilsArray.xsDataToArray(self.xsDataResult.dataQ)
                    # self.frame.err = EDUtilsArray.xsDataToArray(self.xsDataResult.dataStdErr)
                    with self.buffer_sem:
                        self.hplc_run.for_buffer.append(self.frameId)
                        if self.hplc_run.for_buffer_sum_I is None:
                            self.hplc_run.for_buffer_sum_I = self.intensity
                        else:
                            self.hplc_run.for_buffer_sum_I = self.hplc_run.for_buffer_sum_I + self.intensity
                        if self.hplc_run.for_buffer_sum_sigma2 is None:
                            self.hplc_run.for_buffer_sum_sigma2 = self.stdError ** 2
                        else:
                            self.hplc_run.for_buffer_sum_sigma2 += self.stdError ** 2
            else:
                with self.buffer_sem:
                    if self.hplc_run.buffer is None:
                        self.average_buffers()
        elif fidelity > self.SIMILARITY_THRESHOLD_SAMPLE:
            self.isBuffer = True



#     def doFailureDatCmp(self, _edPlugin=None):
#         self.DEBUG("EDPluginBioSaxsHPLCv1_5.doFailureDatCmp")
#         self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsHPLCv1_5.doFailureDatCmp")
#         if _edPlugin and _edPlugin.dataOutput and _edPlugin.dataOutput.status and _edPlugin.dataOutput.status.executiveSummary:
#             self.lstExecutiveSummary.append(_edPlugin.dataOutput.status.executiveSummary.value)
#         else:
#             self.lstExecutiveSummary.append("Edna plugin DatCmp failed.")
#         self.setFailure()

