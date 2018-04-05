# coding: utf-8
#
#    Project: BioSaxs
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) ESRF
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
__copyright__ = "ESRF"
__date__ = "05/04/2018"
__status__ = "production"

import os
import time
import copy
try:
    from Queue import LifoQueue
except:
    from queue import LifoQueue
from EDUtilsParallel        import EDUtilsParallel
from EDVerbose              import EDVerbose
from EDPluginControl        import EDPluginControl
from EDFactoryPlugin        import edFactoryPlugin
from EDThreading            import Semaphore
from EDUtilsUnit            import EDUtilsUnit
from EDUtilsArray           import EDUtilsArray
edFactoryPlugin.loadModule("XSDataEdnaSaxs")
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
edFactoryPlugin.loadModule("XSDataWaitFilev1_0")
from XSDataWaitFilev1_0     import XSDataInputWaitFile
from XSDataBioSaxsv1_0      import XSDataInputBioSaxsProcessOneFilev1_0, XSDataResultBioSaxsProcessOneFilev1_0
from XSDataCommon           import XSDataStatus, XSDataString, XSDataFile, XSDataInteger, XSDataTime
import matplotlib
matplotlib.use("Agg")  # unless pyFAI initializes another backend !
import fabio
import numpy
import pyFAI

if [int(i) for i in pyFAI.version.split(".")[:2]] < [0, 13]:
    EDVerbose.ERROR("Too old version of pyFAI detected ... expect to fail !")


class EDPluginBioSaxsProcessOneFilev1_6(EDPluginControl):
    """
    Control plugin that does the same as previously without sub-plugin call ...
    except WaitFile which is still called.

    Nota normalization is done AFTER integration not before as previously
    
    New in version 1.4: return I/Q/stderr via XSDataArrays
     
    """
    cpWaitFile = "EDPluginWaitFile"
    integrator = pyFAI.AzimuthalIntegrator()
    integrator.wavelength = 1e-10
    available_integrators = LifoQueue()
    CONF_DUMMY_PIXEL_VALUE = "DummyPixelValue"
    CONF_DUMMY_PIXEL_DELTA = "DummyPixelDelta"
    CONF_OPENCL_DEVICE = "DeviceType"
    CONF_NUMBER_OF_BINS = "NumberOfBins"
    __configured = False
    dummy = -2
    delta_dummy = 1.1
    semaphore = Semaphore()
    if pyFAI.opencl.ocl is None:
        METHOD = "fullsplit_csr"
    else:
        METHOD = "fullsplit_csr_ocl_gpu"
    number_of_bins = None

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsProcessOneFilev1_0)
        self.__edPluginWaitFile = None
        self.rawImage = None
        self.rawImageSize = XSDataInteger(1024)
        self.normalizedImage = None
        self.integratedCurve = None
        self.integratedImage = None
        self.lstExecutiveSummary = []
        self.sample = None
        self.experimentSetup = None
        self.integrator_config = {}
        self.normalization_factor = None
        self.detector = None
        self.xsDataResult = XSDataResultBioSaxsProcessOneFilev1_0()

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.rawImage, "No raw image provided")
        self.checkMandatoryParameters(self.dataInput.sample, "No sample information provided")
        self.checkMandatoryParameters(self.dataInput.experimentSetup, "No experimental setup provided")
        self.checkMandatoryParameters(self.dataInput.integratedCurve, "Missing integratedCurve")

    def configure(self):
        """
        Configures the plugin from the configuration file with the following parameters:
         - DummyPixelValue: the value be assigned to dummy pixels.
         - DummyPixelDelta: the value be assigned to delta dummy.
         - DeviceType: "lut_ocl_1,3" will select device #3 on first platform #1 
        """
        EDPluginControl.configure(self)
        if not self.__configured:
            with self.semaphore:
                if not self.__configured:
                    self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.configure")
                    dummy = self.config.get(self.CONF_DUMMY_PIXEL_VALUE)
                    if dummy is None:
                        strMessage = 'EDPluginBioSaxsProcessOneFilev1_6.configure: %s Configuration parameter missing: \
            %s, defaulting to "%s"' % (self.getBaseName(), self.CONF_DUMMY_PIXEL_VALUE, self.dummy)
                        self.WARNING(strMessage)
                        self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                    else:
                        self.__class__.dummy = float(dummy)
                    ddummy = self.config.get(self.CONF_DUMMY_PIXEL_DELTA)
                    if ddummy is None:
                        strMessage = 'EDPluginBioSaxsProcessOneFilev1_6.configure: %s Configuration parameter missing: \
            %s, defaulting to "%s"' % (self.getBaseName(), self.CONF_DUMMY_PIXEL_DELTA, self.delta_dummy)
                        self.WARNING(strMessage)
                        self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                    else:
                        self.__class__.delta_dummy = float(ddummy)
                    method = self.config.get(self.CONF_OPENCL_DEVICE)
                    if method is None:
                        strMessage = 'EDPluginBioSaxsProcessOneFilev1_6.configure: %s Configuration parameter missing: \
            %s, defaulting to "%s"' % (self.getBaseName(), self.CONF_OPENCL_DEVICE, self.METHOD)
                        self.WARNING(strMessage)
                        self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                    else:
                        self.__class__.METHOD = method
                    number_of_bins = self.config.get(self.CONF_NUMBER_OF_BINS)
                    if number_of_bins is None:
                        strMessage = 'EDPluginBioSaxsProcessOneFilev1_6.configure: %s Configuration parameter missing: \
            %s, defaulting to max(image.shape)' % (self.getBaseName(), self.CONF_NUMBER_OF_BINS)
                        self.WARNING(strMessage)
                        self.addErrorWarningMessagesToExecutiveSummary(strMessage)
                    else:
                        self.__class__.number_of_bins = number_of_bins
                    self.__class__.__configured = True

    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.preProcess")
        self.__edPluginWaitFile = self.loadPlugin(self.cpWaitFile)
        if self.dataInput.rawImageSize is not None:
            self.rawImageSize = self.dataInput.rawImageSize
        if self.dataInput.rawImage is not None:
            self.rawImage = self.dataInput.rawImage.path.value
        self.integratedCurve = self.dataInput.integratedCurve.path.value
        curveDir = os.path.dirname(self.integratedCurve)
        if not os.path.isdir(curveDir):
            try:
                os.mkdir(curveDir)
            except OSError:
                # could occur in race condition ...
                pass

        self.sample = self.dataInput.sample
        self.experimentSetup = self.dataInput.experimentSetup
        detector = self.experimentSetup.detector.value
        if detector.lower() == "pilatus":
            self.detector = pyFAI.detector_factory("Pilatus1M")
        else:
            self.detector = pyFAI.detector_factory(detector)
            try:
                self.detector.pixel1 = self.experimentSetup.pixelSize_2.value,  # flip X,Y
                self.detector.pixel2 = self.experimentSetup.pixelSize_1.value,  # flip X,Y
            except Exception as err:
                self.WARNING("in setting pixel size: %s" % err)

        i0 = self.experimentSetup.beamStopDiode.value
        normalization_factor = self.experimentSetup.normalizationFactor.value
        if normalization_factor == 0:
            warn = "normalization_factor is Null --> If we are testing, this is OK, else investigate !!!"
            self.lstExecutiveSummary.append(warn)
            self.warning(warn)
            normalization_factor = 1.
        if i0 == 0:
            warn = "beamStopDiode is Null --> If we are testing, this is OK, else investigate !!!"
            self.lstExecutiveSummary.append(warn)
            self.warning(warn)
            i0 = 1.0

        self.normalization_factor = i0 / normalization_factor

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.process")

        xsd = XSDataInputWaitFile(expectedFile=XSDataFile(XSDataString(self.rawImage)),
                                  expectedSize=self.rawImageSize,
                                  timeOut=XSDataTime(30))
        self.__edPluginWaitFile.setDataInput(xsd)
        self.__edPluginWaitFile.connectSUCCESS(self.doSuccessWaitFile)
        self.__edPluginWaitFile.connectFAILURE(self.doFailureWaitFile)
        self.__edPluginWaitFile.executeSynchronous()
        if self.isFailure():
            return

        self.xsDataResult.sample = self.sample
        self.xsDataResult.experimentSetup = self.experimentSetup

        res = self.integrate()
        self.write3ColumnAscii(res, self.integratedCurve)
        self.xsDataResult.dataQ = EDUtilsArray.arrayToXSData(res.radial)
        self.xsDataResult.dataI = EDUtilsArray.arrayToXSData(res.intensity)
        self.xsDataResult.dataStdErr = EDUtilsArray.arrayToXSData(res.sigma)

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.postProcess")
        # Create some output data
        if os.path.exists(self.integratedCurve):
            self.xsDataResult.integratedCurve = XSDataFile(XSDataString(self.integratedCurve))

        self.xsDataResult.status = XSDataStatus(executiveSummary=XSDataString(os.linesep.join(self.lstExecutiveSummary)))
        self.setDataOutput(self.xsDataResult)

    def integrate(self):
        with fabio.open(self.rawImage) as img:
            # img = fabio.open(self.rawImage)
            number_of_bins = self.number_of_bins or max(img.dim1, img.dim2)
            if "Date" in img.header:
                self.experimentSetup.timeOfFrame = XSDataTime(time.mktime(time.strptime(img.header["Date"], "%a %b %d %H:%M:%S %Y")))

            new_integrator = pyFAI.AzimuthalIntegrator(detector=self.detector)
            new_integrator.setFit2D(self.experimentSetup.detectorDistance.value * 1000,
                                    self.experimentSetup.beamCenter_1.value,
                                    self.experimentSetup.beamCenter_2.value)
            new_integrator.wavelength = EDUtilsUnit.getSIValue(self.experimentSetup.wavelength)

            with self.__class__.semaphore:
                if self.available_integrators.empty():
                    integrator = copy.copy(self.__class__.integrator)
                else:
                    integrator = self.available_integrators.get()
                if ((str(new_integrator) != str(integrator)) or
                        (integrator.detector._maskfile != self.experimentSetup.maskFile.path.value)):
                    self.screen("Resetting PyFAI integrator")
                    new_integrator.detector.mask = self.calc_mask()
                    new_integrator.detector._maskfile = self.experimentSetup.maskFile.path.value
                    integrator = self.__class__.integrator = new_integrator

            res_tuple = integrator.integrate1d(img.data, number_of_bins,
                                               correctSolidAngle=True,
                                               dummy=self.dummy, delta_dummy=self.delta_dummy,
                                               filename=None,
                                               error_model="poisson",
                                               radial_range=None, azimuth_range=None,
                                               polarization_factor=0.99, dark=None, flat=None,
                                               method=self.METHOD, unit="q_nm^-1", safe=False,
                                               normalization_factor=self.normalization_factor
                                               )
            self.available_integrators.put(integrator)
            self.lstExecutiveSummary.append("Azimuthal integration of raw image '%s'-->'%s'." % (self.rawImage, self.integratedCurve))
            valid_bins = abs(res_tuple.intensity - self.dummy) > self.delta_dummy
            short_tuple = res_tuple.__class__(res_tuple.radial[valid_bins],
                                              res_tuple.intensity[valid_bins],
                                              res_tuple.sigma[valid_bins])
            return short_tuple

    def calc_mask(self, maskfile=None):
        """
        Merge the natural mask from the detector with the user proided one.
        @return: numpy array with the mask
        """
        if maskfile is None:
            maskfile = self.experimentSetup.maskFile.path.value
        mask = fabio.open(maskfile).data

        detector_mask = self.detector.calc_mask()
        shape0, shape1 = detector_mask.shape
        if detector_mask.shape == mask.shape:
            mask = numpy.logical_or(mask, detector_mask)
        else:
            # crop the user defined mask
            mask = numpy.logical_or(mask[:shape0, :shape1], detector_mask)
        return mask

    def write3ColumnAscii(self, res, outputCurve="output.dat", hdr="#", linesep=os.linesep):
        """
        @param res: named tuple of numpy array containing Scattering vector, Intensity and deviation
        @param outputCurve: name of the 3-column ascii file to be written
        @param hdr: header mark, usually '#'


Adam Round explicitelly asked for (email from Date: Tue, 04 Oct 2011 15:22:29 +0200) :
Modification from:

# BSA buffer
# Sample c= 0.0 mg/ml (these two lines are required for current DOS pipeline and can be cleaned up once we use EDNA to get to ab-initio models)
#
# Sample environment:
# Detector = Pilatus 1M
# PixelSize_1 = 0.000172
# PixelSize_2 = 6.283185 (I think it could avoid confusion if we give teh actual pixel size as 0.000172 for X and Y and not to give the integrated sizes. Also could there also be a modification for PixelSize_1 as on the diagonal wont it be the hypotenuse (0.000243)? and thus will be on average a bit bigger than 0.000172)
#
# title = BSA buffer
# Frame 7 of 10
# Time per frame (s) = 10
# SampleDistance = 2.43
# WaveLength = 9.31e-11
# Normalization = 0.0004885
# History-1 = saxs_angle +pass -omod n -rsys normal -da 360_deg -odim = 1 /data/id14eh3/inhouse/saxs_pilatus/Adam/EDNAtests/2d/dumdum_008_07.edf/data/id14eh3/inhouse/saxs_pilatus/Adam/EDNAtests/misc/dumdum_008_07.ang
# DiodeCurr = 0.0001592934
# MachCurr = 163.3938
# Mask = /data/id14eh3/archive/CALIBRATION/MASK/Pcon_01Jun_msk.edf
# SaxsDataVersion = 2.40
#
# N 3
# L q*nm  I_BSA buffer  stddev
#
# Sample Information:
# Storage Temperature (degrees C): 4
# Measurement Temperature (degrees C): 10
# Concentration: 0.0
# Code: BSA
s-vector Intensity Error
s-vector Intensity Error
s-vector Intensity Error
s-vector Intensity Error
        """
        hdr = str(hdr)
        headers = []
        if self.sample.comments is not None:
            headers.append(hdr + " " + self.sample.comments.value)
        else:
            headers.append(hdr)
        if self.sample.concentration is not None:
            headers.append(hdr + " Sample c= %s mg/ml" % self.sample.concentration.value)
        else:
            headers.append(hdr + " Sample c= -1  mg/ml")
        headers += [hdr, hdr + " Sample environment:"]
        if self.experimentSetup.detector is not None:
            headers.append(hdr + " Detector = %s" % self.experimentSetup.detector.value)
        if self.experimentSetup.pixelSize_1 is not None:
            headers.append(hdr + " PixelSize_1 = %s" % self.experimentSetup.pixelSize_1.value)
        if self.experimentSetup.pixelSize_2 is not None:
            headers.append(hdr + " PixelSize_2 = %s" % self.experimentSetup.pixelSize_2.value)
        headers.append(hdr)
        if self.sample.comments is not None:
            headers.append(hdr + " title = %s" % self.sample.comments.value)
        if (self.experimentSetup.frameNumber is not None) and\
           (self.experimentSetup.frameMax is not None):
            headers.append(hdr + " Frame %s of %s" % (self.experimentSetup.frameNumber.value, self.experimentSetup.frameMax.value))
        if self.experimentSetup.exposureTime is not None:
            headers.append(hdr + " Time per frame (s) = %s" % self.experimentSetup.exposureTime.value)
        if self.experimentSetup.detectorDistance is not None:
            headers.append(hdr + " SampleDistance = %s" % self.experimentSetup.detectorDistance.value)
        if self.experimentSetup.wavelength is not None:
            headers.append(hdr + " WaveLength = %s" % self.experimentSetup.wavelength.value)
        if self.experimentSetup.normalizationFactor is not None:
            headers.append(hdr + " Normalization = %s" % self.experimentSetup.normalizationFactor.value)
        if self.experimentSetup.beamStopDiode is not None:
            headers.append(hdr + " DiodeCurr = %s" % self.experimentSetup.beamStopDiode.value)
        if self.experimentSetup.machineCurrent is not None:
            headers.append(hdr + " MachCurr = %s" % self.experimentSetup.machineCurrent.value)
        if self.experimentSetup.maskFile is not None:
            headers.append(hdr + " Mask = %s" % self.experimentSetup.maskFile.path.value)
        headers.append(hdr)
        headers.append(hdr + " N 3")
        if self.sample.comments is not None:
            headers.append(hdr + " L q*nm  I_%s  stddev" % self.sample.comments.value)
        else:
            headers.append(hdr + " L q*nm  I_  stddev")
        headers.append(hdr)
        headers.append(hdr + " Sample Information:")
        if self.experimentSetup.storageTemperature is not None:
            headers.append(hdr + " Storage Temperature (degrees C): %s" % self.experimentSetup.storageTemperature.value)
        if self.experimentSetup.exposureTemperature is not None:
            headers.append(hdr + " Measurement Temperature (degrees C): %s" % self.experimentSetup.exposureTemperature.value)

        if self.sample.concentration is not None:
            headers.append(hdr + " Concentration: %s" % self.sample.concentration.value)
        else:
            headers.append(hdr + " Concentration: -1")
        if self.sample.code is not None:
            headers.append(hdr + " Code: %s" % self.sample.code.value)
        else:
            headers.append(hdr + " Code: ")

        with open(outputCurve, "w") as f:
            f.writelines(linesep.join(headers))
            f.write(linesep)

            if res.sigma is None:
                data = ["%14.6e %14.6e " % (q, I)
                        for q, I in zip(res.radial, res.intensity)]
                        # 3if abs(I - self.dummy) > self.delta_dummy]
            else:
                data = ["%14.6e %14.6e %14.6e" % (q, I, std)
                        for q, I, std in zip(res.radial, res.intensity, res.sigma)]
                        # if abs(I - self.dummy) > self.delta_dummy]
            data.append("")
            f.writelines(linesep.join(data))
            f.flush()

    def doSuccessWaitFile(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.doSuccessWaitFile")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginBioSaxsProcessOneFilev1_6.doSuccessWaitFile")

    def doFailureWaitFile(self, _edPlugin=None):
        self.DEBUG("EDPluginBioSaxsProcessOneFilev1_6.doFailureWaitFile")
        self.retrieveFailureMessages(_edPlugin, "EDPluginBioSaxsProcessOneFilev1_6.doFailureWaitFile")
        self.lstExecutiveSummary.append("Timeout in waiting for file '%s'" % (self.rawImage))
        self.setFailure()
