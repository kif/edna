# coding: utf-8
#
#    Project: Edna Saxs
#             http://www.edna-site.org
#
#    Copyright (C) 2012-2015 ESRF
#
#    Principal author: Jerome Kieffer
#
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

from __future__ import print_function, with_statement

__authors__ = ["Jérôme Kieffer"]
__license__ = "GPLv3+"
__copyright__ = "ESRF"
__date__ = "28/08/2015"
__status__ = "Development"

import os, gc
import numpy
try: 
    import matplotlib
except ImportError:
    plt = None
else:    
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

from EDPluginControl import EDPluginControl
from EDFactoryPlugin import edFactoryPlugin
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
edFactoryPlugin.loadModule("XSDataEdnaSaxs")
from XSDataEdnaSaxs import XSDataInputSaxsAnalysis, XSDataResultSaxsAnalysis, \
                           XSDataInputAutoRg, XSDataInputDatGnom, XSDataInputDatPorod#, XSDataRamboTainer
from XSDataBioSaxsv1_0 import XSDataRamboTainer                           
from XSDataCommon import XSDataString, XSDataFile, XSDataInteger, XSDataStatus, XSDataDouble
from saxs_plotting import scatterPlot, guinierPlot, kartkyPlot, densityPlot
from EDUtilsBioSaxs import RamboTainerInvariant

class EDPluginControlSaxsAnalysisv1_0(EDPluginControl):
    """
    Executes the pipeline:
    * AutoRg -> Extract the Guinier region and measure Rg, I0
    * DatGnom -> transformation from reciprocal to direct space. measure Dmax.
    * DatPorod -> calculates the volume of the protein using the porod formula.
    """
    cpAutoRg = "EDPluginExecAutoRgv1_1"
    cpDatGnom = "EDPluginExecDatGnomv1_1"
    cpDatPorod = "EDPluginExecDatPorodv1_0"

    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputSaxsAnalysis)
        self.edPluginAutoRg = None
        self.edPluginDatPorod = None
        self.edPluginDatGnom = None
        self.scatterFile = None
        self.gnomFile = None
        self.autoRg = None
        self.gnom = None
        self.rti = None
        self.xVolume = None
        self.xsDataResult = XSDataResultSaxsAnalysis()

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.scatterCurve, "No scattering curve provided")


    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.preProcess")
        self.scatterFile = self.dataInput.scatterCurve.path.value
        if self.dataInput.gnomFile is not None:
            self.gnomFile = self.dataInput.gnomFile.path.value
        else:
            self.gnomFile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + ".out")
        self.autoRg = self.dataInput.autoRg


    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.process")
        if self.autoRg is None:
            self.edPluginAutoRg = self.loadPlugin(self.cpAutoRg)
            self.edPluginAutoRg.dataInput = XSDataInputAutoRg(inputCurve=[self.dataInput.scatterCurve])
            self.edPluginAutoRg.connectSUCCESS(self.doSuccessRg)
            self.edPluginAutoRg.connectFAILURE(self.doFailureRg)
            self.edPluginAutoRg.executeSynchronous()

        if self.autoRg is None or self.autoRg.rg.value == 0.0:
            self.setFailure()
        else:
            #self.screen("autorg %s scatter %s"%(self.autoRg, self.scatterFile))
            self.calculateRTI(self.autoRg, self.scatterFile)
            #self.screen("RTI: %s"%self.rti)
        
        if self.isFailure():
            return

        self.edPluginDatGnom = self.loadPlugin(self.cpDatGnom)
        self.edPluginDatGnom.dataInput = XSDataInputDatGnom(inputCurve=self.dataInput.scatterCurve,
                                             output=XSDataFile(XSDataString(self.gnomFile)),
                                             rg=self.autoRg.rg,
                                             skip=XSDataInteger(self.autoRg.firstPointUsed.value - 1))
        self.edPluginDatGnom.connectSUCCESS(self.doSuccessGnom)
        self.edPluginDatGnom.connectFAILURE(self.doFailureGnom)
        self.edPluginDatGnom.executeSynchronous()

        if self.gnom is None:
            return

        self.edPluginDatPorod = self.loadPlugin(self.cpDatPorod)
        self.edPluginDatPorod.dataInput = XSDataInputDatPorod(gnomFile=XSDataFile(XSDataString(self.gnomFile)))
        self.edPluginDatPorod.connectSUCCESS(self.doSuccessPorod)
        self.edPluginDatPorod.connectFAILURE(self.doFailurePorod)
        self.edPluginDatPorod.execute()

        if self.dataInput.graphFormat:
            ext = self.dataInput.graphFormat.value
            if not ext.startswith("."):
                ext = "." + ext
            try:
                guinierfile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + "-Guinier" + ext)
#                 guinierplot = guinierPlot(self.scatterFile, unit="nm",
#                                        filename=guinierfile, format=ext[1:])
                guinierplot = guinierPlot(self.scatterFile, unit="nm",
                                       filename=guinierfile, format=ext[1:], first_point=self.autoRg.firstPointUsed.value,
                                       last_point=self.autoRg.lastPointUsed.value)
                guinierplot.clf()
                if plt:
                    plt.close(guinierplot)
            except Exception as error:
                self.ERROR("EDPluginControlSaxsAnalysisv1_0 in guinierplot: %s"%error)
            else:
                self.xsDataResult.guinierPlot = XSDataFile(XSDataString(guinierfile))

            try:
                kratkyfile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + "-Kratky" + ext)
                kratkyplot = kartkyPlot(self.scatterFile, unit="nm",
                                           filename=kratkyfile, format=ext[1:], Rg=self.autoRg.rg.value, I0=self.autoRg.i0.value)
                kratkyplot.clf()
                if plt:
                    plt.close(kratkyplot)
            except Exception as error:
                self.ERROR("EDPluginControlSaxsAnalysisv1_0 in kratkyplot: %s"%error)
            else:
                self.xsDataResult.kratkyPlot = XSDataFile(XSDataString(kratkyfile))
#             if self.autoRg is not None:    
#                 if self.autoRg.i0.value > 0 and self.autoRg.rg.value > 0:
#                     try:
#                         kratkyRgfile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + "-KratkyRg" + ext)
#                         kratkyRgplot = kratkyRgPlot(self.scatterFile, self.autoRg.i0.value, self.autoRg.rg.value,
#                                                    filename=kratkyRgfile, format=ext[1:])
#                         kratkyRgplot.clf()
#                         if plt:
#                             plt.close(kratkyRgplot)
#                     except Exception as error:
#                         self.ERROR("EDPluginControlSaxsAnalysisv1_1 in kratkyRgplot: %s" % error)
#                     else:
#                         self.xsDataResult.kratkyRgPlot = XSDataFile(XSDataString(kratkyRgfile))       
#             if self.autoRg is not None and self.rti is not None:            
#                 if self.autoRg.i0.value > 0 and self.rti.vc.value > 0:
#                     try:
#                         kratkyVcfile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + "-KratkyVc" + ext)
#                         kratkyVcplot = kratkyVcPlot(self.scatterFile, self.autoRg.i0.value, self.rti.vc.value,
#                                                         filename=kratkyVcfile, format=ext[1:])
#                         kratkyVcplot.clf()
#                         if plt:
#                             plt.close(kratkyVcplot)
#                     except Exception as error:
#                         self.ERROR("EDPluginControlSaxsAnalysisv1_1 in kratkyVcplot: %s" % error)
#                     else:
#                         self.xsDataResult.kratkyVcPlot = XSDataFile(XSDataString(kratkyVcfile))

            try:
                scatterplotfile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + "-scattering" + ext)
                scatterplot = scatterPlot(self.scatterFile, unit="nm", gnomfile=self.gnomFile,
                                           filename=scatterplotfile, format=ext[1:])
                scatterplot.clf()
                if plt:
                    plt.close(scatterplot)
            except Exception as error:
                self.ERROR("EDPluginControlSaxsAnalysisv1_0 in scatterplot: %s"%error)
            else:
                self.xsDataResult.scatterPlot = XSDataFile(XSDataString(scatterplotfile))
            try:
                densityplotfile = os.path.join(self.getWorkingDirectory(), os.path.basename(self.scatterFile).split(".")[0] + "-density" + ext)
                densityplot = densityPlot(gnomfile=self.gnomFile, unit="nm",
                                           filename=densityplotfile, format=ext[1:])
                densityplot.clf()
                if plt:
                    plt.close(densityplot)
            except Exception as error:
                self.ERROR("EDPluginControlSaxsAnalysisv1_0 in scatterplot: %s"%error)
            else:
                self.xsDataResult.densityPlot = XSDataFile(XSDataString(densityplotfile))
            gc.collect()
        self.synchronizePlugins()

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.postProcess")
        # Create some output data
        strLog = """Rg   =   %.2f +/- %2f
I(0) =   %.2e +/- %.2e
Points   %i to %i
Quality: %4.2f%%     Aggregated: %s""" % (self.autoRg.rg.value, self.autoRg.rgStdev.value,
                        self.autoRg.i0.value, self.autoRg.i0Stdev.value,
                        self.autoRg.firstPointUsed.value, self.autoRg.lastPointUsed.value,
                        self.autoRg.quality.value * 100., self.autoRg.isagregated.value)
        if self.gnom is None:
            strLog += """
datGnom failed"""
        else:
            strLog += """
Dmax    =    %12.2f       Total =   %12.2f
Guinier =    %12.2f       Gnom =    %12.2f""" % (self.gnom.dmax.value, self.gnom.total.value,
                        self.gnom.rgGuinier.value, self.gnom.rgGnom.value)
        if self.xVolume is None:
            strLog += """
datPorod failed"""
        else:
            strLog += """
Volume  =    %12.2f""" % (self.xVolume.value)

        self.xsDataResult.autoRg = self.autoRg
        self.xsDataResult.gnom = self.gnom
        self.xsDataResult.volume = self.xVolume
        self.xsDataResult.rti = self.rti
        self.xsDataResult.status = XSDataStatus(executiveSummary=XSDataString(strLog),
                                                message=self.getXSDataMessage())
        self.setDataOutput(self.xsDataResult)


    def doSuccessRg(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.doSuccessRg")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsAnalysisv1_0.doSuccessRg")
        self.retrieveMessages(_edPlugin)
        self.autoRg = _edPlugin.dataOutput.autoRgOut[0]

    def doFailureRg(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.doFailureRg")
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsAnalysisv1_0.doFailureRg")
        self.retrieveMessages(_edPlugin)
        self.setFailure()

    def doSuccessGnom(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.doSuccessGnom")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsAnalysisv1_0.doSuccessGnom")
        self.retrieveMessages(_edPlugin)
        self.gnom = _edPlugin.dataOutput.gnom
        self.gnomFile = self.gnom.gnomFile.path.value

    def doFailureGnom(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.doFailureGnom")
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsAnalysisv1_0.doFailureGnom")
        self.retrieveMessages(_edPlugin)
        #self.setFailure()

    def doSuccessPorod(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.doSuccessPorod")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsAnalysisv1_0.doSuccessPorod")
        self.retrieveMessages(_edPlugin)
        self.xVolume = _edPlugin.dataOutput.volume

    def doFailurePorod(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsAnalysisv1_0.doFailurePorod")
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsAnalysisv1_0.doFailurePorod")
        self.retrieveMessages(_edPlugin)
        #self.setFailure()

    def calculateRTI(self,autorg, scatterfile):
        if scatterfile is not None and\
            autorg.rg and autorg.rgStdev and autorg.i0.value > 0 and autorg.i0Stdev:
            #self.screen("Trying to calculate RTI")
            dictRTI = RamboTainerInvariant(numpy.loadtxt(scatterfile), autorg.rg.value,
                                           autorg.rgStdev.value, autorg.i0.value,
                                           autorg.i0Stdev.value, autorg.firstPointUsed.value)
            Vc = dictRTI.get("Vc")
            Vc_Stdev = dictRTI.get("dVc")
            Qr = dictRTI.get("Qr")
            Qr_Stdev = dictRTI.get("dQ")
            mass = dictRTI.get("mass")
            mass_Stdev = dictRTI.get("dmass")
           
        else:
            Vc = Vc_Stdev = Qr = Qr_Stdev = mass =  mass_Stdev = 0.0
        xsdRTI = XSDataRamboTainer(vc=XSDataDouble(Vc),
                                       qr=XSDataDouble(Qr),
                                       mass=XSDataDouble(mass),
                                       dvc=XSDataDouble(Vc_Stdev),
                                       dqr=XSDataDouble(Qr_Stdev),
                                       dmass=XSDataDouble(mass_Stdev))
        self.rti = xsdRTI