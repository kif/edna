# coding: utf-8
#
#    Project: BioSaxs
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) 2012 ESRF
#
#    Principal author:        Al. de Maria
#                             Jerome Kieffer
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
__author__ = "Al. de Maria"
__license__ = "GPLv3+"
__copyright__ = "2012 ESRF"
__status__ = "Development"
__date__ = "20130124"

import os, shutil
from EDPluginControl        import EDPluginControl
from EDFactoryPlugin        import edFactoryPlugin
# from EDConfiguration        import EDConfiguration
from suds.client            import Client
from suds.transport.https   import HttpAuthenticated
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
from XSDataBioSaxsv1_0      import XSDataInputBioSaxsISPyBv1_0, XSDataResultBioSaxsISPyBv1_0
# , XSDataBioSaxsSample, XSDataGnom
# from XSDataEdnaSaxs         import XSDataAutoRg
from XSDataCommon           import  XSDataString, XSDataStatus
# XSDataInteger, XSDataDouble, XSDataString, XSDataFile, XSPluginItem, XSDataLength, XSDataBoolean, XSDataStatus


class EDPluginBioSaxsISPyBv1_0(EDPluginControl):

    URL = None
    CONF_URL_KEY = "url"
    CONF_URL_DEFAULT = "http://pcantolinos:8080/ispyb-ejb3/ispybWS/ToolsForBiosaxsWebService?wsdl"


    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsISPyBv1_0)
        self.dataBioSaxsSample = None
        self.dataAutoRg = None
        self.dataGnom = None

        # Params to be sent and I dont know them
        self.volume = None
        self.framesAverage = None
        self.framesMerged = None
        self.rg = None
        self.rgStdev = None
        self.i0 = None
        self.i0Stdev = None
        self.firstPointUsed = None
        self.lastPointUsed = None
        self.quality = None
        self.isagregated = None
        self.code = None
        self.concentration = None
        self.gnomFile = None
        self.rgGuinier = None
        self.rgGnom = None
        self.dmax = None
        self.total = None
        self.pyarchcurves = []
        self.pyarchgraph = {}
        self.lstError = []
        self.bestBuffer = None
        self.scatterPlot = None
        self.guinierPlot = None
        self.kratkyPlot = None
        self.densityPlot = None
        self.averageSample = None

        #For averages
        self.lstDiscarded = []
        self.lstAveraged = []
        self.averageFilePath = None
	self.visitorFilePath = None

        # For subtraction
        self.bufferFrames = []
        self.sampleFrames = [] 
        self.subtractedFilePath = None

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsISPyBv1_0.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.sample, "Sample is None")

    def configure(self):
        """
        Configures the ISPyB webservice access with the following parameters:
         - The "url" key from config file        
        """
        EDPluginControl.configure(self)
        if self.URL is None:
            self.DEBUG("EDPluginBioSaxsISPyBv1_0.configure")
            url = self.config.get(self.CONF_URL_KEY, None)
            if url:
                self.__class__.URL = url
            else:
                self.__class__.URL = self.CONF_URL_DEFAULT


    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)

        # Initializing webservices
        self.DEBUG("EDPluginBioSaxsISPyBv1_0.preProcess")
        self.dataBioSaxsSample = self.dataInput.sample
        user = None
        password = ""
        if self.dataBioSaxsSample:
            if self.dataBioSaxsSample.login:
                user = self.dataBioSaxsSample.login.value
                password = self.dataBioSaxsSample.passwd.value
                if self.dataBioSaxsSample.ispybURL:
                    self.URL = self.dataBioSaxsSample.ispybURL.value
        if not user:
            self.ERROR("No login/password information in sample configuration. Giving up.")
            self.setFailure()
            return

        if "http:" in  self.URL:
            self.URL = "https://ispyb.esrf.fr/ispyb/ispyb-ws/ispybWS/ToolsForBiosaxsWebService?wsdl"
        self.httpAuthenticatedToolsForBiosaxsWebService = HttpAuthenticated(username=user, password=password)
        self.client = Client(self.URL, transport=self.httpAuthenticatedToolsForBiosaxsWebService, cache=None)

        self.dataAutoRg = self.dataInput.autoRg
        self.dataGnom = self.dataInput.gnom
        if self.dataInput.bestBuffer is not None:
            self.bestBuffer = self.dataInput.bestBuffer.path.value

	if self.dataInput.averageFilePath is not None:
	    #Attention: averageFilePath is changed when copying data! visitorFilePath should stay constant
	    self.averageFilePath = self.dataInput.averageFilePath.path.value
            self.visitorFilePath = self.dataInput.averageFilePath.path.value

	if self.dataInput.sampleFrames is not None:
	    self.sampleFrames = self.dataInput.sampleFrames

	if self.dataInput.bufferFrames is not None:
	    self.bufferFrames = self.dataInput.bufferFrames
	if self.dataInput.subtractedFilePath is not None:
	    self.subtractedFilePath = self.dataInput.subtractedFilePath.path.value

        # Params to be sent and I dont know them
        if self.dataInput.volume:
            self.volume = self.dataInput.volume.value
        if self.dataInput.frameAverage:
            self.framesAverage = self.dataInput.frameAverage.value
        if self.dataInput.frameMerged:
            self.framesMerged = self.dataInput.frameMerged.value
	if self.dataInput.averageSample:
	    self.averageSample = self.dataInput.averageSample.path.value 
        if self.dataAutoRg:
            autoRg = self.dataAutoRg
            #if autoRg.filename and autoRg.filename.path:
            #    self.filename = autoRg.filename.path.value
            if autoRg.rg:
                self.rg = autoRg.rg.value
            if autoRg.rgStdev:
                self.rgStdev = autoRg.rgStdev.value
            if autoRg.i0:
                self.i0 = autoRg.i0.value
            if autoRg.i0Stdev:
                self.i0Stdev = autoRg.i0Stdev.value
            if autoRg.firstPointUsed:
                self.firstPointUsed = autoRg.firstPointUsed.value
            if autoRg.lastPointUsed:
                self.lastPointUsed = autoRg.lastPointUsed.value
            if autoRg.quality:
                self.quality = autoRg.quality.value
            if autoRg.isagregated:
                self.isagregated = autoRg.isagregated.value
        if self.dataBioSaxsSample:
            if self.dataBioSaxsSample.code:
                self.code = self.dataBioSaxsSample.code.value
            if self.dataBioSaxsSample.concentration:
                self.concentration = self.dataBioSaxsSample.concentration.value
        if self.dataGnom:
            if self.dataGnom.gnomFile:
                self.gnomFile = self.dataGnom.gnomFile.path.value
            if self.dataGnom.rgGuinier:
                self.rgGuinier = self.dataGnom.rgGuinier.value
            if self.dataGnom.rgGnom:
                self.rgGnom = self.dataGnom.rgGnom.value
            if self.dataGnom.dmax:
                self.dmax = self.dataGnom.dmax.value
            if self.dataGnom.total:
                self.total = self.dataGnom.total.value
        if self.dataInput.scatterPlot:
            self.scatterPlot = self.dataInput.scatterPlot.path.value
        if self.dataInput.guinierPlot:
            self.guinierPlot = self.dataInput.guinierPlot.path.value
        if self.dataInput.kratkyPlot:
            self.kratkyPlot = self.dataInput.kratkyPlot.path.value
        if self.dataInput.densityPlot:
            self.densityPlot = self.dataInput.densityPlot.path.value

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsISPyBv1_0.process")
        try:
            self.copy_to_pyarch()
        except Exception as error:
            strErrorMessage = "Error while copying to pyarch: %s" % error
            self.ERROR(strErrorMessage)
            self.lstError.append(strErrorMessage)
        if self.dataInput.sample.collectionOrder is not None:
            collectionOrder = str(self.dataInput.sample.collectionOrder.value)
        else:
            collectionOrder = "-1"

        try:
	    frames = []
	    #for i in range(0, len(self.averagedFrames)):
	    for filePath in self.lstAveraged: 
		frames.append({'filePath' : str(filePath)})

	    discarded = []
	    for filePath in self.lstDiscarded:              
	        discarded.append({'filePath' : str(filePath)})

            self.client.service.addAveraged(
                                    self.dataBioSaxsSample.measurementID.value,
	    			    collectionOrder,
                                    str(frames),
                                    str(discarded),
                                    str(self.averageFilePath),
				    str(self.visitorFilePath))
        except Exception, error:
            strError = "ISPyB error: %s" % error
            self.ERROR(strError)
            
        # Adding subtraction into ISPyB
        try:
	    
            if int(collectionOrder) == 2:
		sampleAvgOneDimensionalFiles = [] 
		bufferAvgOneDimensionalFiles = [] 

		for sampleFrameXSFile in self.sampleFrames:
			sampleAvgOneDimensionalFiles.append({'filePath': str(self.copyfile(sampleFrameXSFile.path.value, self.getPyarchFolder()))})

		for bufferFrameXSFile in self.bufferFrames:
			bufferAvgOneDimensionalFiles.append({'filePath': str(self.copyfile(bufferFrameXSFile.path.value, self.getPyarchFolder()))})

		self.subtractedFilePath = self.copyfile(self.subtractedFilePath, self.getPyarchFolder())

		self.client.service.addSubtraction(
						self.dataBioSaxsSample.measurementID.value,
                                                self.rg,
                                                self.rgStdev,
                			        self.i0,
		                                self.i0Stdev,
                		                self.firstPointUsed,
                		                self.lastPointUsed,
                		                self.quality,
                		                self.isagregated,
                                                self.rgGnom,
                                                self.dmax,
                                		self.total,
                		                self.volume,
                                                str(sampleAvgOneDimensionalFiles),
                                                str(bufferAvgOneDimensionalFiles),
                                                self.averageSample, 			#sampleAverageFilePath,
                                                self.bestBuffer,			#bufferAverageFilePath,
                                                self.subtractedFilePath, 		#subtractedFilePath,
                                                self.scatterPlot,  			#experimentalDataPlotFilePath,
                                                self.densityPlot, 			#densityPlotFilePath,
                                                self.guinierPlot, 			#guinierPlotFilePath,
                                                self.kratkyPlot,  			#kratkyPlotFilePath,
                                                self.gnomFile) 				#gnomOutputFilePath
                #self.client.service.storeDataAnalysisResultByMeasurementId(
                #                    self.dataBioSaxsSample.measurementID.value,
                #                    self.filename,
                #                    self.rg,
                #                    self.rgStdev,
                #                    self.i0,
                #                    self.i0Stdev,
                #                    self.firstPointUsed,
                #                    self.lastPointUsed,
                #                    self.quality,
                #                    self.isagregated,
                #                    self.code,
                #                    self.concentration,
                #                    self.gnomFile,
                #                    self.rgGuinier,
                #                    self.rgGnom,
                #                    self.dmax,
                #                    self.total,
                #                    self.volume,
                #                    self.framesAverage,
                #                    self.framesMerged,
                #                    ", ".join(self.pyarchcurves),
                #                    collectionOrder,
                #                    self.bestBuffer,
                #                    self.pyarchgraph.get("scatterPlot",""),
                #                    self.pyarchgraph.get("guinierPlot",""),
                #                    self.pyarchgraph.get("kratkyPlot",""),
                #                    self.pyarchgraph.get("densityPlot",""),
                #                    )
        except Exception, error:
            strError = "ISPyB error: %s" % error
            self.ERROR(strError)
            self.setFailure()


    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsISPyBv1_0.postProcess")
        # Create some output data

        xsDataResult = XSDataResultBioSaxsISPyBv1_0(status=XSDataStatus(executiveSummary=XSDataString(os.linesep.join(self.lstError))))
        self.setDataOutput(xsDataResult)
    
    def getPyarchFolder(self):
	return os.path.join(self.dataInput.sample.ispybDestination.path.value, "1d")
 
    def copy_to_pyarch(self):
        if self.dataInput.sample.ispybDestination:
            pyarch = self.getPyarchFolder()
            try:
                if not os.path.isdir(pyarch):
                    os.makedirs(pyarch)
            except IOError as error:
                ermsg = "Error while directory creation in pyarch: %s " % error
                self.lstError.append(ermsg)
                self.ERROR(ermsg)
            #for xsdfile in self.dataInput.curves:
            for xsdfile in (self.dataInput.discardedFrames):
                if xsdfile:
                    self.lstDiscarded.append(self.copyfile(xsdfile.path.value, pyarch))

            for xsdfile in (self.dataInput.averagedFrames):
	        if xsdfile:
	            self.lstAveraged.append(self.copyfile(xsdfile.path.value, pyarch))

	    # Used only for the average
	    self.averageFilePath = self.copyfile(self.averageFilePath,pyarch)

            # Used for the subtraction	    
            self.bestBuffer = self.copyfile(self.bestBuffer, pyarch)
            self.gnomFile = self.copyfile(self.gnomFile, pyarch)
            self.scatterPlot = self.copyfile(self.scatterPlot, pyarch)
            self.guinierPlot = self.copyfile(self.guinierPlot, pyarch)
            self.kratkyPlot = self.copyfile(self.kratkyPlot, pyarch)
            self.densityPlot = self.copyfile(self.densityPlot, pyarch)
	    self.averageSample = self.copyfile(self.averageSample, pyarch)
	

    def copyfile(self, afile, pyarch):
        if not pyarch:
            self.ERROR("pyArch is %s" % pyarch)
        if afile and os.path.exists(afile) and os.path.isdir(pyarch):
            try:
                shutil.copy(afile, pyarch)
		return os.path.join(pyarch, os.path.basename(afile))
            except IOError as error:
                ermsg = "BioSAXSISPyBv1_0: Error while copying %s to pyarch %s: %s " % (afile, pyarch, error)
                self.lstError.append(ermsg)
                self.WARNING(ermsg)
	return None
