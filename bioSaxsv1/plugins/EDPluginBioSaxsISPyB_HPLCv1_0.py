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
__date__ = "20131003"

import os, shutil, traceback
from EDPluginControl        import EDPluginControl
from EDFactoryPlugin        import edFactoryPlugin
from suds.client            import Client
from suds.transport.http    import HttpAuthenticated
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
from XSDataBioSaxsv1_0      import XSDataResultBioSaxsISPyB_HPLCv1_0, XSDataInteger
from XSDataBioSaxsv1_0      import XSDataInputBioSaxsISPyB_HPLCv1_0
from XSDataCommon           import  XSDataString, XSDataStatus


class EDPluginBioSaxsISPyB_HPLCv1_0(EDPluginControl):

    URL = None
    CONF_URL_KEY = "url"
    CONF_URL_DEFAULT = "http://ispyvalid.esrf.fr:8080/ispyb-ejb3/ispybWS/ToolsForBiosaxsWebService?wsdl"


    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsISPyB_HPLCv1_0)
        self.dataBioSaxsSample = None
        self.pyarchcurves = []
        self.pyarchgraph = {}
        self.lstError = []
        self.code = "opd"
        self.number = "29"
        self.hdf5File = None
        self.jsonFile = None
        self.hplcPlot = None
        self.visitorFilePath = None

        self.xsdResult = XSDataResultBioSaxsISPyB_HPLCv1_0()


    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsISPyB_HPLCv1_0.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.sample, "Sample is None")

    def configure(self):
        """
        Configures the ISPyB_HPLC webservice access with the following parameters:
         - The "url" key from config file        
        """
        EDPluginControl.configure(self)
        if self.URL is None:
            self.DEBUG("EDPluginBioSaxsISPyB_HPLCv1_0.configure")
            url = self.config.get(self.CONF_URL_KEY, None)
            if url:
                self.__class__.URL = url
            else:
                self.__class__.URL = self.CONF_URL_DEFAULT


    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        # Initializing webservices
        self.DEBUG("EDPluginBioSaxsISPyB_HPLCv1_0.preProcess")
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
            self.ERROR("No login/password information in sample configuration. Giving up. falling back on test mode")
#             self.setFailure()
#             return
            user = "opd29"
            password = "tonic29"

        # construct code based on letters + numbers user concept
        self.code = ''.join(i for i in user if i.isalpha())
        self.number = ''.join(i for i in user if i.isdigit())


        # I don't trust in this authentication.... but it is going to work soon

        self.httpAuthenticatedToolsForBiosaxsWebService = HttpAuthenticated(username=user, password=password)
        self.client = Client(self.URL, transport=self.httpAuthenticatedToolsForBiosaxsWebService, cache=None)
        
        
        if (self.dataInput.hdf5File is not None):
            if (self.dataInput.hdf5File.path is not None):
                #Attenttion: self.hdf5File will be changed later on
                self.hdf5File = self.dataInput.hdf5File.path.value
                self.visitorFilePath = self.dataInput.hdf5File.path.value
               
        if (self.dataInput.jsonFile is not None):
            if (self.dataInput.jsonFile.path is not None):        
                self.jsonFile = self.dataInput.jsonFile.path.value
                
        if (self.dataInput.hplcPlot is not None):
            if (self.dataInput.hplcPlot.path is not None):   
                self.hplcPlot = self.dataInput.hplcPlot.path.value    

    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsISPyB_HPLCv1_0.process")
        self.experimentId = None
        try:
            experimentId = self.client.service.createHPLC(
                                                            self.code,
                                                            self.number,
                                                            'name')
        except Exception, error:
            strError = "ISPyB error: %s" % error
            self.ERROR(strError)
            self.setFailure()
        try:
            self.copy_to_pyarch(str(experimentId))
            
        except Exception as error:
            strErrorMessage = "In EDPluginBioSaxsISPyB_HPLCv1_0.process: Error while copying to pyarch: %s" % error
            self.ERROR(strErrorMessage)
            self.lstError.append(strErrorMessage)
            
        try:
            if (experimentId is not None):
                self.xsdResult.setExperimentId(XSDataInteger(experimentId))
                #Updating folder
                if (self.dataInput.sample.ispybDestination is not None):
                    if (self.dataInput.sample.ispybDestination.path is not None):
                        self.dataInput.sample.ispybDestination.path.value = os.path.join(self.dataInput.sample.ispybDestination.path.value, str(experimentId))
                        self.xsdResult.setSample(self.dataInput.sample)
                        self.client.service.storeHPLC(
                                            experimentId,
                                            self.hdf5File,
                                            self.jsonFile, self.visitorFilePath)
        except Exception as error:
            traceback.print_exc()
            strErrorMessage = "ISPyB storeHPLC error: %s" % error
            self.ERROR(strErrorMessage)
            self.lstError.append(strErrorMessage)


    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsISPyB_HPLCv1_0.postProcess")
        # Create some output data

        xsDataResult = XSDataResultBioSaxsISPyB_HPLCv1_0(status=XSDataStatus(executiveSummary=XSDataString(os.linesep.join(self.lstError))))
        self.setDataOutput(xsDataResult)

    def copy_to_pyarch(self, experimentId):
        if self.dataInput.sample.ispybDestination:
            pyarch = os.path.join(self.dataInput.sample.ispybDestination.path.value)
            if (experimentId is not None):
                pyarch = os.path.join(self.dataInput.sample.ispybDestination.path.value, experimentId)
            try:
                if not os.path.isdir(pyarch):
                    os.makedirs(pyarch)
            except IOError as error:
                ermsg = "Error while directory creation in pyarch: %s " % error
                self.lstError.append(ermsg)
                self.ERROR(ermsg)
            self.hdf5File = self.copyfile(self.hdf5File, pyarch)
            self.jsonFile = self.copyfile(self.jsonFile, pyarch)
            self.hplcPlot = self.copyfile(self.hplcPlot, pyarch)

    def copyfile(self, afile, pyarch, dest="curve"):
        if not pyarch:
            self.ERROR("pyArch is %s" % pyarch)
        if afile and os.path.exists(afile) and os.path.isdir(pyarch):
            try:
                shutil.copy(afile, pyarch)
            except IOError as error:
                ermsg = "BiosaxsISPyB HPLC Error while copying %s to pyarch %s: %s " % (afile, pyarch, error)
                self.lstError.append(ermsg)
                self.WARNING(ermsg)
            else:
                return os.path.join(pyarch, os.path.basename(afile))  # os.path.join(pyarch, afile)
