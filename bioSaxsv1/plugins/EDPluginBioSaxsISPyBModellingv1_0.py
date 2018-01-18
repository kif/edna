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

import os, shutil, json
from EDPluginControl        import EDPluginControl
from EDFactoryPlugin        import edFactoryPlugin
from suds.client            import Client
from suds.transport.https   import HttpAuthenticated
edFactoryPlugin.loadModule("XSDataBioSaxsv1_0")
from XSDataBioSaxsv1_0      import XSDataInputBioSaxsISPyBModellingv1_0, XSDataResultBioSaxsISPyBModellingv1_0
from XSDataCommon           import  XSDataString, XSDataStatus

def sensibleDict(dico):
    """
    Convert an EDNA XSData dictionary into a normal python dictionary for human beings
    
    This is likely to be called recursively 
    """
    sdic = {}
    if not isinstance(dico, dict):
        return dico
    mytype = dico.get("__XSDataName", "")
    if mytype == 'XSDataCommon.XSDataFile':
        return  dico.get("_path", {}).get("_value", "")
    elif mytype == 'XSDataCommon.XSDataDouble':
        return  dico.get("_value", "")
    elif mytype == 'XSDataCommon.XSDataString':
        return  dico.get("_value", "")

    else:
        for key in dico:
            if key.startswith("__"):
                continue
            if key.startswith("_"):
                newkey = key[1:]
            else:
                newkey = key
            val = dico[key]
            if isinstance(val, dict):
                sdic[newkey] = sensibleDict(val) or {}
            elif isinstance(val, (list, tuple)):
                sdic[newkey] = [sensibleDict(i) for i in val]
            else:
                sdic[newkey] = dico[key]
        return sdic


class EDPluginBioSaxsISPyBModellingv1_0(EDPluginControl):

    URL = None
    CONF_URL_KEY = "url"
    CONF_URL_DEFAULT = "http://pcantolinos:8080/ispyb-ejb3/ispybWS/ToolsForBiosaxsWebService?wsdl"


    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputBioSaxsISPyBModellingv1_0)

        # Params to be sent and I dont know them
#         self.modellingResult = None
        self.models = []
        self.damaver = {}
        self.damfilt = {}
        self.dammin = {}
        self.nsdPlot = ""
        self.chi2plot = ""
        self.pyarchmodels = []
        self.lstError = []

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginBioSaxsISPyBModellingv1_0.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.sample, "Sample is None")
#         self.checkMandatoryParameters(self.dataInput.saxsModelingResult, "SaxsModelingResult is None")


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
        self.DEBUG("EDPluginBioSaxsISPyBModellingv1_0.preProcess")
        self.dataBioSaxsSample = self.dataInput.sample


        user = None
        password = ""
        if self.dataBioSaxsSample:
            if self.dataBioSaxsSample.login:
                user = self.dataBioSaxsSample.login.value
                password = self.dataBioSaxsSample.passwd.value
                self.URL = self.dataBioSaxsSample.ispybURL
        if not user:
            self.ERROR("No login/password information in sample configuration. Giving up.")
            self.setFailure()
            return
#         self.modellingResult = self.dataInput.saxsModelingResult
        if "http:" in  self.URL:
            #print "OLD URL"
            self.URL = self.URL.replace("http:", "https:",1)
       
            self.URL = self.URL.replace("8080","",1)
        # I don't trust in this authentication.... but it is going to work soon
        self.httpAuthenticatedToolsForBiosaxsWebService = HttpAuthenticated(username=user, password=password)
        self.client = Client(self.dataBioSaxsSample.ispybURL.value, transport=self.httpAuthenticatedToolsForBiosaxsWebService, cache=None)


    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginBioSaxsISPyBModellingv1_0.process")
        self.addExecutiveSummaryLine("Registering analysis in ISPyB")
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
        dico = sensibleDict(self.dataInput.exportToDict())
        self.models = dico.get("dammifModels", []) or []
        self.damaver = dico.get("damaverModel", {}) or {}
        self.damfilt = dico.get("damfiltModel", {}) or {}
        self.dammin = dico.get("damminModel", {}) or {}
        self.nsdPlot = dico.get("nsdPlot", "") or ""
        self.chi2plot = dico.get("chiRfactorPlot", "") or ""
        try:
            self.copy_to_pyarch()
        except Exception as error:
            strErrorMessage = "Error while copying to pyarch: %s" % error
            self.ERROR(strErrorMessage)
            self.lstError.append(strErrorMessage)
            self.writeErrorTrace()
        try:
            self.id = self.dataInput.sample.measurementID.value
            return self.client.service.storeAbInitioModels("[%s]" % self.id ,
                                                           json.dumps(self.models),
                                                           json.dumps(self.damaver),
                                                           json.dumps(self.damfilt),
                                                           json.dumps(self.dammin),
                                                           self.nsdPlot,
                                                           self.chi2plot)
        except Exception, error:
            strError = "ISPyB error: %s" % error
            self.ERROR(strError)
            self.setFailure()


    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginBioSaxsISPyBModellingv1_0.postProcess")
        # Create some output data

        xsDataResult = XSDataResultBioSaxsISPyBModellingv1_0(status=XSDataStatus(executiveSummary=XSDataString(os.linesep.join(self.lstError))))
        self.setDataOutput(xsDataResult)

    def copyCustomFile(self, model, key, targetFileName):
        pyarch = os.path.join(self.dataInput.sample.ispybDestination.path.value, str(self.dataInput.sample.measurementID.value))
        if model.get(key) is not None:
            afile = model.get(key)
            model[key] = self.copyfile(afile, pyarch, targetFileName)
    
    def copyModelFile(self, model, i):
    	#afile = model.get("pdbFile")I
    	#if afile:
	#    amodel["pdbFile"] = self.copyfile(afile, pyarch, "model_%02i.pdb" % i)
        self.copyCustomFile(model, "pdbFile", "model_%02i.pdb" % i)
        self.copyCustomFile(model, "fitFile", "model_%02i.fit" % i)
        self.copyCustomFile(model, "firFile", "model_%02i.fir" % i)
        self.copyCustomFile(model, "logFile", "model_%02i.log" % i)
	

    def copy_to_pyarch(self):
        if self.dataInput.sample.ispybDestination:
            pyarch = os.path.join(self.dataInput.sample.ispybDestination.path.value, str(self.dataInput.sample.measurementID.value))
            pyarch = self.updatePyArch(pyarch)
            try:
                if not os.path.isdir(pyarch):
                    os.makedirs(pyarch)
            except IOError as error:
                ermsg = "Error while directory creation in pyarch: %s " % error
                self.lstError.append(ermsg)
                self.ERROR(ermsg)

            for i in range(0, len(self.models)):
	    	    self.copyModelFile(self.models[i], i)
		


            self.damaver["pdbFile"] = self.copyfile(self.damaver.get("pdbFile"), pyarch, "damaver.pdb")
            self.damfilt["pdbFile"] = self.copyfile(self.damfilt.get("pdbFile"), pyarch, "damfilt.pdb")
            self.dammin["pdbFile"] = self.copyfile(self.dammin.get("pdbFile"), pyarch, "dammin.pdb")
            self.dammin["fitFile"] = self.copyfile(self.dammin.get("fitFile"), pyarch, "dammin.fit")
            self.dammin["firFile"] = self.copyfile(self.dammin.get("firFile"), pyarch, "dammin.fir")
            self.dammin["logFile"] = self.copyfile(self.dammin.get("logFile"), pyarch, "dammin.log")
            self.nsdPlot = self.copyfile(self.nsdPlot, pyarch)
            self.chi2plot = self.copyfile(self.chi2plot, pyarch)




    def copyfile(self, afile, pyarch, dest=None):
        fullname = None
        if not pyarch:
            self.ERROR("pyArch is %s" % pyarch)
            return

        if afile and os.path.exists(afile) and os.path.isdir(pyarch):
            if dest == "model":
                fullname = os.path.join(pyarch, "model_%02i.pdb" % (len(self.pyarchmodels) + 1))
                self.pyarchmodels.append(fullname)
            elif dest:
                fullname = os.path.join(pyarch, dest)
            else:
                fullname = os.path.join(pyarch, os.path.basename(afile))
            try:
                shutil.copy(afile, fullname)
            except IOError as error:
                ermsg = "Error while copying %s to pyarch %s: %s " % (afile, pyarch, error)
                self.lstError.append(ermsg)
                self.WARNING(ermsg)
        return fullname
    
    def updatePyArch(self,originalPyArch):
        now = datetime.datetime.now()
        yearString = str(now.year)
        archDatePattern = re.compile(r'/pyarch/20[0-9][0-9]')
        newPyArch = archDatePattern.sub("/pyarch/" + yearString,originalPyArch)
        return newPyArch
