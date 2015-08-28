# coding: utf-8
#
#    Project: ExecPlugins/GroupAtsas
#             http://www.edna-site.org
#
#    Copyright (C) 2011, ESRF, Grenoble
#
#    Principal author:       Jérôme Kieffer
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

__author__ = "Jérôme Kieffer, Martha Brennich"
__license__ = "GPLv3+"
__copyright__ = "2014, ESRF, Grenoble"
__date__ = "28/08/2015"
__status__ = "production"

import os
from EDPluginExecProcessScript import EDPluginExecProcessScript
from XSDataEdnaSaxs import XSDataInputDatcmp, XSDataResultDatcmp, XSDataDouble


class EDPluginExecDatcmpv2_0(EDPluginExecProcessScript):
    """
    Execution plugin to run datcmp: evaluation of the Similarity of two curves (in a 2-3 column ascii file)
    Updated for ATSAS 2.6.0 CORMAP type test
    Currently only supports pair-wise test
    """


    def __init__(self):
        """
        Constructor of EDPluginExecDatcmpv2_0
        """
        EDPluginExecProcessScript.__init__(self)
        self.setXSDataInputClass(XSDataInputDatcmp)
        self.listFiles = []
        self.fChi = None
        self.fFidelity = None
        self.naFidelity = None
        self.atsasVersion = "2.6.1"
        self.testType = 'CHI-SQUARE' #CORMAP'
        self.outputType = 'FULL'
        self.commandString = ""


    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginExecDatcmpv2_0.checkParameters")
        self.checkMandatoryParameters(self.getDataInput(), "Data Input is None")
        self.checkMandatoryParameters(self.getDataInput().inputCurve, "No input 1D curves file provided")


    def preProcess(self, _edObject=None):
        EDPluginExecProcessScript.preProcess(self)
        self.DEBUG("EDPluginExecDatcmpv2_0.preProcess")
        self.listFiles = [i.path.value for i in self.getDataInput().inputCurve if os.path.isfile(i.path.value)]
        if len(self.listFiles) != 2:
            self.WARNING("You did not provide the right number of valid files !!! %s" % " ".join(self.listFiles))
        self.generateDatcmpScript()
        self.atsasVersion = self.config.get("atsasVersion", self.atsasVersion)

    def postProcess(self, _edObject=None):
        EDPluginExecProcessScript.postProcess(self)
        self.DEBUG("EDPluginExecDatcmpv2_0.postProcess")

        strResultFile = os.path.join(os.path.dirname(self.getScriptFilePath()), self.getScriptLogFileName())
        if os.path.isfile(strResultFile):
            for line in open(strResultFile):
                words = line.split()
                if (self.atsasVersion == "2.5.2") and (len(words) == 5):
                    try:
                        self.fChi = float(words[-2])
                        self.fFidelity = float(words[-1])
                    except ValueError:
                        self.WARNING("Strange ouptut from %s:%s %s" % (strResultFile, os.linesep, line))
                    else:
                        break
                elif (self.atsasVersion == "2.6.1") and (len(words) == 6) and ('vs.' in words):
                    if self.testType == 'CHI-SQUARE':
		                try:
		                    self.fChi = float(words[-3].strip('*')) ** 0.5
		                    self.fFidelity = float(words[-1].strip('*'))
		                except ValueError:
		                    self.WARNING("Strange ouptut from %s:%s %s" % (strResultFile, os.linesep, line))
		                else:
		                    break
					else:
		                try:
		                    self.fFidelity = float(words[-1].strip('*'))
		                    self.naFidelity = float(words[-2].strip('*'))
		                except ValueError:
		                    self.WARNING("Strange ouptut from %s:%s %s" % (strResultFile, os.linesep, line))
		                else:
		                    break

        # Create some output data
        xsDataResult = XSDataResultDatcmp()
        if self.fChi is not None:
            xsDataResult.chi = XSDataDouble(self.fChi)
        if self.fFidelity is not None:
            xsDataResult.fidelity = XSDataDouble(self.fFidelity)
        if self.naFidelity is not None:
            xsDataResult.nonadjustedFidelity = XSDataDouble(self.naFidelity)
        self.setDataOutput(xsDataResult)

    def generateDatcmpScript(self):
        self.DEBUG("EDPluginExecDatcmpv2_0.generateScript")
        self.commamdString = ""
        self.commandString += "--test=" + self.testType
        self.commandString += " --format="  +  self.outputType
        self.commandString += " " + " ".join(self.listFiles[:])
        self.setScriptCommandline(self.commandString)
        
