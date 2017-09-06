# coding: utf-8
#
#    Project: EdnaSaxs/Atsas
#             http://www.edna-site.org
#
#    File: "$Id$"
#
#    Copyright (C) 2012, ESRF Grenoble
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
__license__ = "MIT"
__copyright__ = "2012, ESRF Grenoble"
__date__ = "2017-07-26"
__status__ = "Development"

import os, shutil
from EDPluginExecProcessScript import EDPluginExecProcessScript
from XSDataEdnaSaxs import XSDataInputDatGnom, XSDataResultDatGnom, XSDataGnom
from XSDataCommon import XSDataString, XSDataDouble, XSDataFile, XSDataLength, XSDataStatus

class EDPluginExecDatGnomv1_1(EDPluginExecProcessScript):
    """
    Plugin that simply runs datgnom; a clever version of gnom; on a 1D curve. For ATSAS 2.8.2  
    """


    def __init__(self):
        """
        """
        EDPluginExecProcessScript.__init__(self)
        self.setXSDataInputClass(XSDataInputDatGnom)
        self.datFile = None
        self.outFile = None
        self.Rg = None
        self.skip = None

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginExecDatGnomv1_1.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.inputCurve, "No input Curve file provided")

    def preProcess(self, _edObject=None):
        EDPluginExecProcessScript.preProcess(self)
        self.DEBUG("EDPluginExecDatGnomv1_1.preProcess")
        self.datFile = self.dataInput.inputCurve.path.value
        if self.dataInput.output is not None:
            self.outFile = self.dataInput.output.path.value
        else:
            self.outFile = os.path.join(self.getWorkingDirectory(), os.path.splitext(os.path.basename(self.datFile)))
        if self.dataInput.rg is not None:
            self.rg = self.dataInput.rg.value
        if self.dataInput.skip is not None:
            self.skip = self.dataInput.skip.value
        self.symlink(self.datFile, os.path.basename(self.datFile))
        self.generateCommandLineOptions()


    def postProcess(self, _edObject=None):
        EDPluginExecProcessScript.postProcess(self)
        self.DEBUG("EDPluginExecDatGnomv1_1.postProcess")
        # Create some output data
        cwd = self.getWorkingDirectory()
        outfile = os.path.join(cwd, os.path.basename(self.outFile))
        if not os.path.isfile(outfile):
            self.error("EDPluginExecDatGnomv1_1 did not produce output file %s as expected !" % self.outFile)
            self.setFailure()
            self.dataOutput = XSDataResultDatGnom()
            return
        try:
            os.rename(outfile, self.outFile)
        except OSError:  # may fail if src and dst on different filesystem
            shutil.copy(outfile, self.outFile)
        gnom = XSDataGnom(gnomFile=XSDataFile(XSDataString(self.outFile)))
        logfile = os.path.join(self.getWorkingDirectory(), self.getScriptLogFileName())
        out = open(logfile, "r").read().split()
        for key, val, typ in (("dmax:", "dmax", XSDataLength),
                            ("Guinier:", "rgGuinier", XSDataLength),
                            ("Gnom:", "rgGnom", XSDataLength),
                            ("Total:", "total", XSDataDouble)):
            if key in out:
                idx = out.index(key)
                res = out[idx + 1]
                gnom.__setattr__(val, typ(float(res)))
            else:
                self.error("EDPluginExecDatGnomv1_1.postProcess No key %s in file %s" % (key, logfile))
                self.setFailure()
        self.dataOutput = XSDataResultDatGnom(gnom=gnom)
        self.dataOutput.status = XSDataStatus(message=self.getXSDataMessage())

    def generateCommandLineOptions(self):
        lstArg = [self.datFile, "-o", os.path.basename(self.outFile)]
        if self.rg:
            lstArg += ["-r", str(self.rg)]
        if self.skip and self.skip > 0:
            lstArg += ["--first", str(self.skip)] 
        self.setScriptCommandline(" ".join(lstArg))

    def symlink(self, filen, link):
        """
        Create a symlink to CWD with relative path
        """
        src = os.path.abspath(filen)
        cwd = self.getWorkingDirectory()
        dest = os.path.join(cwd, link)
        os.symlink(os.path.relpath(src, cwd), dest)
