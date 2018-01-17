# coding: utf-8
#
#    Project: BioSaxs PROJECT
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
__status__ = "development"

import os, gc
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from EDThreading import Semaphore
from EDPluginControl import EDPluginControl
from XSDataCommon import XSDataStatus, XSDataString, XSDataBoolean, XSDataInteger, XSDataFile
from XSDataEdnaSaxs import XSDataInputSaxsModeling, XSDataResultSaxsModeling, \
                            XSDataInputDammif, XSDataInputDamaver, \
                            XSDataInputDamstart, XSDataInputDamfilt, XSDataInputDammin

from freesas.align import AlignModels

class EDPluginControlSaxsModelingv1_1(EDPluginControl):
    """
    Basically this is a re-implementation of EDPluginControlSolutionScattering starting after Gnom and without web page generation
    
    new in 1.1: replace supcomb with FreeSAS implementation
    """
    classlock = Semaphore()
    configured = False
    cluster_size = 2  # duplicate from ControlPlugin
    dammif_jobs = 16  # number of dammif job to run
    unit = "NANOMETER"  # unit of the GNOM file
    symmetry = "P1"  #
    mode = "fast"  #
    # constants:  plugin names
    strPluginExecDammif = "EDPluginExecDammifv0_2"
    strPluginExecSupcomb = "EDPluginExecSupcombv0_3"
    strPluginExecDamaver = "EDPluginExecDamaverv0_3"
    strPluginExecDamfilt = "EDPluginExecDamfiltv0_3"
    strPluginExecDamstart = "EDPluginExecDamstartv0_3"
    strPluginExecDammin = "EDPluginExecDamminv0_2"
    Rg_min = 0.5  # nm
    def __init__(self):
        """
        """
        EDPluginControl.__init__(self)
        self.setXSDataInputClass(XSDataInputSaxsModeling)
        self.edPlugin = None
        self.edPlugin = None
        self.xsGnomFile = None
        self.result = XSDataResultSaxsModeling()
        self.result.dammifModels = []
        self.graph_format = "png"
        self.dammif_plugins = []
        self.dammif = None
        self.supcomb_plugins = {}
        self.actclust_supcomb = None
        self.valid = None  # index of valid damif models
        self.mask2d = None
        self.arrayNSD = None
        self.ref = None  # reference frame number (starting ar 0)

    def checkParameters(self):
        """
        Checks the mandatory parameters.
        """
        self.DEBUG("EDPluginControlSaxsModelingv1_1.checkParameters")
        self.checkMandatoryParameters(self.dataInput, "Data Input is None")
        self.checkMandatoryParameters(self.dataInput.gnomFile, "gnom output is missing")

    def configure(self):
        EDPluginControl.configure(self)
        if not self.configured:
            with self.classlock:
                if not self.configured:
                    EDPluginControl.configure(self)
                    dammif_jobs = self.config.get("dammifJobs", None)
                    if (dammif_jobs != None):
                        self.__class__.dammif_jobs = int(dammif_jobs)
                        self.DEBUG("EDPluginControlSaxsModelingv1_1.configure: setting number of dammif jobs to %d" % self.dammif_jobs)
                    unit = self.config.get("unit", None)
                    if (unit != None):
                        self.__class__.unit = unit.upper()
                        self.DEBUG("EDPluginControlSaxsModelingv1_1.configure: setting input units to %s" % self.unit)
                    symmetry = self.config.get("symmetry", None)
                    if (symmetry != None):
                        self.__class__.symmetry = symmetry
                        self.DEBUG("EDPluginControlSaxsModelingv1_1.configure: setting symmetry to %s" % self.symmetry)
                    mode = self.config.get("mode", None)
                    if (mode != None):
                        self.__class__.mode = mode
                        self.DEBUG("EDPluginControlSaxsModelingv1_1.configure: setting dammif mode to %s" % self.mode)
                    clusterSize = self.config.get("clusterSize", None)
                    if (clusterSize != None):
                        self.__class__.cluster_size = int(clusterSize)
                        self.DEBUG("EDPluginControl.configure: setting cluster size to %d" % self.cluster_size)
                    self.__class__.configured = True


    def preProcess(self, _edObject=None):
        EDPluginControl.preProcess(self)
        self.DEBUG("EDPluginControlSaxsModelingv1_1.preProcess")
        self.xsGnomFile = self.dataInput.gnomFile
        if self.dataInput.graphFormat:
            self.graph_format = self.dataInput.graphFormat.value
        self.checkRg()


    def checkRg(self):
        """
        If there is nothing in the sample, Rg = 0.1 nm 
        damaver is likely to produce log files of many GB  
        """
        last_line = open(self.xsGnomFile.path.value).readlines()[-1]
#         self.WARNING("last Gnom file line is %s" % last_line)
        key = "Rg ="
        start = last_line.find(key) + len(key)
        val = last_line[start:].split()[0]
        try:
            rg = float(val)
        except ValueError:
            rg = 0.0
        if rg < self.Rg_min:
            str_err = "Radius of Giration is too small (%s<%s). Stop processing !!!!" % (rg, self.Rg_min)
            self.ERROR(str_err)
            self.setFailure()
            raise RuntimeError(str_err)
 
    def process(self, _edObject=None):
        EDPluginControl.process(self)
        self.DEBUG("EDPluginControlSaxsModelingv1_1.process")
        xsDataInputDammif = XSDataInputDammif(gnomOutputFile=self.xsGnomFile,
                                              unit=XSDataString(self.unit),
                                              symmetry=XSDataString(self.symmetry),
                                              mode=XSDataString(self.mode))
        for i in range(self.dammif_jobs):
            dammif = self.loadPlugin(self.strPluginExecDammif)
            dammif.connectSUCCESS(self.doSuccessExecDammif)
            dammif.connectFAILURE(self.doFailureExecDammif)
            xsd = xsDataInputDammif.copyViaDict()
            xsd.order = XSDataInteger(i + 1)
            dammif.dataInput = xsd
            self.addPluginToActionCluster(dammif)
            self.dammif_plugins.append(dammif)
        self.executeActionCluster()
        self.synchronizeActionCluster()
        for plugin in self.dammif_plugins:
            if plugin.isFailure():
                self.ERROR("dammif plugin %s-%08i failed" % (plugin.getName(), plugin.getId()))
                self.setFailure()
            self.retrieveMessages(plugin)
        if self.isFailure():
            return

        # retrieve results from best dammif
        self.dammif = self.bestDammif()

        self.chi2plot("chi2_R.png")
        self.result.chiRfactorPlot = XSDataFile(XSDataString(os.path.join(self.getWorkingDirectory(), "chi2_R.png")))

        # temporary results: use best dammif
        self.result.fitFile = self.dammif.dataOutput.fitFile
        self.result.logFile = self.dammif.dataOutput.logFile
        self.result.pdbMoleculeFile = self.dammif.dataOutput.pdbMoleculeFile
        self.result.pdbSolventFile = self.dammif.dataOutput.pdbSolventFile

        # align models, compute NSD and choose the reference model
        inputfiles = [self.dammif_plugins[idx].dataOutput.pdbMoleculeFile.path.value for idx in range(self.dammif_jobs)]
        align = AlignModels(inputfiles, slow=False)
        
        outputfiles = []
        for i in range(self.dammif_jobs):
            outputfiles.append(os.path.join(self.getWorkingDirectory(), "model-%02i.pdb" % (i+1)))
        align.outputfiles = outputfiles
        #align.outputfiles = ["model-%02i.pdb" % (i+1) for i in range(self.dammif_jobs)]
        align.validmodels = self.valid
        align.assign_models()
        align.makeNSDarray()
        align.alignment_reference()
        self.ref = align.reference
        
        pngfile = os.path.join(self.getWorkingDirectory(), "nsd.png")
        align.plotNSDarray(filename=pngfile ,save=True)
        self.result.nsdPlot = XSDataFile(XSDataString(pngfile))

#        Now that all (valid) models are aligned we can combine them using damaver
        pdbFiles = [XSDataFile(XSDataString(align.outputfiles[self.ref]))]

        for idx in range(self.dammif_jobs):
            if self.valid[idx] and idx != self.ref:
                pdbFiles.append(XSDataFile(XSDataString(align.outputfiles[idx])))

        damaver = self.loadPlugin(self.strPluginExecDamaver)
        damaver.dataInput = XSDataInputDamaver(pdbInputFiles=pdbFiles,
                                                automatic=XSDataBoolean(False))
        damaver.connectSUCCESS(self.doSuccessExecDamaver)
        damaver.connectFAILURE(self.doFailureExecDamaver)
        damaver.executeSynchronous()

        if self.isFailure():
            return

        damfilt = self.loadPlugin(self.strPluginExecDamfilt)
        damfilt.dataInput = XSDataInputDamfilt(inputPdbFile=damaver.dataOutput.damaverPdbFile)
        damfilt.connectSUCCESS(self.doSuccessExecDamfilt)
        damfilt.connectFAILURE(self.doFailureExecDamfilt)
        damfilt.execute()
        ########################################################################
        # TODO: This is a dead end : do it in parallel
        ########################################################################

        if self.isFailure():
            return

        damstart = self.loadPlugin(self.strPluginExecDamstart)
        damstart.dataInput = XSDataInputDamstart(inputPdbFile=damaver.dataOutput.damaverPdbFile)
        damstart.connectSUCCESS(self.doSuccessExecDamstart)
        damstart.connectFAILURE(self.doFailureExecDamstart)
        damstart.executeSynchronous()

        if self.isFailure():
            return
        ########################################################################
        # Finally call dammin
        ########################################################################
        if self.config.get("do_dammin") in ["False", "0", False, 0]:
            pass
            #return
        dammin = self.loadPlugin(self.strPluginExecDammin)
        dammin.dataInput = XSDataInputDammin(pdbInputFile=damstart.dataOutput.outputPdbFile,
                                             gnomOutputFile=self.xsGnomFile,
                                             symmetry=XSDataString(self.symmetry),
                                             mode=XSDataString(self.mode))
        dammin.connectSUCCESS(self.doSuccessExecDammin)
        dammin.connectFAILURE(self.doFailureExecDammin)
        dammin.executeSynchronous()
        # Dammin takes as lot of time ... wait here for completion

    def postProcess(self, _edObject=None):
        EDPluginControl.postProcess(self)
        self.DEBUG("EDPluginControlSaxsModelingv1_1.postProcess")

        self.synchronizePlugins()
        # Create some output data
#        self.result.

    def finallyProcess(self, _edObject=None):
        EDPluginControl.finallyProcess(self, _edObject=_edObject)
        self.result.status = XSDataStatus(message=self.getXSDataMessage(),
                                          executiveSummary=XSDataString(os.linesep.join(self.getListExecutiveSummaryLines())))

        self.setDataOutput(self.result)
        # clean up memory
        self.dammif_plugins = []
        self.dammif = None
        self.emptyListOfLoadedPlugin()
        self.supcomb_plugins = {}
        self.actclust_supcomb = None
        gc.collect()

    def doSuccessExecDammif(self, _edPlugin=None):
        """
        Locked as dammif is called many times in parallel
        """
        with self.locked():
            self.DEBUG("EDPluginControlSaxsModelingv1_1.doSuccessExecDammif")
            self.retrieveMessages(_edPlugin)
            self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDammif")
            try:
                self.result.dammifModels.append(_edPlugin.dataOutput.model)
                # this has to be done only for the best model (once determined) !
#                self.result.pdbMoleculeFile = _edPlugin.dataOutput.pdbMoleculeFile
#                self.result.pdbSolventFile = _edPlugin.dataOutput.pdbSolventFile
#                self.result.fitFile = _edPlugin.dataOutput.fitFile
#                self.result.logFile = _edPlugin.dataOutput.logFile
            except Exception as error:
                self.ERROR("Error in doSuccessExecDammif: %s" % error)

    def doFailureExecDammif(self, _edPlugin=None):
        """
        Locked as dammif is called many times in parallel
        """

        with self.locked():
            self.DEBUG("EDPluginControlSaxsModelingv1_1.doFailureExecDammif")
            self.retrieveMessages(_edPlugin)
            self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDammif")
            self.setFailure()


    def doSuccessExecDamaver(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doSuccessExecDamaver")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doSuccessExecDamaver")
        self.retrieveMessages(_edPlugin)
        try:
            self.result.damaverModel = _edPlugin.dataOutput.model
            self.symlink(_edPlugin.dataOutput.model.pdbFile.path.value, _edPlugin.dataOutput.model.name.value + ".pdb")
        except Exception as error:
            self.ERROR("Error in doSuccessExecDamaver: %s" % error)

    def doFailureExecDamaver(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doFailureExecDamaver")
        self.retrieveMessages(_edPlugin)
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDamaver")
        self.setFailure()


    def doSuccessExecDamfilt(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doSuccessExecDamfilt")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doSuccessExecDamfilt")
        self.retrieveMessages(_edPlugin)
        try:
            self.result.damfiltModel = _edPlugin.dataOutput.model
            self.symlink(_edPlugin.dataOutput.model.pdbFile.path.value, _edPlugin.dataOutput.model.name.value + ".pdb")
        except Exception as error:
            self.ERROR("Error in doSuccessExecDamfilt: %s" % error)



    def doFailureExecDamfilt(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doFailureExecDamfilt")
        self.retrieveMessages(_edPlugin)
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDamfilt")
        self.setFailure()


    def doSuccessExecDamstart(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doSuccessExecDamstart")
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doSuccessExecDamstart")
        self.retrieveMessages(_edPlugin)
        try:
            self.result.damstartModel = _edPlugin.dataOutput.model
            self.symlink(_edPlugin.dataOutput.model.pdbFile.path.value, _edPlugin.dataOutput.model.name.value + ".pdb")
        except Exception as error:
            self.ERROR("Error in doSuccessExecDamstart: %s" % error)


    def doFailureExecDamstart(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doFailureExecDamstart")
        self.retrieveMessages(_edPlugin)
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDamstart")
        self.setFailure()


    def doSuccessExecDammin(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doSuccessExecDammin")
        self.retrieveMessages(_edPlugin)
        self.retrieveSuccessMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDammin")
        try:
            self.result.pdbMoleculeFile = _edPlugin.dataOutput.pdbMoleculeFile
            self.result.pdbSolventFile = _edPlugin.dataOutput.pdbSolventFile
            self.result.fitFile = _edPlugin.dataOutput.fitFile
            self.result.firFile = _edPlugin.dataOutput.model.firFile
            self.result.logFile = _edPlugin.dataOutput.logFile
            self.result.damminModel = _edPlugin.dataOutput.model
            self.symlink(_edPlugin.dataOutput.model.pdbFile.path.value, _edPlugin.dataOutput.model.name.value + ".pdb")
        except Exception as error:
            self.ERROR("Error in doSuccessExecDammin: %s" % error)


    def doFailureExecDammin(self, _edPlugin=None):
        self.DEBUG("EDPluginControlSaxsModelingv1_1.doFailureExecDammin")
        self.retrieveMessages(_edPlugin)
        self.retrieveFailureMessages(_edPlugin, "EDPluginControlSaxsModelingv1_1.doFailureExecDammin")
        self.setFailure()


    def bestDammif(self):
        """
        Find DAMMIF run with best chi-square value
        """
        fitResultDict = dict([(plg.dataOutput.chiSqrt.value, plg)
                              for plg in self.dammif_plugins
                              if (plg.dataOutput is not None) and (plg.dataOutput.chiSqrt is not None)])
        fitResultList = fitResultDict.keys()
        fitResultList.sort()

        return fitResultDict[fitResultList[0]]

    def symlink(self, filen, link):
        """
        Create a symlink to CWD with relative path
        """
        src = os.path.abspath(filen)
        cwd = self.getWorkingDirectory()
        dest = os.path.join(cwd, link)
        os.symlink(os.path.relpath(src, cwd), dest)


    def chi2plot(self, filename=None, close=True):

        chi2 = numpy.array([ plg.dataOutput.chiSqrt.value if plg.dataOutput is not None else numpy.NaN for plg in self.dammif_plugins])
        chi2max = chi2.mean() + 2 * chi2.std()

        xticks = 1 + numpy.arange(self.dammif_jobs)
        fig = plt.figure(figsize=(15, 10))
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.bar(xticks - 0.5, chi2)
        ax1.set_ylabel(u"$\sqrt{\u03C7}$")
        ax1.set_xlabel(u"Model number")
        ax1.plot([0.5, self.dammif_jobs + 0.5], [chi2max, chi2max], "-r", label=u"${\u03C7}^2$$_{max}$ = %.3f" % chi2max)
        ax1.set_xticks(xticks)
        ax1.legend(loc=8)
        R = numpy.array([ plg.dataOutput.rfactor.value if plg.dataOutput is not None else numpy.NaN for plg in self.dammif_plugins])
        Rmax = R.mean() + 2 * R.std()
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.bar(xticks - 0.5, R)
        ax2.plot([0.5, self.dammif_jobs + 0.5], [Rmax, Rmax], "-r", label=u"R$_{max}$ = %.3f" % Rmax)
        ax2.set_ylabel(u"R factor")
        ax2.set_xlabel(u"Model number")
        ax2.set_xticks(xticks)
        ax2.legend(loc=8)
#        fig.set_title("Selection of dammif models based on \u03C7$^2$")
        self.valid = (chi2 < chi2max) * (R < Rmax)
        self.mask2d = (1 - numpy.identity(self.dammif_jobs)) * numpy.outer(self.valid, self.valid)
#        print self.valid
        bbox_props = dict(fc="pink", ec="r", lw=1)
        for i in range(self.dammif_jobs):
            if not self.valid[i]:
                ax1.text(i + 0.95, chi2max / 2, "Discarded", ha="center", va="center", rotation=90, size=10, bbox=bbox_props)
                ax2.text(i + 0.95, Rmax / 2, "Discarded", ha="center", va="center", rotation=90, size=10, bbox=bbox_props)
        if filename:
            filename = os.path.join(self.getWorkingDirectory(), filename)
            self.log("Wrote %s" % filename)
            fig.savefig(filename)
        if close:
            fig.clf()
            plt.close(fig)
        else:
            return fig

    def makeNSDarray(self, filename=None, close=True):
        self.arrayNSD = numpy.zeros(self.mask2d.shape, numpy.float32)
        fig = plt.figure(figsize=(15, 10))
        ax1 = fig.add_subplot(1, 2, 1)
        # for now just an empty figure but a placeholder
        ax1.imshow(self.arrayNSD, interpolation="nearest", origin="upper")

        xticks = 1 + numpy.arange(self.dammif_jobs)
        lnsd = []
        for key, plugin in self.supcomb_plugins.items():
            i0, i1 = key
            nsd = plugin.dataOutput.NSD.value
            self.arrayNSD[i0, i1] = nsd
            self.arrayNSD[i1, i0] = nsd
            lnsd.append(nsd)
            ax1.text(i0, i1, "%.2f" % nsd, ha="center", va="center", size=12 * 8 // self.dammif_jobs)
            ax1.text(i1, i0, "%.2f" % nsd, ha="center", va="center", size=12 * 8 // self.dammif_jobs)
        if len(lnsd) == 0: 
            return
        lnsd = numpy.array(lnsd)
#        print lnsd
#        print lnsd.mean() , lnsd.std(), lnsd.mean() + 2 * lnsd.std()
        nsd_max = lnsd.mean() + lnsd.std()
        data = self.arrayNSD.sum(axis=-1) / self.mask2d.sum(axis=-1)
        best_val = data[data > 0].min()
#        print data
#        print best_val
#        print numpy.where(data == best_val)
        self.ref = int(numpy.where(data == best_val)[0][-1])
#        print self.ref
        ax1.imshow(self.arrayNSD, interpolation="nearest", origin="upper")
        ax1.set_title(u"NSD correlation table")
        ax1.set_xticks(range(self.dammif_jobs))
        ax1.set_xticklabels([str(i) for i in range(1, 1 + self.dammif_jobs)])
        ax1.set_xlim(-0.5, self.dammif_jobs - 0.5)
        ax1.set_ylim(-0.5, self.dammif_jobs - 0.5)
        ax1.set_yticks(range(self.dammif_jobs))
        ax1.set_yticklabels([str(i) for i in range(1, 1 + self.dammif_jobs)])
        ax1.set_xlabel(u"Model number")
        ax1.set_ylabel(u"Model number")
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.bar(xticks - 0.5, data)
        ax2.plot([0.5, self.dammif_jobs + 0.5], [nsd_max, nsd_max], "-r", label=u"NSD$_{max}$ = %.2f" % nsd_max)
        ax2.set_title(u"NSD between any model and all others")
        ax2.set_ylabel("Normalized Spatial Discrepancy")
        ax2.set_xlabel(u"Model number")
        ax2.set_xticks(xticks)
        bbox_props = dict(fc="cyan", ec="b", lw=1)
        ax2.text(self.ref + 0.95, data[self.ref] / 2, "Reference", ha="center", va="center", rotation=90, size=10, bbox=bbox_props)
        ax2.legend(loc=8)
        self.valid *= (data < nsd_max)
        bbox_props = dict(fc="pink", ec="r", lw=1)
        for i in range(self.dammif_jobs):
            if not self.valid[i]:
                ax2.text(i + 0.95, data[self.ref] / 2, "Discarded", ha="center", va="center", rotation=90, size=10, bbox=bbox_props)
#        print self.valid
#        print self.ref
        if filename:
            filename = os.path.join(self.getWorkingDirectory(), filename)
            self.log("Wrote %s" % filename)
            fig.savefig(filename)
        if close:
            fig.clf()
            plt.close(fig)
        else:
            return fig
