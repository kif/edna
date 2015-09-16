#!/usr/bin/python
# coding: utf-8
#
#    Project: Edna Saxs
#             http://www.edna-site.org
#
#
#    Copyright (C) 2012-2015 ESRF
#
#    Principal author: Jérôme Kieffer
#
#
#    The MIT License (MIT)
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from __future__ import with_statement

__authors__ = ["Jérôme Kieffer"]
__license__ = "MIT"
__copyright__ = "ESRF"
__date__ = "15/09/2015"
__status__ = "Development"
__version__ = "0.2"

import os, sys, time, logging
from StringIO import  StringIO
from optparse import OptionParser
import numpy
from scipy import stats
import matplotlib
exe = sys.argv[0].lower()
if "autorg" in exe  or "testall" in exe or "plotting" in exe:
    matplotlib.use('Qt4Agg')
else:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.ndimage
from scipy.cluster.vq import kmeans, vq
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger("saxs")
timelog = logging.getLogger("timeit")
import collections
import functools
from threading import Semaphore


def timeit(func):
    def wrapper(*arg, **kw):
        '''This is the docstring of timeit:
        a decorator that logs the execution time'''
        t1 = time.time()
        res = func(*arg, **kw)
        timelog.warning("%s took %.3fs" % (func.func_name, time.time() - t1))
        return res
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


class memoized(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''
    def __init__(self, func):
        self.func = func
        self.cache = {}
        self.args = []
        self.sem = Semaphore()
    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            value = self.cache[args]
        else:
            value = self.func(*args)
        with self.sem:
            self.cache[args] = value
            if args in self.args:
                self.args.remove(args)
            self.args.append(args)
            if len(args) > 100: # Keep only the
                rm = self.args.pop(0)
                self.cache.pop(rm)
            return value
    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__
    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)


@memoized
def load_saxs(filename):
    """
    return q, I, stderr as a dict
    """
    data = None
    res = {}
    for i in range(10):  # skip up to 10 comments lines
        try:
            data = numpy.loadtxt(filename, skiprows=i)
        except Exception:
            pass
        else:
            break
    if data == None:
        raise RuntimeError("Unable to read input file")
    if data.ndim == 2 and data.shape[1] == 2:
        res["q"] = data[:, 0]
        res["I"] = data[:, 1]
    elif  data.ndim == 2 and data.shape[1] == 3:
        res["q"] = data[:, 0]
        res["I"] = data[:, 1]
        res["std"] = data[:, 2]
    else:
        raise RuntimeError("Unable to find columns in data file")

    return res

@memoized
def load_gnom(filename):
    """

    @param filename: path of the Gnom output File
    @return: dict with many parameters: gnomRg, gnomRg_err, gnomI0, gnomI0_err, q_fit, I_fit, r, P(r), P(r)_err

    """
    if isinstance(filename, dict):
        return filename
    pr = StringIO("")
    reg = StringIO("")
    do_pr = False
    do_reg = False
    out = {}
    with open(filename, "r") as logLines:
        for idx, line in enumerate(logLines):
            words = line.split()
            if "Total  estimate" in line:
                out["fit_quality"] = float(words[3])
            if "Reciprocal space:" in line:
                do_pr = False
            if "Distance distribution" in line :
                do_reg = False
            if words:
                if do_reg:
                    reg.write("%s %s\n" % (words[0], words[-1]))
                if do_pr:
                    pr.write(line)
            if "I REG" in line:
                do_reg = True
            if "P(R)" in line:
                do_pr = True
    out["gnomRg"] = float(words[4])
    out["gnomRg_err"] = float(words[6])
    out["gnomI0"] = float(words[9])
    out["gnomI0_err"] = float(words[11])
    reg.seek(0)
    pr.seek(0)
    out["q_fit"], out["I_fit"] = numpy.loadtxt(reg, unpack=True, dtype="float32")
    out["r"], out["P(r)"], out["P(r)_err"] = numpy.loadtxt(pr, unpack=True, dtype="float32")
    return out
loadGnomFile = load_gnom

def density_plot(gnomfile, filename=None, format="png", unit="nm", ax=None):
    """
    Generate a density plot P = f(r) 

    @param gnomfile: name of the GNOM output file
    @param  filename: name of the file where the cuve should be saved
    @param format: image format
    @param ax: subplotib where to plot in
    @return: the matplotlib figure
    """
    out = load_gnom(gnomfile)
    if ax:
        fig = ax.figure
    else:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1)
    ax.errorbar(out["r"], out["P(r)"], out["P(r)_err"], label="Density", capsize=0, color="blue", ecolor="lightblue")
    ax.set_ylabel('$\\rho (r)$')
    ax.set_xlabel('$r$ (%s)' % unit)
    ax.set_title("Pair distribution function")
#    ax.set_yscale("log")
    ax.legend()
#    ax.legend(loc=3)
    if filename:
        if format:
            fig.savefig(filename, format=format)
        else:
            fig.savefig(filename)
    return fig
densityPlot = density_plot


def scatter_plot(curve_file, first_point=None, last_point=None, filename=None, format="png", unit="nm", gnomfile=None, ax=None):
    """
    Generate a scattering plot I = f(q) in semi log.

    @param curve_file: name of the saxs curve file
    @param: first_point,last point: integers, by default 0 and -1
    @param  filename: name of the file where the cuve should be saved
    @param format: image format
    @param ax: subplot where to plot in
    @return: the matplotlib figure
    """
    data = load_saxs(curve_file)
    q = data.get("q")
    I = data.get("I")
    std = data.get("std")
    if (first_point is None) and (last_point is None):
        for line in open(curve_file):
            if "# AutoRg: Points" in line:
                d = [int(i) for i in line.split() if i.isdigit()]
                if len(d) >= 2:
                    first_point = d[0]
                    last_point = d[1] + 1
    if first_point is None:
        first_point = 0
    if last_point is None:
        last_point = -1

    if ax:
        fig = ax.figure
    else:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1)
    if std is not None:
        ax.errorbar(q, I, std, label="Experimental curve", capsize=0, color="blue", ecolor="lightblue")
    else:
        ax.plot(q, I, label="Experimental curve", color = "blue")
    if first_point is not None:
        if std is not None:
            ax.errorbar(q[first_point:], I[first_point:], std[first_point:], label="Experimental curve (cropped)", capsize=0, color="green", ecolor="lightgreen")
        else:
            ax.plot(q[first_point:], I[first_point:], label="Experimental curve (cropped)", color = "green")
    if gnomfile:
        gnom = load_gnom(gnomfile)
        ax.plot(gnom["q_fit"], gnom["I_fit"], label="GNOM fitted curve", color="red")
    ax.set_ylabel('$I(q)$')
    ax.set_xlabel('$q$ (%s$^{-1}$)' % unit)
    ax.set_title("Scattering curve")
    ax.set_yscale("log")
    ax.legend()
#    ax.legend(loc=3)
    if filename:
        if format:
            fig.savefig(filename, format=format)
        else:
            fig.savefig(filename)
    return fig
scatterPlot = scatter_plot

def guinier_plot(curve_file, first_point=None, last_point=None, filename=None, format="png", unit="nm", ax=None):
    """
    Generate a Guinier plot Ln(I) vs q**2
    @param curve_file: name of the saxs curve file
    @param: first_point,last point: integers, by default 0 and -1
    @param  filename: name of the file where the cuve should be saved
    @param format: image format
    @param: ax: subplot where to plot in
    @return: the matplotlib figure
    """
    data = load_saxs(curve_file)
    q = data.get("q")
    I = data.get("I")
    std = data.get("std")
    if (first_point is None) and (last_point is None):
        for line in open(curve_file):
            if "# AutoRg: Points" in line:
                d = [int(i) for i in line.split() if i.isdigit()]
                if len(d) >= 2:
                    first_point = d[0]
                    last_point = d[1] + 1
    if first_point is None:
        first_point = 0
    if last_point is None:
        last_point = -1

    q2 = q * q
    logI = numpy.log(I)

    slope, intercept, r_value, p_value, std_err = stats.linregress(q2[first_point:last_point], logI[first_point:last_point])
    Rg = numpy.sqrt(-3 * slope)
    I0 = numpy.exp(intercept)
    end = min(q.size, (-1.5 / slope > q).sum())
    q = q[:end]
    I = I[:end]
    q2 = q2[:end]
    logI = logI[:end]

    if ax:
        fig = ax.figure
    else:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1)
    ax.plot(q2, logI, label="Experimental curve")
    ax.plot(q2[first_point:last_point], logI[first_point:last_point], marker='D', markersize=5, label="Guinier region")
    ax.annotate("qRg$_{min}$=%.1f" % (Rg * q[first_point]), (q2[first_point], logI[first_point]), xytext=None, xycoords='data',
         textcoords='data')
    ax.annotate("qRg$_{max}$=%.1f" % (Rg * q[last_point]), (q2[last_point], logI[last_point]), xytext=None, xycoords='data',
         textcoords='data')

    ax.plot(q2, intercept + slope * q2, label="ln[$I(q)$] = %.2f %.2f * $q^2$" % (intercept, slope))
    ax.set_ylabel('ln[$I(q)$]')
    ax.set_xlabel('$q^2$ (%s$^{-2}$)' % unit)
    ax.set_title("Guinier plot: $Rg=$%.1f %s $I0=$%.1f" % (Rg, unit, I0))
    ax.legend(loc=3)
    if filename:
        if format:
            fig.savefig(filename, format=format)
        else:
            fig.savefig(filename)
    return fig
guinierPlot = guinier_plot


def kartky_plot(curve_file, filename=None, format="png", unit="nm", ax = None):
    """
    Generate a Kratky: q2I(q) vs q
    @param curve_file: name of the saxs curve file
    @param: first_point,last point: integers, by default 0 and -1
    @param  filename: name of the file where the cuve should be saved
    @param format: image format
    @param ax: subplot to display in 
    @return: the matplotlib figure
    """
    data = load_saxs(curve_file)
    q = data.get("q")
    I = data.get("I")
    std = data.get("std")
    q2I = q * q * I
    if ax:
        has_subplot = True
    if has_subplot:
        fig = ax.figure
    else:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1)
    ax.plot(q, q2I, label="Experimental curve")
    ax.set_ylabel('$q^2I (%s^2)$' % unit)
    ax.set_xlabel('$q$ (%s$^{-1}$)' % unit)
    ax.set_title("Kratky plot")
    ax.legend(loc=0)
    if filename:
        if format:
            fig.savefig(filename, format=format)
        else:
            fig.savefig(filename)
    return fig
kartkyPlot = kartky_plot

class AutoRg(object):
    """
    a class to calculate Automatically the Radius of Giration based on Guinier approximation.
    """
    def __init__(self, q=None, I=None, std=None, datfile=None, mininterval=10, qminRg=1.0, qmaxRg=1.3):
        self.q = q
        self.I = I
        self.std = std
        if (q is None) or (I is None) and datfile:
            dico = load_saxs(datfile)
            self.q = dico.get("q") 
            self.I = dico.get("I")
            self.std = dico.get("std")

        self.mininterval = mininterval
        self.qminRg = qminRg
        self.qmaxRg = qmaxRg
        self.results = {}
        self.start_search = 0
        self.stop_search = len(self.q)
        self.len_search = len(self.q)
        self.result = {}
        self.n = None
        self.start = None
        self.stop = None
        self.Sx = None
        self.Sy = None
        self.Sxx = None
        self.Sxy = None
        self.Syy = None
        self.Sw = None
        self.slope = None
        self.intercept = None
        self.I0 = None
        self.correlationR = None
        self.sterrest = None
        self.best = None
        self.big_dim = None
        self.dI0 = None
        self.dRg = None

    @timeit
    def select_range(self):
        """
        First step: limit the range of search:

        * remove all point before maximum I
        * keep only up to Imax/10
        * if some points have I<0 segment region and keep the largest sub-region

        """
        self.start_search = self.I.argmax()
        Imax = self.I[self.start_search]
        keep = (self.I > (Imax / 10.0))
        keep[:self.start_search] = False
        if self.I[keep].min() <= 0:
            logger.debug("Negatives values in search range: refining")
            keep[self.I <= 0] = False
            label = scipy.ndimage.label(keep)
            lab_max = label.max()
            res = [ 0 ]
            for idx in range(1, lab_max + 1):
                res.append((label == idx).sum() / idx)
            largest_region = numpy.array(res).argmax()
            keep = (label == largest_region)
        self.start_search = keep.argmax()
        self.len_search = keep.sum()
        self.stop_search = self.start_search + self.len_search
        logger.debug("Searching range: %i -> %i (%i points)" % (self.start_search, self.stop_search, self.len_search))

    @timeit
    def allocate(self):
        """
        Allocate 3 big buffers:
        x = q*q
        y = log(I)
        w = 1/dy = I/std

        calculate the sums: xw, wy, w, wxy, wxx and wyy
        """
        self.big_dim = (self.len_search - self.mininterval + 1) * (self.len_search - self.mininterval) / 2  # + len_search * mininterval
        array_size = self.big_dim * self.len_search * 8 / 1e6
        if array_size > 1000:
            print("Allocating large array: %.3f MB !!!! expect to fail" % array_size)
#        try:
        x = numpy.zeros((self.big_dim, self.len_search), dtype="float64")
        y = numpy.zeros((self.big_dim, self.len_search), dtype="float64")
        w = numpy.zeros((self.big_dim, self.len_search), dtype="float64")  # (1/dy = 1/(d(logI)=I/std)
        self.n = numpy.zeros(self.big_dim, dtype="int16")
        self.start = numpy.zeros(self.big_dim, dtype="int16")
        self.stop = numpy.zeros(self.big_dim, dtype="int16")
#        except MemoryError as error:q2
        q2 = self.q * self.q
        logI = numpy.log(self.I)

        if self.std is not None:
            I_over_std = self.I / self.std
        else:
            I_over_std = numpy.ones_like(self.I)

        idx = 0
        for sta in range(self.start_search, self.start_search + self.len_search - self.mininterval):
            for sto in range(sta + self.mininterval, self.start_search + self.len_search):
                x[idx, sta - self.start_search:sto - self.start_search] = q2[sta :sto]
                y[idx, sta - self.start_search:sto - self.start_search] = logI[sta :sto]
                w[idx, sta - self.start_search:sto - self.start_search] = I_over_std[sta :sto]
                self.n[idx] = sto - sta
                self.start[idx] = sta
                self.stop[idx] = sto
                idx += 1
        del q2, logI, I_over_std
        self.Sx = (w * x).sum(axis= -1)
        self.Sy = (w * y).sum(axis= -1)
        self.Sxx = (w * x * x).sum(axis= -1)
        self.Sxy = (w * y * x).sum(axis= -1)
        self.Syy = (w * y * y).sum(axis= -1)
        self.Sw = w.sum(axis= -1)
        del x, y, w

    @timeit
    def refine(self):
        """
        Keep only ranges with valid qminRg and qmaxRg.

        Calculate Rg, I0 and the linear regression quality fit.
        """

        self.slope = (self.Sw * self.Sxy - self.Sx * self.Sy) / (self.Sw * self.Sxx - self.Sx * self.Sx)
        self.Rg = numpy.sqrt(-self.slope * 3)
        valid = numpy.logical_and((self.Rg * self.q[self.start] <= self.qminRg) , (self.Rg * self.q[self.stop - 1] <= self.qmaxRg))
        valid[self.slope > 0] = False
        nvalid = valid.sum()
        if nvalid > 0:
            for ds in ("start", "stop", "n", "slope", "Rg", "Sx", "Sy", "Sw", "Sxx", "Sxy", "Syy"):
                setattr(self, ds, getattr(self, ds)[valid])
            self.intercept = (self.Sy - self.Sx * self.slope) / self.Sw
            self.I0 = numpy.exp(self.intercept)
            df = self.n - 2
            r_num = ssxym = (self.Sw * self.Sxy) - (self.Sx * self.Sy)
            ssxm = self.Sw * self.Sxx - self.Sx * self.Sx
            ssym = self.Sw * self.Syy - self.Sy * self.Sy
            r_den = numpy.sqrt(ssxm * ssym)
            self.correlationR = r_num / r_den
            self.correlationR[r_den == 0] = 0.0
            self.correlationR[self.correlationR > 1.0] = 1.0  # Numerical errors
            self.correlationR[self.correlationR < -1.0] = -1.0  # Numerical errors
            self.sterrest = numpy.sqrt((1.0 - self.correlationR * self.correlationR) * ssym / ssxm / df)
            var_slope = self.Sw / ssxm
            var_intercept = self.Sxx / ssxm
            self.dI0 = numpy.sqrt(var_intercept) * self.I0
            self.dRg = 1.5 * numpy.sqrt(var_slope) / self.Rg
#            valid2 = (self.sterrest < (10 * self.sterrest.min()))
#            if valid2.sum() > 10:
#                logger.info("Cut off std_max = 4*std_min")
#                for ds in ("start", "stop", "n", "slope", "Rg", "Sx", "Sy", "Sw", "Sxx", "Sxy", "Syy", "sterrest", "dI0", "dRg", "intercept", "I0", "correlationR"):
#                    setattr(self, ds, getattr(self, ds)[valid2])


    @timeit
    def cluster(self):
#        ratio = self.Rg.std() / self.Rg.mean()
#        logger.warning("Rg.std: %s ratio %s" % (self.Rg.std(), ratio))
        features = numpy.hstack((self.Rg.reshape(-1, 1), self.I0.reshape(-1, 1)))
        centroids, variance = kmeans(features, 2)
        code, distance = vq(features, centroids)
        if logger.level <= logging.INFO:
            fig = matplotlib.pyplot.figure()
            ax = fig.add_subplot(111)
            code0 = numpy.where(code == 0)
            code1 = numpy.where(code == 1)
            ax.plot(self.Rg[code0], self.I0[code0], 'b*')
            ax.plot(self.Rg[code1], self.I0[code1], 'r*')
            ax.plot([p[0] for p in centroids], [p[1] for p in centroids], 'go')
            fig.show()
        code_min_start = numpy.array(centroids[:, 0]).argmin()  # code with largest Rg
        valid = numpy.where(code == code_min_start)
        for ds in ("start", "stop", "n", "slope", "Rg", "Sx", "Sy", "Sw", "Sxx", "Sxy", "Syy", "sterrest", "dI0", "dRg", "intercept", "I0", "correlationR"):
            setattr(self, ds, getattr(self, ds)[valid])




    @timeit
    def finish(self):
        if self.sterrest is not None:
            best = self.best = self.sterrest.argmin()
            sta = self.start[best]
            sto = self.stop[best]
            res = {"start":sta, "end":sto,
                   "Rg":self.Rg[best], "logI0":self.I0[best],
                   "R":self.correlationR[best], "stderr":self.sterrest[best],
                   "len":sto - sta,
                   "I0":self.I0[best],
                   "qminRg":self.Rg[best] * self.q[sta],
                   "qmaxRg":self.Rg[best] * self.q[sto - 1],
                   }
#            if (sto - sta) > self.mininterval:
#                shift = numpy.where(numpy.logical_and(self.start >= (sta), self.stop <= (sto)))[0]
#            else:
#                shift = numpy.where(numpy.logical_and(self.start >= (sta - 1), self.stop <= (sto + 1)))[0]
            res["deltaRg"] = self.dRg[best]
            res["deltaI0"] = self.dI0[best]
            res["start_search"] = self.start_search
            res["stop_search"] = self.start_search + self.len_search
            res["intervals"] = self.big_dim

            parab = lambda p, x, y: p[0] * x * x + p[1] * x + p[2] - y
            out = scipy.optimize.leastsq(parab, [0, self.slope[best], self.intercept[best]], (self.q[sta:sto] * self.q[sta:sto], numpy.log(self.I[sta:sto])))
            if out[0][0] > 0:
                res["Aggregated"] = True
            else:
                res["Aggregated"] = False
            self.result = res
            return res


def autoRg(q=None, I=None, std=None, datfile=None, mininterval=10, qminRg=1.0, qmaxRg=1.3):
    ag = AutoRg(q, I, std, datfile, mininterval, qminRg, qmaxRg)
    ag.select_range()
    ag.allocate()
    ag.refine()
    ag.cluster()
#    if logger.level <= logging.INFO:
#            raw_input("Enter to quit")

#    fig.show()

    return ag.finish()

def gnom(fname, Rg, I0):
    """
    Exploratory work for Gnom replacement
    """
    data = load_saxs(curve_file)
    q = data.get("q")
    I = data.get("I")
    err = data.get("std")
    qstep = q[1] - q[0]
    q.max()
    q.min()
    rmax = 1. / q.min()
    rmin = 1 / q.max()
    rstep = rmin
    lostep = int(numpy.ceil((rmin - 1e-08) / rstep))
    histep = int(numpy.floor((rmax + 1e-08) / rstep)) + 1
    qmaxrstep = numpy.pi / rstep
    nbase = max((len(q), histep, qmaxrstep / qstep))
    nlog2 = int(numpy.ceil(numpy.log2(nbase)))
    nout = 2 ** nlog2
    qmaxdb = 2 * nout * qstep
    q_full = numpy.linspace(0, 6, nout)
    guinier = I0 * numpy.exp(-Rg * Rg * q_full * q_full / 3)
    I_interp = numpy.interp(q_full, q, I, 0, 0)
    (numpy.log(I_interp) - numpy.log(guinier))
    I_interp = numpy.interp(q_full, q, I, 0, 0)
    delta = (numpy.log(I_interp) - numpy.log(guinier))
    delta.min()
    delta.max()
    delta[delta == -numpy.inf] = delta.max()
    switch = numpy.argmin(delta)
    I_merge = guinier
    I_merge[switch:] = I_interp[switch:]
    crho = numpy.fft.ifft(I_merge) ** 2
    plot(numpy.imag(crho))

#        from numpy.fft import ifft
#        lostep = int(numpy.ceil((self.rmin - 1e-08) / self.rstep))
#        histep = int(numpy.floor((self.rmax + 1e-08) / self.rstep)) + 1
#        self.xout = numpy.arange(lostep, histep) * self.rstep
#        self.qstep = self.xin[1] - self.xin[0]
#        self.qmaxrstep = numpy.pi / self.rstep
#        nin = len(self.xin)
#        nbase = max([nin, histep, self.qmaxrstep / self.qstep])
#        nlog2 = int(numpy.ceil(numpy.log2(nbase)))
#        nout = 2 ** nlog2
#        qmaxdb = 2 * nout * self.qstep
#        yindb = numpy.concatenate((self.yin, numpy.zeros(2 * nout - nin)))
#        cyoutdb = ifft(yindb) * 2 / numpy.pi * qmaxdb
#        youtdb = numpy.imag(cyoutdb)
#        xstepfine = 2 * numpy.pi / qmaxdb
#        xoutfine = numpy.arange(nout) * xstepfine
#        youtfine = youtdb[:nout]
#        self.yout = numpy.interp(self.xout, xoutfine, youtfine)


def plotting(curve_file, gnom_file=None, filename=None, format="png", unit="nm"):
    """
    Generate a plot with scattering, Guinier plot, Kratky, anf Gnom plot.
    
    @param curve_file: name of the saxs curve file
    @param filename: name of the file where the cuve should be saved
    @param format: image format
    @return: the matplotlib figure
    """
    if gnom_file is None and os.path.exists(curve_file[:-3]+"out"):
        gnom_file = curve_file[:-3]+"out"
    fig = plt.figure(figsize=(12,10))
    sp_scat = fig.add_subplot(2,2,1)
    scatter_plot(curve_file, gnomfile=gnom_file, ax=sp_scat)
    sp_krat = fig.add_subplot(2,2,2)
    kartky_plot(curve_file, ax=sp_krat)
    sp_guin = fig.add_subplot(2,2,3)
    guinier_plot(curve_file, ax=sp_guin)
    sp_gnom = fig.add_subplot(2,2,4)
    density_plot(gnomfile=gnom_file, ax=sp_gnom)
    if filename:
        if format:
            fig.savefig(filename, format=format)
        else:
            fig.savefig(filename)
    elif matplotlib.get_backend()!="Agg":
        fig.show()
    return fig


def main(exe):
    if "autorg" in exe:
            usage = """autorg.py [OPTIONS] <DATAFILE(S)>

Estimation of radius of gyration from SAS data by Guinier approximation.

Output: estimated Rg, standard deviation of Rg, estimated I(0), standard
deviation of I(0), first point of chosen Guinier interval, last point of chosen
Guinier interval, estimated data quality, aggregation check.

Known options:
  -h, --help                print this help, then exit
  -v, --version             print version information, then exit
  -o, --output <FILE>       relative or absolute file path to save result
  -f, --format <FORMAT>     output format, one of: csv, ssv, table
      --mininterval <VALUE> minimum interval length in points (default: 10)
      --smaxrg <VALUE>      maximum Smax*Rg value (default: 1.3)
      --sminrg <VALUE>      maximum Smin*Rg value (default: 1.0)

Report bugs to <jerome.kieffer@esrf.fr>.
"""
            version = "%prog " + __version__
            parser = OptionParser(usage=usage, version=version)
            parser.add_option("-o", "--output", dest="output",
                              type='string', default=None,
                              help="relative or absolute file path to save result")
            parser.add_option("-f", "--format", dest="format",
                              type="string", default=None,
                              help="output format, one of: csv, ssv, table")
            parser.add_option("--mininterval", dest="mininterval",
                              type="int", default=10,
                              help="minimum interval length in points (default: 10)")
            parser.add_option("--smaxrg", dest="smaxrg",
                              type="float", default=1.3,
                              help="maximum Smax*Rg value (default: 1.3)")
            parser.add_option("--sminrg", dest="sminrg",
                              type="float", default=1.0,
                              help="maximum Smin*Rg value (default: 1.0)")
            (options, args) = parser.parse_args()
            if len(args) < 1:
                parser.error("incorrect number of arguments: use -h to get help")
            for afile in args:
                if os.path.isfile(afile):
                    r = autoRg(datfile=afile, mininterval=options.mininterval, qminRg=options.sminrg, qmaxRg=options.smaxrg)
                    if r:
                        print """Rg   =  %5.2f  +/- %.2f (%i%%)
I(0) =  %5.1f +/- %.2f
Points   %i to %i (%i total)""" % (r["Rg"], r["deltaRg"], 100.0 * r["deltaRg"] / r["Rg"], r["I0"], r["deltaI0"], r["start"] + 1, r["end"] , r["len"])
                        if r.get("Aggregated", None):
                            print "Aggregated."
                        print """(Searched from point %i to %i, %i intervals analysed)""" % (r["start_search"] + 1, r["stop_search"], r["intervals"])
                    else:
                        print "No Rg found for '%s'." % afile


#                    fig = guinier_plot(curve_file=afile)

                else:
                    print("No such file %s" % afile)
#            plt.show()
    elif "testall" in exe:
        for afile in sys.argv[1:]:
            for line in open(afile):
                if "# AutoRg: Points" in line:
                    d = [int(i) for i in line.split() if i.isdigit()]
                    if len(d) >= 2:
                        first_point = d[0]
                        last_point = d[1]
                        break
                elif "# AutoRg: Rg" in line:
                    Rg = float(line.split()[4])
                elif "# AutoRg: I(0)" in line:
                    I0 = float(line.split()[4])
            t0 = time.time()
            r = autoRg(datfile=afile)
            t1 = time.time()
            if r:
                print("%s: Rg=%.2f(%.2f) I0=%.2f(%.2f) on %i(%i) -> %i(%i). took %.3fs" % (afile, r["Rg"], Rg, r["I0"], I0, r["start"] + 1, first_point, r["end"], last_point, t1 - t0))
            else:
                if first_point == 0 and 0 == last_point:
                    print("OK %s: No Rg, was %i %i.\t took %.3fs" % (afile, first_point, last_point, t1 - t0))
                else:
                    print("!! %s: No Rg, was %i %i.\t took %.3fs" % (afile, first_point, last_point, t1 - t0))
    elif "plotting" in exe:
        for afile in sys.argv[1:]:
            plotting(afile)
        raw_input("Enter to quit")

if __name__ == "__main__":
    main(exe)
    