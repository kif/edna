#!/usr/bin/python
# Written by Jerome Kieffer <jerome.kieffer@esrf.fr> 18/07/2017
#Compress the logs fro edna in /nobackup to free inodes

#import sys
import os
import glob
import shutil
import tarfile
import logging
import numpy

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("compress_log")

for root in glob.glob("/nobackup/*"):
    if not os.path.isdir(root):
        continue
    logdirs = glob.glob(root+"/edna-????????T??????")
    logdirs.sort()
    contains_data =[ bool(os.listdir(i)) for i in logdirs]
    cumulative = numpy.cumsum(contains_data)
    last = cumulative[-1]
    if last<2:
        logger.info("Nothing to compress, only found %s", logdirs)
        continue
    actual_last = numpy.where(cumulative == (last - 1))[0][0]
    for logdir in logdirs[:actual_last]:
        dest = logdir + ".tar.bz2"
        logger.debug("Compressing %s", logdir)
        with tarfile.open(dest, mode="w:bz2") as tar:
            tar.add(logdir, arcname=logdir[len(root):])
        logger.debug("Remove directory %s", logdir)
        shutil.rmtree(logdir)

