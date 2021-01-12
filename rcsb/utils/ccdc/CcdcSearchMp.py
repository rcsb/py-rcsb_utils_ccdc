##
#
# File:    CcdcSearhMp.py
# Author:  J. Westbrook
# Date:    15-Dec-2020
# Version: 0.001
#
# Updated:
#
##
"""
MP execution wrapper for substructure and similarity search against the CCDC local Python API -

NOTE:  The latest ccdc API appears not compatible with rcsb.utils.multiproc and Python multiprocessing modules.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

# pylint: disable=redefined-outer-name

import logging
import time
import os
import os.path

from rcsb.utils.ccdc.CcdcSearch import CcdcSearch
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)


class CcdcSearchWorker(object):
    def __init__(self, verbose=True):
        self.__verbose = verbose

    def __checkStop(self, path):
        try:
            if os.access(path, os.F_OK):
                return True
        except Exception:
            pass
        return False

    def search(self, dataList, procName, optionsD, workingDir):
        """"""
        resultPath = optionsD["resultPath"]
        searchType = optionsD["searchType"]
        _ = workingDir
        ky = None
        startTime = time.time()
        logger.info("starting %s at %s", procName, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        resultList = []
        try:
            stopPath = os.path.join(resultPath, "STOP")
            logger.info("%s search list length %d", procName, len(dataList))

            for ii, queryTargetPath in enumerate(dataList):
                _, fn = os.path.split(queryTargetPath)
                queryTargetId, _ = os.path.splitext(fn)
                parentId = queryTargetId.split("|")[0]
                logger.info("%s starting %s (%s) (%d of %d)", procName, parentId, queryTargetId, ii + 1, len(dataList))
                if self.__checkStop(stopPath):
                    logger.info("%s stopping at %d of %d", procName, ii + 1, len(dataList))
                    break
                logger.info("%s search list %d of %d", procName, ii + 1, len(dataList))
                try:
                    vS = CcdcSearch(verbose=self.__verbose)
                    numHits = vS.search(queryTargetId, queryTargetPath, resultPath, searchType=searchType)
                    if numHits:
                        resultList.append(queryTargetPath)
                except Exception as e:
                    logger.exception("Failing for %r with %s", ky, str(e))

        except Exception as e:
            logger.exception("Failing for %r with %s", ky, str(e))

        endTime = time.time()
        logger.info(" %s completed at %s (%.2f seconds)", procName, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return resultList, [], []


class CcdcSearchMp(object):
    def __init__(self, verbose=True):
        """MP execution wrapper for variant search against the CCDC local Python API -"""
        self.__verbose = verbose
        #

    def runSearch(self, molFilePathList, resultPath, searchType="similarity", numProc=4, chunkSize=10):
        """"""
        logger.info("Starting with molfile path list length %d", len(molFilePathList))
        try:
            pU = CcdcSearchWorker(verbose=self.__verbose)
            mpu = MultiProcUtil(verbose=True)
            # mpu.setWorkingDir(resultPath)
            mpu.setOptions(optionsD={"resultPath": resultPath, "searchType": searchType})
            #
            mpu.set(workerObj=pU, workerMethod="search")

            ok, failList, resultList, _ = mpu.runMulti(dataList=molFilePathList, numProc=numProc, numResults=1, chunkSize=chunkSize)
            logger.info("run ended status %r success count %d failures %r", ok, len(resultList[0]), len(failList))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
