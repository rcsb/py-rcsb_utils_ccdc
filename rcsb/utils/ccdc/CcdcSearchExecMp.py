##
#
# File:    CcdcSearchExecMp.py
# Author:  J. Westbrook
# Date:    15-Jan-2021
# Version: 0.001
#
# Updated:
#
##
"""
Multiprocessing wrapper for substructure and similarity search against the CCDC local Python API -

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

from rcsb.utils.io.ExecUtils import ExecUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)


class CcdcSearchExecWorker(object):
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
        """Worker method to execute a shell to search CCDC for the input mol2 path list.

        Args:
            dataList (list): list of mol2 file paths to be searched
            procName (str): processName
            optionsD (dict): dictionary of options
            workingDir (str): path to working directory (not used)

        Returns:
            (successList, resultList, []): success and result lists of mol2 paths with CCDC matches
        """
        resultPath = optionsD["resultPath"]
        searchType = optionsD["searchType"]
        pythonRootPath = optionsD["pythonRootPath"]
        csdHome = optionsD["csdHome"]
        _ = workingDir
        resultList = []
        startTime = time.time()
        logger.info("starting %s at %s", procName, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        #
        try:
            stopPath = os.path.join(resultPath, "STOP")
            logger.info("%s search list length %d", procName, len(dataList))
            if self.__checkStop(stopPath):
                logger.info("%s stopping", procName)
                return resultList, resultList, []
            #
            queryListFilePath = os.path.join(resultPath, procName, "queryFileList.list")
            mU = MarshalUtil()
            ok = mU.doExport(queryListFilePath, dataList, fmt="list")
            if not ok:
                return resultList, resultList, []
            #
            exU = ExecUtils()
            logger.info("%s executing shell for %s", procName, queryListFilePath)
            cmdPath = os.path.join(pythonRootPath, "bin", "ccdc_search_cli")
            hitListPath = os.path.join(resultPath, procName, "hitList.list")
            logPath = os.path.join(resultPath, procName, "execlog.log")

            logger.info("cmdPath %r", cmdPath)
            ok = exU.runShell(
                "%s --mol_list_path %s --result_path %s --search_type %s --csdhome %s --hit_list_path %s"
                % (cmdPath, queryListFilePath, resultPath, searchType, csdHome, hitListPath),
                outPath=logPath,
                outAppend=False,
                timeOut=60,
                suppressStderr=False,
            )
            #
            if ok and mU.exists(hitListPath):
                resultList = mU.doImport(hitListPath, fmt="list")
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        endTime = time.time()
        logger.info("%s (result len %d) completed at %s (%.2f seconds)", procName, len(resultList), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return resultList, resultList, []


class CcdcSearchExecMp(object):
    def __init__(self, pythonRootPath, csdHome, verbose=True):
        """MP execution wrapper for search against the CCDC local Python API -"""
        self.__verbose = verbose
        self.__pythonRootPath = pythonRootPath
        self.__csdHome = csdHome
        #

    def runSearch(self, molFilePathList, resultPath, searchType="similarity", numProc=4, chunkSize=10):
        """Run CCDC search in multiprocess mode.

        Args:
            molFilePathList (list): input mol2/sdf path list to search
            resultPath (str): directory path to store results
            searchType (str, optional): search type (substructure|similarity). Defaults to "similarity".
            numProc (int, optional): number of processes to invoke. Defaults to 4.
            chunkSize (int, optional): work chunksize. Defaults to 10.
        """
        logger.info("Starting with molfile path list length %d", len(molFilePathList))
        try:
            pU = CcdcSearchExecWorker(verbose=self.__verbose)
            mpu = MultiProcUtil(verbose=True)
            mpu.setWorkingDir(resultPath)
            mpu.setOptions(optionsD={"resultPath": resultPath, "searchType": searchType, "pythonRootPath": self.__pythonRootPath, "csdHome": self.__csdHome})
            #
            mpu.set(workerObj=pU, workerMethod="search")

            ok, failList, resultList, _ = mpu.runMulti(dataList=molFilePathList, numProc=numProc, numResults=1, chunkSize=chunkSize)
            logger.info("Run ended with status %r success count %d failures %r", ok, len(resultList[0]), len(failList))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
