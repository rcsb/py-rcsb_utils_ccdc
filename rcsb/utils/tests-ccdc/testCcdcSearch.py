##
#
# File:    testCcdcSearch.py
# Author:  J. Westbrook
# Date:    13-Dec-2020
# Version: 0.001
#
# Updated:
#
##
"""
Test cases for chemical component search against the CCDC local Python API -

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import glob
import logging
import unittest
import time
import os
import os.path
import platform
import resource

from rcsb.utils.ccdc.CcdcSearch import CcdcSearch
from rcsb.utils.ccdc import __version__


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class CcdcSearchTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__debug = True
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        # Path to a set of test mol2 files ...
        self.__molFilePath = os.path.join(self.__dataPath, "molfiles")
        # Test output paths
        self.__simResultPath = os.path.join(self.__workPath, "ccdc_sim")
        self.__ssResultPath = os.path.join(self.__workPath, "ccdc_ss_mol")
        self.__smartsResultPath = os.path.join(self.__workPath, "ccdc_ss_smarts")
        #
        self.__smartsList = [("000", "COC(=O)O")]
        self.__startTime = time.time()
        logger.info("Starting %s (%s) at %s", self.id(), __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSmartsSearch(self):
        """Test case:  CCDC SMARTS substructure"""
        try:
            #
            for queryTargetId, smarts in self.__smartsList:
                logger.info("search for %r", queryTargetId)
                resultPath = self.__smartsResultPath
                vS = CcdcSearch(verbose=self.__verbose)
                vS.searchSmarts(queryTargetId, smarts, resultPath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSimilaritySearch(self):
        """Test case:  CCDC similarity search"""
        try:
            pL = glob.glob(os.path.join(self.__molFilePath, "*.mol2"))
            logger.info("search list length %d", len(pL))
            #
            for queryTargetPath in pL:
                _, fn = os.path.split(queryTargetPath)
                queryTargetId, _ = os.path.splitext(fn)
                logger.info("search for %r", queryTargetId)
                resultPath = self.__simResultPath
                vS = CcdcSearch(verbose=self.__verbose)
                vS.search(queryTargetId, queryTargetPath, resultPath, searchType="similarity")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSubStructureSearch(self):
        """Test case:  CCDC substructure search"""
        try:
            pL = glob.glob(os.path.join(self.__molFilePath, "*.mol2"))
            logger.info("search list length %d", len(pL))
            #
            for queryTargetPath in pL:
                _, fn = os.path.split(queryTargetPath)
                queryTargetId, _ = os.path.splitext(fn)
                logger.info("search for %r", queryTargetId)
                resultPath = self.__ssResultPath
                vS = CcdcSearch(verbose=self.__verbose)
                vS.search(queryTargetId, queryTargetPath, resultPath, searchType="substructure")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CcdcSearchTests("testSimilaritySearch"))
    suiteSelect.addTest(CcdcSearchTests("testSubStructureSearch"))
    suiteSelect.addTest(CcdcSearchTests("testSmartsSearch"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
