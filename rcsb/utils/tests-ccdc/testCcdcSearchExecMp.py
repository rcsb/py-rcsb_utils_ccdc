##
#
# File:    testCcdcSearchMp.py
# Author:  J. Westbrook
# Date:    15-Jan-2021
# Version: 0.001
#
# Updated:
#
##
"""
Test cases for chemical component search (mp) against the CCDC local Python API -
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

from rcsb.utils.ccdc.CcdcSearchExecMp import CcdcSearchExecMp

from rcsb.utils.ccdc import __version__


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class CcdcSearchMpTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__molFilePath = os.path.join(self.__dataPath, "molfiles")
        self.__pythonRootPath = os.path.join(os.environ["CSD_PYTHON_ROOT_PATH"])
        self.__csdHome = os.environ["CSDHOME"]
        #
        self.__simResultPath = os.path.join(self.__workPath, "test_chem_comp_ccdc_sim")
        self.__ssResultPath = os.path.join(self.__workPath, "test_chem_comp_ccdc_ss_exec")
        #
        self.__startTime = time.time()
        logger.info("Starting %s (%s) at %s", self.id(), __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSubStructureSearchExecMp(self):
        """Test case:  CCDC substructure search"""
        try:
            pL = glob.glob(os.path.join(self.__molFilePath, "*.mol2"), recursive=True)
            logger.info("search list length %d", len(pL))
            #
            csmp = CcdcSearchExecMp(pythonRootPath=self.__pythonRootPath, csdHome=self.__csdHome)
            csmp.runSearch(pL, self.__ssResultPath, searchType="substructure", numProc=2, chunkSize=2)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CcdcSearchMpTests("testSubStructureSearchExecMp"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
