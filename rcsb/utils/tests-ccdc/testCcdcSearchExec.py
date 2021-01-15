##
#
# File:    testCcdcSearchExec.py
# Author:  J. Westbrook
# Date:    16-Dec-2020
# Version: 0.001
#
# Updated:
#
##
"""
Test cases for execution of CCDC search cli -

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

from rcsb.utils.io.ExecUtils import ExecUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.ccdc import __version__


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class CcdcSearchExecTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__molFileDirPath = os.path.join(self.__dataPath, "molfiles")
        self.__ssResultPath = os.path.join(self.__workPath, "test_chem_comp_ccdc_ss_cli")
        self.__logPath = os.path.join(self.__workPath, "execStdout.log")
        # path to the CSD and python interpreter where the CCDC api is installed.
        self.__pythonBinPath = os.path.join(os.environ["CSD_PYTHON_ROOT_PATH"], "bin", "python")
        self.__csdHome = os.environ["CSDHOME"]
        #
        self.__queryListFilePath = os.path.join(self.__workPath, "query_list.txt")

        self.__startTime = time.time()
        logger.info("Starting %s (%s) at %s", self.id(), __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSearchExec(self):
        """Test case:  search cli"""
        try:
            mL = glob.glob(os.path.join(self.__molFileDirPath, "*"))
            logger.info("search list length %d", len(mL))
            mU = MarshalUtil()
            ok = mU.doExport(self.__queryListFilePath, mL, fmt="list")
            exU = ExecUtils()
            logger.info("Executing shell for %s", self.__queryListFilePath)
            cmdPath = os.path.join(TOPDIR, "rcsb", "utils", "ccdc", "CcdcSearchExec.py")

            logger.info("cmdPath %r", cmdPath)
            ok = exU.runShell(
                "%s %s --mol_list_path %s --result_path %s --search_type %s --csdhome %s"
                % (self.__pythonBinPath, cmdPath, self.__queryListFilePath, self.__ssResultPath, "substructure", self.__csdHome),
                outPath=self.__logPath,
                outAppend=False,
                timeOut=60,
                suppressStderr=False,
            )
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteSearchExecTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CcdcSearchExecTests("testSearchExec"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteSearchExecTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
