##
#
# File:    testCcdcGeomAnal.py
# Author:  J. Westbrook
# Date:    17-Dec-2020
# Version: 0.001
#
# Updated:
#
##
"""
Test cases for chemical component geometrical analysis using CCDC Python API-
"""

import glob
import json
import logging
import os
import platform
import resource
import sys
import time
import unittest

from rcsb.utils.ccdc.CcdcGeomAnal import CcdcGeomAnal
from rcsb.utils.ccdc import __version__

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class CcdcGeomAnalTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stderr
        self.__verbose = True
        self.__debug = True
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__molFilePath = os.path.join(self.__dataPath, "molfiles-xyz")
        self.__resultPath = os.path.join(self.__workPath, "ccdc_anal")
        self.__useAtomMap = False
        self.__startTime = time.time()
        logger.info("Starting %s (%s) at %s", self.id(), __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def __getTargetIdList(self, fn="../data/ccid_exercise_list.txt"):
        idList = []
        try:
            fh = open(fn, "r")
            for line in fh:
                idList.append(str(line[:-1]).strip())
            fh.close()
        except Exception as e:
            logger.exception("FAILING with %s", str(e))

        return idList

    def testGeomAnalSettings(self):
        try:
            cga = CcdcGeomAnal(verbose=self.__verbose, log=self.__lfh)
            logger.info("Starting/default settings\n%s\n", cga.settings())
            cga.globalSettings(rfactor_filter=0.05, generalisation=True, solvent_filter="exclude_solvent", organometallic_filter="organics_only", heaviest_atom="Ni")
            cga.featureSettings("bond", zscore_threshold=4.0, few_hits_threshold=10)
            cga.featureSettings("angle", zscore_threshold=4.0, few_hits_threshold=10)
            cga.featureSettings("torsion", local_density_threshold=10.0, local_density_tolerance=15.0, few_hits_threshold=35)
            cga.featureSettings("ring", local_density_threshold=10.0, local_density_tolerance=15.0, few_hits_threshold=35)
            #
            logger.info("Modified settings:\n%s\n", cga.settings())
        except Exception as e:
            logger.exception("FAILING with %s", str(e))
            self.fail()

    def testGeomAnalFromTargetList(self):
        """"""
        try:
            pL = glob.glob(os.path.join(self.__molFilePath, "*.mol2"))
            logger.info("anal list length %d", len(pL))

            for queryTargetPath in pL:
                _, fn = os.path.split(queryTargetPath)
                queryTargetId, _ = os.path.splitext(fn)
                logger.info("search for %r", queryTargetId)
                atomMap = self.__getAtomMap(queryTargetId, self.__resultPath) if self.__useAtomMap else {}
                cga = CcdcGeomAnal(verbose=self.__verbose, log=self.__lfh)
                rD = cga.anal(queryTargetPath)
                self.__printSummary(queryTargetId, rD, atomMap)
        except Exception as e:
            logger.exception("FAILING with %s", str(e))
            self.fail()

    def __printSummary(self, queryTargetId, rD, atomMap):
        """rD - dictionary of analysis results -
        atomMap - dictionary with mol2 to cc atom name mapping (optional)
        """
        logger.info("\n---------------------------- %s -----------------------", queryTargetId)
        outN = ["bond_outliers", "angle_outliers", "torsion_outliers", "ring_outliers"]
        for ind in outN:
            logger.info("Type: %-20s  Outlier count: %4d", ind, rD[ind])
        #
        outL = ["bond_list", "angle_list", "torsion_list", "ring_list"]
        for ind in outL:
            ll = rD[ind]
            logger.info("Feature: %-20s total count: %4d", ind, len(ll))
            for dD in ll:
                if dD["unusual"]:
                    mappedAtomL = self.__mapAtomNames(dD["atom_labels"], atomMap) if atomMap else dD["atom_labels"]
                    if dD["type"] in ["bond", "angle"]:
                        logger.info("%20s %20s %.4f %.4f %.4f %.4f", dD["atom_labels"], mappedAtomL, dD["value"], dD["mean"], dD["standard_deviation"], dD["z_score"])
                    else:
                        logger.info("%20s %20s %.4f %.4f %.4f %.4f", dD["atom_labels"], mappedAtomL, dD["value"], dD["mean"], dD["standard_deviation"], dD["local_density"])

    def __mapAtomNames(self, atNameL, atomMap):
        mL = []
        try:
            for atName in atNameL:
                mL.append(atomMap[str(atName)])
        except Exception as e:
            logger.exception("Failing for %r and map %r with %s", atNameL, atomMap.items(), str(e))
        return mL

    def __getAtomMap(self, ccId, pth):
        atomMap = {}
        try:
            fp = os.path.join(pth, ccId + "-atom-map.json")
            with open(fp, "r") as ifh:
                atomMap = json.load(ifh)
        except Exception as e:
            logger.exception("FAILING ccid %r path %r with %s", ccId, pth, str(e))

        return atomMap


def suiteGeoAnalTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CcdcGeomAnalTests("testGeomAnalSettings"))
    suiteSelect.addTest(CcdcGeomAnalTests("testGeomAnalFromTargetList"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteGeoAnalTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
