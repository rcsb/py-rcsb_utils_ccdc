##
# File: CcdcSearchExec.py
# Date: 15-Dec-2020  jdw
#
#  Execution wrapper  --  for CCDC search operations (wraps up the py37 environment)
#
#  Updates:
#   15-Jan-2021 jdw add option to export search hit list.
#
##
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import argparse
import logging
import os
import sys

from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


def main():
    parser = argparse.ArgumentParser()
    #
    parser.add_argument("--mol_list_path", default=None, help="Molecule file list path")
    parser.add_argument("--result_path", default=None, help="Molecule file list path")
    parser.add_argument("--search_type", default=None, help="Search type (similarity|substructure)")
    parser.add_argument("--start_record", default=None, help="Starting record")
    parser.add_argument("--end_record", default=None, help="End record")
    parser.add_argument("--csdhome", default=None, help="Path to the CSD release (path to CSD_202x)")
    parser.add_argument("--python_lib_path", default=None, help="Path to Python library")
    parser.add_argument("--python_version", default=None, help="Python library version (default: 3.7)")
    parser.add_argument("--hit_list_path", default=None, help="Path to list of molecule identifers with search results")
    #
    args = parser.parse_args()
    #
    try:
        pyLib = args.python_lib_path if args.python_lib_path else os.path.join(os.environ["PYENV_ROOT"], "versions", "3.7.9", "lib")
        pyVer = args.python_version if args.python_version else "3.7"

        csdHome = args.csdhome
        molFilePath = args.mol_list_path
        resultPath = args.result_path
        searchType = args.search_type
        startRecord = args.start_record
        endRecord = args.end_record
        hitListPath = args.hit_list_path
    except Exception as e:
        logger.exception("Argument processing problem %s", str(e))
        parser.print_help(sys.stderr)
        exit(1)
    #
    try:
        os.environ["CSDHOME"] = csdHome
        os.environ["LD_LIBRARY_PATH"] = "%s:%s/python%s/site-packages/ccdc/_lib:$LD_LIBRARY_PATH" % (pyLib, pyLib, pyVer)
        os.environ["DYLD_LIBRARY_PATH"] = "%s/python%s/site-packages/ccdc/_lib" % (pyLib, pyVer)
        os.environ["DYLD_FRAMEWORK_PATH"] = "%s/python%s/site-packages/ccdc/_lib" % (pyLib, pyVer)

        logger.info("Using CSDHOME %s", os.environ["CSDHOME"])
        logger.info("Using DYLD_LIBRARY_PATH %s", os.environ["DYLD_LIBRARY_PATH"])
        logger.info("Using DYLD_FRAMEWORK_PATH %s", os.environ["DYLD_FRAMEWORK_PATH"])

        from rcsb.utils.ccdc.CcdcSearch import CcdcSearch  # pylint: disable=import-outside-toplevel

        ccdcS = CcdcSearch(verbose=True)
        pL = ccdcS.getList(molFilePath, startRecord=startRecord, endRecord=endRecord)
        logger.info("Search file %s record length %r", molFilePath, len(pL) if pL else [])
        #
        hitL = []
        for ii, queryTargetPath in enumerate(pL, 1):
            _, fn = os.path.split(queryTargetPath)
            queryTargetId, _ = os.path.splitext(fn)
            #
            logger.info("(%d/%d) Start search for %r %r", ii, len(pL), queryTargetId, queryTargetPath)
            numHits = ccdcS.search(queryTargetId, queryTargetPath, resultPath, searchType=searchType)
            if numHits:
                hitL.append(queryTargetId)
        logger.info("%d searches completed - matched %d", len(pL), len(hitL))
        if hitListPath:
            mU = MarshalUtil()
            ok = mU.doExport(hitListPath, hitL, fmt="list")
            logger.info("Wrote hit list (%r) to %s", ok, hitListPath)
    except Exception as e:
        logger.exception("Failing with %s", str(e))

    # ----------------------- - ----------------------- - ----------------------- - ----------------------- - ----------------------- -


if __name__ == "__main__":
    main()
