##
# File:    PdbxChemCompVariantIndex.py
# Author:  jdw
# Date:    20-June-2016
# Version: 0.001
#
# Updates:
#    22-June-2016 jdw add CcdcMatchIndex/Inst() classes -
#    28-July-2017 jdw add PdbxChemCompIndex and PdbxChemCompIndexInst classes -
#    28-July-2017 jdw change variant_id to target_id in the ccdc results index
##
"""
Accessor methods for index of chemical component variants

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import sys
import os

from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class IndexBase(object):
    """Base category iterator class."""

    def __init__(self, indexFilePath, func, verbose=True):
        self._verbose = verbose
        self._debug = False
        self._rL = []
        self.__func = func
        self.__fmt = "json"
        self._indexFilePath = indexFilePath
        self._indexPath, self._indexFileName = os.path.split(self._indexFilePath)

    def get(self, index=0):
        try:
            return self._rL[index]
        except Exception:
            return []

    def __iter__(self):
        return self.forward()

    def forward(self):
        # Forward generator
        currentRow = 0
        while currentRow < len(self._rL):
            row = self._rL[currentRow]
            currentRow += 1
            yield self.__func(row)

    def reverse(self):
        # The reverse generator
        currentRow = len(self._rL)
        while currentRow > 0:
            currentRow -= 1
            yield self.__func(self._rL[currentRow])

    def clear(self):
        self._rL = []
        return True

    def writeIndex(self):
        try:
            mU = MarshalUtil()
            ok = mU.doExport(self._indexFilePath, self._rL, fmt=self.__fmt, indent=3)
            return ok
        except Exception as e:
            logger.error("Failing with %s", str(e))

        return False

    def readIndex(self):
        try:
            mU = MarshalUtil()
            if not mU.exists(self._indexFilePath):
                return False
            indexObj = mU.doImport(self._indexFilePath, fmt=self.__fmt)
            if indexObj is not None and len(indexObj) > 0:
                self._rL.extend(indexObj)
            return True
        except Exception as e:
            logger.error("Failing with %s", str(e))

        return False

    def load(self, rowList):
        self._rL.extend(rowList)

    def dump(self):
        for ii, dD in enumerate(self._rL):
            logger.info("%4d: %r", ii, dD)
        logger.info("Completed dump")


class PdbxChemCompIndex(IndexBase):
    def __init__(self, indexFilePath="../results/test_comp/index.pic", verbose=True):
        obj = PdbxChemCompIndexInst(None, verbose=verbose)
        super(PdbxChemCompIndex, self).__init__(indexFilePath, obj.set, verbose)
        #
        self.readIndex()
        logger.info("Chemical component index length %d", len(self._rL))
        if self._debug:
            self.dump()

    def makeChemCompPathIndex(self):
        """Index is a list of tuples of (ccID, InChIKey)

        Return index  d[cc_id] = [ {<component details>}, {<component details},...]  sorted by component index.
        """
        iD = {}
        for dd in self._rL:
            ccId = dd["cc_id"]
            if ccId not in iD:
                iD[ccId] = []
            iD[ccId].append(dd)

        pD = {}
        for ky in iD:
            pList = []
            dL = iD[ky]
            for dD in dL:
                dD["cif_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".cif")
                dD["cif_mapped_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + "-mapped.cif")
                dD["mol2_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".mol2")
                dD["mol_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".sdf")
                dD["sort_key"] = dD["cc_id"]
                dI = PdbxChemCompIndexInst(dObj=dD)
                pList.append(dI)
            tL = sorted(pList, key=lambda x: x.getSortKey())
            pD[ky] = tL
        return pD

    def makeChemCompPathList(self):
        """Index is a list of tuples of (ccID, InChIKey)

        Return List   pL = [(cc_id,vObj), (cc_id,vObj), ... ]
        """
        iD = {}
        for dd in self._rL:
            ccId = dd["cc_id"]
            if ccId not in iD:
                iD[ccId] = []
            iD[ccId].append(dd)
        #
        #
        pL = []
        for ky in iD:
            pList = []
            dL = iD[ky]
            for dD in dL:
                dD["cif_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".cif")
                dD["cif_mapped_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + "-mapped.cif")
                dD["mol2_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".mol2")
                dD["mol_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".sdf")
                dD["sort_key"] = dD["cc_id"]
                dI = PdbxChemCompIndexInst(dObj=dD)
                pList.append(dI)
            tL = sorted(pList, key=lambda x: x.getSortKey())
            for tt in tL:
                pL.append((ky, tt))

        return pL

    def makeChemCompPathPriorityList(self, priorityListPath="../data/ccid_target_list.txt", filterListPath=None):
        """Index is a list of tuples of (ccID, InChIKey)

        Return List   pL = [(cc_id,vObj), (cc_id,vObj), ... ]
        """
        iD = {}
        for dd in self._rL:
            ccId = dd["cc_id"]
            if ccId not in iD:
                iD[ccId] = []
            iD[ccId].append(dd)
        #
        ##
        #
        filterList = []
        filterD = {}
        idList = self.__getTargetIdList(fn=priorityListPath)
        logger.info("Priority list length %d", len(idList))
        try:
            if filterListPath is not None and os.access(filterListPath, os.R_OK):
                filterList = self.__getTargetIdList(fn=filterListPath)
                logger.info("Filterlist length %d", len(filterList))
                for tid in filterList:
                    filterD[tid] = tid
            else:
                logger.info("No filterlist")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        fCount = 0
        pL = []
        for ky in idList:
            if ky not in iD:
                continue
            pList = []
            dL = iD[ky]
            for dD in dL:
                if dD["cc_id"] in filterD:
                    fCount += 1
                    continue
                dD["cif_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".cif")
                dD["cif_mapped_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + "-mapped.cif")
                dD["mol2_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".mol2")
                dD["mol_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_id"] + ".sdf")
                dD["sort_key"] = dD["cc_id"]
                dI = PdbxChemCompIndexInst(dObj=dD)
                pList.append(dI)
            tL = sorted(pList, key=lambda x: x.getSortKey())
            for tt in tL:
                pL.append((ky, tt))
        #
        logger.info("Filtered length %d", fCount)
        return pL

    def makeChemCompDepictionList(self):
        """Make index with meta details for depiction functions -

        Return index d[cc_id] = [(ccId, PDBx path, cc title), ... ]

        """
        rD = {}
        pD = self.makeChemCompPathIndex()
        for ky in pD:
            if ky not in rD:
                rD[ky] = []
            for pp in pD[ky]:
                rD[ky].append((pp.getComponentId(), pp.getCifMappedPath(), pp.getComponentId()))
        return rD

    def __getTargetIdList(self, fn="../data/ccid_target_list.txt"):
        idList = []
        try:
            fh = open(fn, "r")
            for line in fh:
                idList.append(str(line[:-1]).strip())
            fh.close()
        except Exception:
            pass
        return idList


class PdbxChemCompVariantIndex(IndexBase):
    def __init__(self, indexFilePath="../results/test_prot_taut/index.pic", verbose=True):
        obj = PdbxChemCompVariantIndexInst(None, verbose=verbose)
        super(PdbxChemCompVariantIndex, self).__init__(indexFilePath, obj.set, verbose)
        #
        self.readIndex()
        logger.info("Variant index length %d", len(self._rL))
        if self._debug:
            self.dump()

    def makeChemCompPathIndex(self):
        """Index is a list of tuples of (Base ccID, Protomer index (1-n), Tautomer index (1-n), Variant ccIdU, InChIKey)

        Return index  d[cc_id] = [ {<variant details>}, {<variant details},...]  sorted by protomer then tautomer index.
        """
        iD = {}
        for dd in self._rL:
            ccId = dd["cc_id"]
            if ccId not in iD:
                iD[ccId] = []
            iD[ccId].append(dd)

        pD = {}
        for ky in iD:
            pList = []
            dL = iD[ky]
            for dD in dL:
                dD["cif_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".cif")
                dD["cif_mapped_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + "-mapped.cif")
                dD["mol2_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".mol2")
                dD["mol_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".sdf")
                dD["sort_key"] = "%03d%03d" % (int(dD["protomer_index"]), int(dD["tautomer_index"]))
                dI = PdbxChemCompVariantIndexInst(dObj=dD)
                pList.append(dI)
            tL = sorted(pList, key=lambda x: x.getSortKey())
            pD[ky] = tL
        return pD

    def makeChemCompPathList(self):
        """Index is a list of tuples of (Base ccID, Protomer index (1-n), Tautomer index (1-n), Variant ccIdU, InChIKey)

        Return List   pL = [(cc_id,vObj), (cc_id,vObj), ... ]
        """
        iD = {}
        for dd in self._rL:
            ccId = dd["cc_id"]
            if ccId not in iD:
                iD[ccId] = []
            iD[ccId].append(dd)
        #
        #
        pL = []
        for ky in iD:
            pList = []
            dL = iD[ky]
            for dD in dL:
                dD["cif_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".cif")
                dD["cif_mapped_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + "-mapped.cif")
                dD["mol2_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".mol2")
                dD["mol_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".sdf")
                dD["sort_key"] = "%03d%03d" % (int(dD["protomer_index"]), int(dD["tautomer_index"]))
                dI = PdbxChemCompVariantIndexInst(dObj=dD)
                pList.append(dI)
            tL = sorted(pList, key=lambda x: x.getSortKey())
            for tt in tL:
                pL.append((ky, tt))

        return pL

    def makeChemCompPathPriorityList(self, priorityListPath="../data/ccid_target_list.txt", filterListPath=None):
        """Index is a list of tuples of (Base ccID, Protomer index (1-n), Tautomer index (1-n), Variant ccIdU, InChIKey)

        Return List   pL = [(cc_id,vObj), (cc_id,vObj), ... ]
        """
        iD = {}
        for dd in self._rL:
            ccId = dd["cc_id"]
            if ccId not in iD:
                iD[ccId] = []
            iD[ccId].append(dd)
        #
        ##
        #
        filterList = []
        filterD = {}
        idList = self.__getTargetIdList(fn=priorityListPath)
        logger.info("Priority list length %d", len(idList))
        try:
            if filterListPath is not None and os.access(filterListPath, os.R_OK):
                filterList = self.__getTargetIdList(fn=filterListPath)
                logger.info("Filterlist length %d", len(filterList))
                for tid in filterList:
                    filterD[tid] = tid
            else:
                logger.info("No filterlist")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        fCount = 0
        pL = []
        for ky in idList:
            if ky not in iD:
                continue
            pList = []
            dL = iD[ky]
            for dD in dL:
                if dD["cc_variant_id"] in filterD:
                    fCount += 1
                    continue
                dD["cif_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".cif")
                dD["cif_mapped_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + "-mapped.cif")
                dD["mol2_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".mol2")
                dD["mol_variant_path"] = os.path.join(self._indexPath, dD["cc_id"], dD["cc_variant_id"] + ".sdf")
                dD["sort_key"] = "%03d%03d" % (int(dD["protomer_index"]), int(dD["tautomer_index"]))
                dI = PdbxChemCompVariantIndexInst(dObj=dD)
                pList.append(dI)
            tL = sorted(pList, key=lambda x: x.getSortKey())
            for tt in tL:
                pL.append((ky, tt))
        #
        logger.info("Filtered length %d", fCount)
        return pL

    def makeChemCompDepictionList(self):
        """Make index with meta details for depiction functions -

        Return index d[cc_id] = [(variantId, PDBx variant path, variant title), ... ]

        """
        rD = {}
        pD = self.makeChemCompPathIndex()
        for ky in pD:
            if ky not in rD:
                rD[ky] = []
            for pp in pD[ky]:
                rD[ky].append((pp.getVariantId(), pp.getCifMappedVariantPath(), pp.getVariantId()))
        return rD

    def __getTargetIdList(self, fn="../data/ccid_target_list.txt"):
        idList = []
        try:
            fh = open(fn, "r")
            for line in fh:
                idList.append(str(line[:-1]).strip())
            fh.close()
        except Exception:
            pass
        return idList


class PdbxChemCompIndexInst(object):
    """Accessor methods for chemical component instance index elements ."""

    def __init__(self, dObj=None, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        self._lfh = log
        self.__dObj = dObj if dObj is not None else {}

    def _getAttribute(self, name):
        try:
            return self.__dObj[name]
        except Exception:
            return None

    def set(self, dObj=None):
        self.__dObj = dObj
        return self

    def get(self):
        return self.__dObj

    def _setAttribute(self, name, value):
        try:
            self.__dObj[name] = value
        except Exception:
            return False

    def getChemCompId(self):
        return self._getAttribute("cc_id")

    def getSortKey(self):
        return self._getAttribute("sort_key")

    def getCifPath(self):
        return self._getAttribute("cif_path")

    def getCifMappedPath(self):
        return self._getAttribute("cif_mapped_path")

    def getMol2Path(self):
        return self._getAttribute("mol2_path")

    def getMolPath(self):
        return self._getAttribute("mol_path")

    def getInChIKey(self):
        return self._getAttribute("inchikey")

    #
    def setChemCompId(self, ccId):
        self.__dObj["cc_id"] = ccId

    def setInChIKey(self, ky):
        self.__dObj["inchikey"] = ky

    def setSmiles(self, smi):
        self.__dObj["smiles"] = smi

    def setStereoSmiles(self, smi):
        self.__dObj["stereo_smiles"] = smi

    def __repr__(self):
        return str(self.__dObj)

    def __str__(self):
        return str(self.__dObj)


class PdbxChemCompVariantIndexInst(PdbxChemCompIndexInst):
    """Accessor methods for chemical component variant instance index elements ."""

    def __init__(self, dObj=None, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        self._lfh = log
        super(PdbxChemCompVariantIndexInst, self).__init__(dObj, verbose, log)

    def getVariantId(self):
        return self._getAttribute("cc_variant_id")

    def getCifVariantPath(self):
        return self._getAttribute("cif_variant_path")

    def getCifMappedVariantPath(self):
        return self._getAttribute("cif_mapped_variant_path")

    def getMol2VariantPath(self):
        return self._getAttribute("mol2_variant_path")

    def getMolVariantPath(self):
        return self._getAttribute("mol_variant_path")

    #
    def setTautomerIndex(self, i):
        self._setAttribute("tautomer_index", i)

    def setProtomerIndex(self, i):
        self._setAttribute("protomer_index", i)

    def setChemCompVariantId(self, ccId):
        self._setAttribute("cc_variant_id", ccId)


class CcdcMatchIndex(IndexBase):
    def __init__(self, indexFilePath=None, verbose=True):
        obj = CcdcMatchIndexInst(None, verbose=verbose)
        super(CcdcMatchIndex, self).__init__(indexFilePath, obj.set, verbose)
        #
        self.readIndex()
        if self._debug:
            self.dump()


class CcdcMatchIndexInst(object):
    """Accessor methods for index elements ."""

    def __init__(self, dObj=None, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        self._lfh = log
        self.__dObj = dObj if dObj is not None else {}

    def __getAttribute(self, name):
        try:
            return self.__dObj[name]
        except Exception:
            return None

    def set(self, dObj=None):
        self.__dObj = dObj
        return self

    def get(self):
        return self.__dObj

    def getIdentifier(self):
        return self.__getAttribute("identifier")

    def getChemicalName(self):
        return self.__getAttribute("chemical_name")

    def getTargetId(self):
        return self.__getAttribute("target_id")

    def getMatchType(self):
        return self.__getAttribute("match_type")

    def getMatchNumber(self):
        return self.__getAttribute("match_number")

    def getMol2Path(self):
        return self.__getAttribute("mol2_file_path")

    def getMolPath(self):
        return self.__getAttribute("mol_file_path")

    #

    def getRFactor(self):
        return self.__getAttribute("r_factor")

    def getTemperature(self):
        return self.__getAttribute("temperature")

    def getRadiationSource(self):
        return self.__getAttribute("radiation_source")

    def getCitationDOI(self):
        return self.__getAttribute("doi")

    def getMatchedAtomLength(self):
        return self.__getAttribute("match_atoms")

    def getHasDisorder(self):
        return self.__getAttribute("has_disorder")

    #

    def getSimilarityScore(self):
        return self.__getAttribute("similarity")

    def setIdentifier(self, v):
        self.__dObj["identifier"] = v

    def setChemicalName(self, name):
        self.__dObj["chemical_name"] = name

    def setTargetId(self, ccId):
        self.__dObj["target_id"] = ccId

    def setMatchType(self, v):
        self.__dObj["match_type"] = v

    def setMatchNumber(self, v):
        self.__dObj["match_number"] = v

    def setMolPath(self, fp):
        self.__dObj["mol_file_path"] = fp

    def setMol2Path(self, fp):
        self.__dObj["mol2_file_path"] = fp

    def setRFactor(self, rv):
        self.__dObj["r_factor"] = rv

    def setTemperature(self, tv):
        self.__dObj["temperature"] = tv

    def setRadiationSource(self, sv):
        self.__dObj["radiation_source"] = sv

    def setSimilarityScore(self, sv):
        self.__dObj["similarity_score"] = sv

    def setMatchedAtomLength(self, lv):
        self.__dObj["match_atoms"] = lv

    def setHasDisorder(self, lv):
        self.__dObj["has_disorder"] = lv

    def setCitationDOI(self, doi):
        self.__dObj["doi"] = doi

    def __repr__(self):
        return str(self.__dObj)

    def __str__(self):
        return str(self.__dObj)
