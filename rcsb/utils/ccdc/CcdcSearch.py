##
#
# File:    CcdcSearch.py
# Author:  J. Westbrook
# Date:    28-July-2017
# Version: 0.001
#
# Updated:
#   21-Jun-2016   jdw  export additional metadata with each matching structure  -
#   22-Jun-2016   jdw  refactor with general index class CcdcMatchIndex -
#   28-Jul-2017   jdw  Generalize to CcdcSearch.py
#   28-Jul-2017   jdw  remove parentId ---
#
##
"""
Chemical component search against the local CCDC using Python API

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

# pylint: disable=not-context-manager

import copy
import logging
import time
import os

from ccdc.io import EntryReader, MoleculeWriter, csd_version, csd_directory
from ccdc.search import SimilaritySearch, TextNumericSearch, MoleculeSubstructure, SubstructureSearch, SMARTSSubstructure

from rcsb.utils.io.IndexUtils import CcdcMatchIndex, CcdcMatchIndexInst
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class CcdcSearch(object):
    def __init__(self, verbose=True, similarityThreshold=0.95, rValueMaxPercent=10.0):
        self.__verbose = verbose
        self.__similarityThreshold = similarityThreshold
        self.__rValueMaxPercent = rValueMaxPercent

    def getList(self, listPath, startRecord=None, endRecord=None):
        rL = []
        try:
            mU = MarshalUtil()
            logger.debug("Reading path list %r", listPath)
            pL = mU.doImport(listPath, fmt="list")
            logger.debug("Reading path list %r (%r)", listPath, len(pL) if pL else None)
            logger.debug("path list %r", pL)
            if startRecord and endRecord:
                rL = pL[startRecord - 1 : endRecord]
            else:
                # take the full list
                rL = pL
        except Exception as e:
            logger.exception("Failing for %r with %s", listPath, str(e))
        #
        return rL

    def search(self, queryTargetId, queryTargetPath, resultPath, normalizeFlag=True, maxHits=50, searchType="similarity", suppressMetals=False):
        """Search the CCDC database for similar or substructure matches for the input query molecule.

        Args:
            queryTargetId (str): query identifier
            queryTargetPath (str): path to the query molfile (mol, sdf, mol2)
            resultPath (str): output path to match results
            normalizeFlag (bool, optional): do standard perceptions on matching molecules. Defaults to True.
            maxHits (int, optional): maximum number of matches to return. Defaults to 50.
            searchType (str, optional): search mode (substructure, similarity). Defaults to "similarity".
            suppressMetals (bool, optional): filter structures containing metals. Defaults to False.

        Returns:
            (int): number of matches
        """

        mU = MarshalUtil()
        logger.info("Start search for target %s path %s result path %s", queryTargetId, queryTargetPath, resultPath)
        #
        summaryList = []
        #
        targetDirPath = os.path.dirname(queryTargetPath)
        cifTargetPath = os.path.join(targetDirPath, queryTargetId + ".cif")

        #
        targetStructures = EntryReader(queryTargetPath)
        dirPath = os.path.join(resultPath, queryTargetId)
        numHits = 0
        for ii, e in enumerate(targetStructures, 1):
            numHits = 0
            startTime = time.time()
            targetMol = e.molecule
            if normalizeFlag:
                targetMol.assign_bond_types(which="unknown")
                targetMol.standardise_aromatic_bonds()
                targetMol.standardise_delocalised_bonds()
            #
            logger.info("(%d) begin %s search - query id %s", ii, searchType, queryTargetId)
            if searchType == "similarity":
                hits = self.__similaritySearch(targetMol, suppressMetals=suppressMetals)
            elif searchType == "substructure":
                hits = self.__moleculeSubstructureSearch(targetMol, suppressMetals=suppressMetals)
            else:
                hits = []
            logger.info("(%d) completed search query id %s in %.3f seconds", ii, queryTargetId, time.time() - startTime)

            if hits:
                numHits += len(hits)
                logger.info("(%d) search for %s matched %d: %r", ii, queryTargetId, numHits, [targetHit.identifier for targetHit in hits])

                #
                for targetHit in hits[:maxHits]:
                    #
                    hI = CcdcMatchIndexInst()
                    hI.setCsdVersion(csd_version())
                    hI.setCsdDirectory(csd_directory())
                    hI.setTargetId(queryTargetId)
                    hI.setTargetPath(queryTargetPath)
                    if mU.exists(cifTargetPath):
                        hI.setTargetCcPath(cifTargetPath)
                    hI.setIdentifier(targetHit.identifier)
                    hI.setMatchType(searchType)
                    try:
                        hI.setRFactor(targetHit.entry.r_factor)
                        hI.setChemicalName(targetHit.entry.chemical_name)
                        hI.setTemperature(targetHit.entry.temperature)
                        hI.setRadiationSource(targetHit.entry.radiation_source)
                        cit = targetHit.entry.publication
                        if cit.doi is not None:
                            hI.setCitationDOI(cit.doi)
                        if searchType == "similarity":
                            hI.setSimilarityScore(targetHit.similarity)
                        elif searchType == "substructure":
                            hI.setMatchedAtomLength(len(targetHit.match_atoms()))
                    except Exception as e:
                        logger.exception("Failing with %s", str(e))
                        #
                    #
                    mU.mkdir(dirPath)
                    mol2L = []
                    if searchType == "substructure":
                        for jj, mc in enumerate(targetHit.match_components(), 1):
                            fp = os.path.join(dirPath, queryTargetId + "_" + targetHit.identifier + "_%03d" % jj + ".mol2")
                            mol2L.append(fp)
                            with MoleculeWriter(fp) as ofh:
                                ofh.write(mc)

                            #
                            fp = os.path.join(dirPath, queryTargetId + "_" + targetHit.identifier + "_%03d" % jj + ".sdf")
                            with MoleculeWriter(fp) as ofh:
                                ofh.write(mc)
                        #
                        #  Check for multiple generated result files -
                        #
                        for jj, fp in enumerate(mol2L, 1):
                            logger.debug("(%d) adding component fp %s", jj, fp)
                            hI.setMatchNumber(jj)
                            hI.setMol2Path(fp)
                            tt = fp[:-4] + "sdf"
                            hI.setMolPath(tt)
                            summaryList.append(copy.deepcopy(hI.get()))
                            #
                    else:
                        hI.setMatchNumber(1)
                        summaryList.append(copy.deepcopy(hI.get()))
            else:
                logger.info("(%d) search for %s returns no matches", ii, targetMol.identifier)
                hits = None
        #
        if numHits > 0:
            mU.mkdir(dirPath)
            fp = os.path.join(dirPath, queryTargetId + "-index.json")
            cmI = CcdcMatchIndex(indexFilePath=fp, verbose=self.__verbose)
            cmI.load(summaryList)
            cmI.writeIndex()

        return numHits

    def searchSmarts(self, queryTargetId, smarts, resultPath, maxHits=50, suppressMetals=False):
        """Search the CCDC database for substructure matches for the input SMARTS pattern.

        Args:
            queryTargetId (str): query identifier
            smarts (srt): smarts search pattern (NON STEREO)
            resultPath (str): output path to match results
            maxHits (int, optional): maximum number of matches to return. Defaults to 50.
            suppressMetals (bool, optional): filter structures containing metals. Defaults to False.

        Returns:
            (int): number of matches
        """

        mU = MarshalUtil()
        logger.info("Start smarts search for target %s result path %s", queryTargetId, resultPath)
        #
        ii = 1
        searchType = "substructure"
        summaryList = []
        dirPath = os.path.join(resultPath, queryTargetId)
        numHits = 0
        startTime = time.time()
        logger.info("(%d) begin %s search - query id %s", ii, searchType, queryTargetId)

        if searchType == "substructure":
            hits = self.__smartsSubstructureSearch(smarts, suppressMetals=suppressMetals)
        else:
            hits = []
        logger.info("(%d) completed search query id %s in %.3f seconds", ii, queryTargetId, time.time() - startTime)

        if hits:
            numHits += len(hits)
            logger.info("(%d) search for %s matched %d: %r", ii, queryTargetId, numHits, [targetHit.identifier for targetHit in hits])

            #
            for targetHit in hits[:maxHits]:
                #
                hI = CcdcMatchIndexInst()
                hI.setTargetId(queryTargetId)
                hI.setIdentifier(targetHit.identifier)
                hI.setMatchType(searchType)
                try:
                    hI.setRFactor(targetHit.entry.r_factor)
                    hI.setChemicalName(targetHit.entry.chemical_name)
                    hI.setTemperature(targetHit.entry.temperature)
                    hI.setRadiationSource(targetHit.entry.radiation_source)
                    cit = targetHit.entry.publication
                    if cit.doi is not None:
                        hI.setCitationDOI(cit.doi)
                    if searchType == "similarity":
                        hI.setSimilarityScore(targetHit.similarity)
                    elif searchType == "substructure":
                        hI.setMatchedAtomLength(len(targetHit.match_atoms()))
                except Exception as e:
                    logger.exception("Failing with %s", str(e))
                    #
                #
                mU.mkdir(dirPath)
                mol2L = []
                for jj, mc in enumerate(targetHit.molecule.components, 1):
                    fp = os.path.join(dirPath, queryTargetId + "_" + targetHit.identifier + "_%03d" % jj + ".mol2")
                    mol2L.append(fp)
                    with MoleculeWriter(fp) as ofh:
                        ofh.write(mc)

                    #
                    fp = os.path.join(dirPath, queryTargetId + "_" + targetHit.identifier + "_%03d" % jj + ".sdf")
                    with MoleculeWriter(fp) as ofh:
                        ofh.write(mc)
                #
                #  Check for multiple generated result files -
                #
                for jj, fp in enumerate(mol2L, 1):
                    logger.debug("(%d) adding component fp %s", jj, fp)
                    hI.setMatchNumber(jj)
                    hI.setMol2Path(fp)
                    tt = fp[:-4] + "sdf"
                    hI.setMolPath(tt)
                    summaryList.append(copy.deepcopy(hI.get()))
                    #
        else:
            logger.info("(%d) se sarch for %s returns no matches", ii, queryTargetId)
            hits = None
        #
        if numHits > 0:
            mU.mkdir(dirPath)
            fp = os.path.join(dirPath, queryTargetId + "-index.json")
            cmI = CcdcMatchIndex(indexFilePath=fp, verbose=self.__verbose)
            cmI.load(summaryList)
            cmI.writeIndex()

        return numHits

    def __moleculeSubstructureSearch(self, aMol, suppressMetals=False):
        ms = MoleculeSubstructure(aMol)
        search = SubstructureSearch()
        search.add_substructure(ms)
        search.settings.has_3d_coordinates = True
        search.settings.no_disorder = True
        if suppressMetals:
            search.settings.only_organic = True
            search.settings.no_metals = True
        search.settings.max_r_factor = self.__rValueMaxPercent
        hits = search.search(max_hits_per_structure=1)
        return hits

    def __smartsSubstructureSearch(self, smarts, suppressMetals=False):
        ss = SMARTSSubstructure(smarts)
        search = SubstructureSearch()
        search.add_substructure(ss)
        search.settings.has_3d_coordinates = True
        search.settings.no_disorder = True
        if suppressMetals:
            search.settings.no_metals = True
            search.settings.only_organic = True
        search.settings.max_r_factor = self.__rValueMaxPercent
        hits = search.search(max_hits_per_structure=1)
        return hits

    def __similaritySearch(self, aMol, suppressMetals=False):
        # the similarity threshold is a score from 0 to 1 of how 'similar' the structures will be to the input molecule
        search = SimilaritySearch(aMol, threshold=self.__similarityThreshold)
        search.settings.has_3d_coordinates = True
        search.settings.no_disorder = True
        if suppressMetals:
            search.settings.only_organic = True
            search.settings.no_metals = True
        search.settings.max_r_factor = self.__rValueMaxPercent
        hits = search.search(max_hits_per_structure=1)
        return hits

    def __textSearch(self, aMol):
        search = TextNumericSearch()
        search.add_compound_name(aMol.identifier, mode="anywhere", ignore_non_alpha_num=True)
        # search.add_synonym(m.identifier) - no need to add synonym separate (should automatically search both name and synonym)
        textHits = search.search(max_hits_per_structure=1)
        return textHits
