##
# File:    CcdcGeomAnal.py
# Author:  J. Westbrook
# Date:    17-Dec-2020
# Version: 0.001
#
# Updated:
#
##
# pylint: disable=exec-used
"""
Utilities for chemical component geometrical analyis using the local CCDC Python API
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import sys

from ccdc.io import EntryReader
from ccdc import conformer

logger = logging.getLogger(__name__)


class CcdcGeomAnal(object):
    def __init__(self, verbose=True, log=sys.stderr):
        self.__lfh = log
        self.__verbose = verbose
        self.__engine = conformer.GeometryAnalyser()

    def settings(self):
        return self.__engine.settings.summary()

    def featureSettings(self, featureType, **kw):
        """
        Features specific settings -

            featureType = 'bond', 'angle', 'torsion', 'ring'

        Features and values:

            few_hits_threshold: Threshold below which a distribution is considered to have too few hits.

            local_density_threshold: Local density threshold used to classify torsions and rings as unusual.

                                    Note that the local density is irrelevant for bonds and angles.

            local_density_tolerance: The local density tolerance.

            min_obs_exact: Minimum acceptable size of an exact distribution.
                            If there is no distribution containing at least this number of observations the
                            geometry analyser will perform a generalised search according to the criteria specified by other settings.

            min_obs_generalised: Minimum number of observations that the geometry analyser should try to find.
                                If this is 0 then generalised searches will never be performed.
                                Similarly, if generalisation has been turned off this setting will not have an effect.

            min_relevance:  Relevance criterion for a generalised hit to be accepted.

                            The geometry analyser determines how similar a fragment is to the query by calculating
                            a relevance value. The min_relevance setting tells the geometry analyser to accept,
                            in a generalised search, only fragments whose relevance is equal to or greater than this threshold.

            zscore_threshold:  Z-score thresh

            classification_measure: Z-score : Local density, Nearest observation, Mean distance

            classification_measure_threshold:  float


        """
        fTypes = ["bond", "angle", "torsion", "ring"]
        fNames = ["few_hits_threshold", "local_density_threshold", "local_density_tolerance", "min_obs_exact", "min_obs_generalised", "min_relevance", "zscore_threshold"]

        if featureType not in fTypes:
            return False

        for k in kw:
            if k in fNames:
                cm = "e.settings.%s.%s=%s" % (featureType, k, str(kw[k]))
                logger.info("Settings command: %s", cm)
                self.__doExec(self.__engine, cm)

        return True

    def globalSettings(self, **kw):
        """
        Global setting keywords and values:

            rfactor_filter: 0.05,0.075,0.1,any
            generalisation True: False
            organometallic_filter : all,metalorganics_only, organics_only
            solvent_filter: include_solvent,exclude_solvent,only_solvent
            heaviest_element:  <element symbol> (case sensistive)
        """
        fNames = ["rfactor_filter", "generalisation"]
        sNames = ["organometallic_filter", "solvent_filter", "heaviest_element"]

        for k in kw:
            if k in fNames:
                cm = "e.settings.%s=%s" % (k, str(kw[k]))
                logger.info("Settings command: %s", cm)
                self.__doExec(self.__engine, cm)
            if k in sNames:
                cm = "e.settings.%s='%s'" % (k, str(kw[k]))
                logger.info("Settings command: %s", cm)
                self.__doExec(self.__engine, cm)
        return True

    def anal(self, queryTargetPath, normalizeFlag=False):
        """Perform geometrical analysis against the CCDC data source-"""
        retD = {}
        targetStructures = EntryReader(queryTargetPath)

        for e in targetStructures:
            mol = e.molecule
            if normalizeFlag:
                mol.assign_bond_types(which="unknown")
                mol.standardise_aromatic_bonds()
                mol.standardise_delocalised_bonds()
            #
            logger.info("begin analysis - for %s", queryTargetPath)
            gam = self.__engine.analyse_molecule(mol)
            bondOutliers = len([b for b in gam.analysed_bonds if b.unusual and b.enough_hits])
            angleOutliers = len([a for a in gam.analysed_angles if a.unusual and a.enough_hits])
            torsionOutliers = len([t for t in gam.analysed_torsions if t.unusual and t.enough_hits])
            ringOutliers = len([r for r in gam.analysed_rings if r.unusual and r.enough_hits])

            bL = self.__getBondAnalysis(gam)
            aL = self.__getAngleAnalysis(gam)
            tL = self.__getTorsionAnalysis(gam)
            rL = self.__getRingAnalysis(gam)
            retD = {
                "bond_outliers": bondOutliers,
                "angle_outliers": angleOutliers,
                "torsion_outliers": torsionOutliers,
                "ring_outliers": ringOutliers,
                "bond_list": bL,
                "angle_list": aL,
                "torsion_list": tL,
                "ring_list": rL,
            }
        return retD

    def __extractAnalFeatures(self, feature):
        rD = {
            "atom_labels": feature.atom_labels,
            "d_min": feature.d_min,
            "value": feature.value,
            "mean": feature.mean,
            "standard_deviation": feature.standard_deviation,
            "z_score": feature.z_score,
            "type": feature.type,
            "median": feature.median,
            "lower_quartile": feature.lower_quartile,
            "upper_quartile": feature.upper_quartile,
            "minimum": feature.minimum,
            "maximum": feature.maximum,
            "nhits": feature.nhits,
            "unusual": feature.unusual,
            "generalised": feature.generalised,
            "few_hits": feature.few_hits,
            "enough_hits": feature.enough_hits,
            "local_density": feature.local_density,
        }
        return rD

    def __getBondAnalysis(self, mol):
        rL = []
        for feature in mol.analysed_bonds:
            rL.append(self.__extractAnalFeatures(feature))
        return rL

    def __getAngleAnalysis(self, mol):
        rL = []
        for feature in mol.analysed_angles:
            rL.append(self.__extractAnalFeatures(feature))
        return rL

    def __getTorsionAnalysis(self, mol):
        rL = []
        for feature in mol.analysed_torsions:
            rL.append(self.__extractAnalFeatures(feature))
        return rL

    def __getRingAnalysis(self, mol):
        rL = []
        for feature in mol.analysed_rings:
            rL.append(self.__extractAnalFeatures(feature))
        return rL

    def __doExec(self, e, cmd):
        _ = e
        exec(cmd, globals(), locals())
