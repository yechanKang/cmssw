import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from Validation.MuonHits.muonSimHitMatcherPSet import muonSimHitMatcherPSet
from Validation.MuonGEMDigis.muonGEMDigiPSet import muonGEMDigiPSet
from Validation.MuonGEMHits.MuonGEMCommonParameters_cfi import GEMValidationCommonParameters

gemStripValidation = DQMEDAnalyzer('GEMStripDigiValidation',
  GEMValidationCommonParameters,
  gemStripDigi = muonGEMDigiPSet.gemUnpackedStripDigi,
  gemSimHit = muonSimHitMatcherPSet.gemSimHit,
)

gemPadValidation = DQMEDAnalyzer('GEMPadDigiValidation',
  GEMValidationCommonParameters,
  gemPadDigi = muonGEMDigiPSet.gemPadDigi,
)

gemClusterValidation = DQMEDAnalyzer('GEMPadDigiClusterValidation',
  GEMValidationCommonParameters,
  gemPadCluster = muonGEMDigiPSet.gemPadCluster,
)

gemCoPadValidation = DQMEDAnalyzer('GEMCoPadDigiValidation',
  GEMValidationCommonParameters,
  gemCoPadDigi = muonGEMDigiPSet.gemCoPadDigi,
)

gemGeometryChecker = DQMEDAnalyzer('GEMCheckGeometry',
  GEMValidationCommonParameters,
)

gemDigiValidation = cms.Sequence(gemStripValidation +
                                 gemPadValidation +
                                 gemClusterValidation +
                                 gemCoPadValidation +
                                 gemGeometryChecker)

from Configuration.Eras.Modifier_run2_GEM_2017_cff import run2_GEM_2017

run2_GEM_2017.toReplaceWith(
    gemDigiValidation,
    gemDigiValidation.copyAndExclude([
        gemPadValidation,
        gemClusterValidation,
        gemCoPadValidation]))
