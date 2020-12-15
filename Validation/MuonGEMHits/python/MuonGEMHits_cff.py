import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from Validation.MuonHits.muonSimHitMatcherPSet import muonSimHitMatcherPSet
from Validation.MuonGEMHits.MuonGEMCommonParameters_cfi import GEMValidationCommonParameters

gemSimHitValidation = DQMEDAnalyzer('GEMSimHitValidation',
    GEMValidationCommonParameters,
    gemSimHit = muonSimHitMatcherPSet.gemSimHit,
    TOFRange = cms.untracked.vdouble(10, 30, # GE0
                                     16, 24, # GEM11
                                     24, 32), # GE21
)

gemSimValidation = cms.Sequence(gemSimHitValidation)
