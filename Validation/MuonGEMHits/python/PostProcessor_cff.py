import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from Validation.MuonGEMHits.MuonGEMCommonParameters_cfi import GEMValidationCommonParameters

gemSimHarvesting = DQMEDHarvester("MuonGEMHitsHarvestor",
    GEMValidationCommonParameters,
)
MuonGEMHitsPostProcessors = cms.Sequence( gemSimHarvesting ) 
