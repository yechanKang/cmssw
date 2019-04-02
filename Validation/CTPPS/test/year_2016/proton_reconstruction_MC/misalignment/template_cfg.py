import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process('CTPPSProtonReconstructionTest', eras.ctpps_2016)

# minimal logger settings
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# number of events
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(280000)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# particle-data table
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# random seeds
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    sourceSeed = cms.PSet(initialSeed =cms.untracked.uint32(98765)),
    generator = cms.PSet(initialSeed = cms.untracked.uint32(98766)),
    beamDivergenceVtxGenerator = cms.PSet(initialSeed =cms.untracked.uint32(3849))
)

# geometry
process.load("Geometry.VeryForwardGeometry.geometryRPFromDD_2017_cfi")
del(process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles[-1])
process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.append("Validation/CTPPS/test/year_2016/RP_Dist_Beam_Cent.xml")

process.load("CalibPPS.ESProducers.ctppsRPAlignmentCorrectionsDataESSourceXML_cfi")

# beam-smearing settings
process.load("IOMC.EventVertexGenerators.beamDivergenceVtxGenerator_cfi")

# simulation-level settings
from Validation.CTPPS.year_2016.simulation_levels_cff import *

process.ctppsBeamParametersESSource = ctppsBeamParametersESSource
process.generator = generator
process.ctppsDirectProtonSimulation = ctppsDirectProtonSimulation

# strips reco
process.load('RecoCTPPS.TotemRPLocal.totemRPUVPatternFinder_cfi')
process.totemRPUVPatternFinder.tagRecHit = cms.InputTag('ctppsDirectProtonSimulation')

process.load('RecoCTPPS.TotemRPLocal.totemRPLocalTrackFitter_cfi')

# common reco: lite track production
process.load('RecoCTPPS.TotemRPLocal.ctppsLocalTrackLiteProducer_cff')
process.ctppsLocalTrackLiteProducer.includeDiamonds = False
process.ctppsLocalTrackLiteProducer.includePixels = False

# proton reconstruction
process.load("RecoCTPPS.ProtonReconstruction.ctppsProtons_cff")
process.ctppsProtons.tagLocalTrackLite = cms.InputTag('ctppsLocalTrackLiteProducer')
#process.ctppsProtons.fitVtxY = False

# reconstruction validation
process.ctppsProtonReconstructionSimulationValidator = cms.EDAnalyzer("CTPPSProtonReconstructionSimulationValidator",
  tagHepMCBeforeSmearing = cms.InputTag("generator", "unsmeared"),
  tagHepMCAfterSmearing = cms.InputTag("beamDivergenceVtxGenerator"),
  tagRecoProtonsSingleRP = cms.InputTag("ctppsProtons", "singleRP"),
  tagRecoProtonsMultiRP = cms.InputTag("ctppsProtons", "multiRP"),

  lhcInfoLabel = process.ctppsLHCInfoESSource_2016.label,  # TODO: replace with "" once data loaded from DB

  outputFile = cms.string("")
)

# processing path
process.p = cms.Path(
    process.generator
    * process.beamDivergenceVtxGenerator
    * process.ctppsDirectProtonSimulation

    * process.totemRPUVPatternFinder
    * process.totemRPLocalTrackFitter
    * process.ctppsLocalTrackLiteProducer

    * process.ctppsProtons

    * process.ctppsProtonReconstructionSimulationValidator
)

#----------

SetLargeTheta()
SetLevel4()

