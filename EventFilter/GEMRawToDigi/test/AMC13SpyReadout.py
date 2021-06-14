import FWCore.ParameterSet.Config as cms

process = cms.Process("AMC13SpyReadout")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
process.MessageLogger.cout.threshold = cms.untracked.string('INFO')
process.MessageLogger.debugModules = cms.untracked.vstring('*')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.untracked.int32(-1),
)

process.source = cms.Source(
    "FRDStreamSource",
    fileNames = cms.untracked.vstring("file:run000000_Testing_CERN904_2021-06-08_chunk_0.dat"),
    verifyAdler32 = cms.untracked.bool(False),
    verifyChecksum = cms.untracked.bool(False),
    useL1EventID = cms.untracked.bool(False),
    firstLuminosityBlockForEachRun = cms.untracked.VLuminosityBlockID(*[cms.LuminosityBlockID(0,0)]),
)

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("output_raw.root"),
    outputCommands = cms.untracked.vstring("keep *")
)

process.outpath = cms.EndPath(process.output)
