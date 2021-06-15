import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')
options.parseArguments()

process = cms.Process("AMC13SpyReadout")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cout.threshold = cms.untracked.string('INFO')
process.MessageLogger.debugModules = cms.untracked.vstring('*')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.untracked.int32(-1),
)

print(options.inputFiles)
process.source = cms.Source(
    "FRDStreamSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    verifyAdler32 = cms.untracked.bool(False),
    verifyChecksum = cms.untracked.bool(False),
    useL1EventID = cms.untracked.bool(False),
    firstLuminosityBlockForEachRun = cms.untracked.VLuminosityBlockID(*[cms.LuminosityBlockID(1,0)]),
)

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("output_raw.root"),
    outputCommands = cms.untracked.vstring("keep *")
)

process.outpath = cms.EndPath(process.output)
