import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()


process = cms.Process("dumpLHE")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.dummy = cms.EDAnalyzer("WeightLHEAnalyzer",
    src = cms.InputTag("source"),
    histoutputFile= cms.untracked.string(options.outputFile)          
)


process.p = cms.Path(process.dummy)


