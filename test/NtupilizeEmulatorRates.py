import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.outputFile = "NtupilizeEmulatorRates.root"
options.parseArguments()
process = cms.Process('EmulatorNtupilizer')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')\

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#Max Events
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(options.maxEvents)
	#input = cms.untracked.int32(10000)
)

#Input
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles)
#        fileNames = cms.untracked.vstring("/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/C2D9064D-FB7F-E311-BFCC-003048679228.root")
)

#Output
#options.outputFile = "L1EmulatorTree.root"
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string(options.outputFile)
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG')
#process.GlobalTag.globaltag = cms.string('POSTLS162_V2::All')
process.GlobalTag.globaltag = cms.string('POSTLS172_V2::All')


process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')

#from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
#process.newRCTConfig = cms.ESSource("PoolDBESSource",
#    CondDBSetup,
#    connect = cms.string('frontier://FrontierPrep/CMS_COND_L1T'),
#    DumpStat=cms.untracked.bool(True),
#    toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string('L1RCTParametersRcd'),
#	    tag = cms.string('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactorsV1')
#        )
#    )
#)
#process.prefer("newRCTConfig")
#process.load("L1Trigger.L1TCalorimeter.caloStage1RCTLuts_cff")
process.load("L1Trigger.L1TCalorimeter.caloStage1RCTLuts_NewRCTCalib_cff")
#process.RCTConfigProducers.hActivityCut = options.hActivityCut
#process.RCTConfigProducers.eActivityCut = options.eActivityCut

process.p1 = cms.Path(
	process.L1TCaloStage1_PPFromRaw
	+process.l1ExtraLayer2
)

process.TauEmul = cms.EDAnalyzer(
	"RateTreeEmul",
	#src = cms.VInputTag(cms.InputTag("L1TCaloRCTToUpgradeConverter","tauJets"))
       # src = cms.VInputTag(cms.InputTag("CaloStage1FinalDigis", "rlxTaus"))
        src = cms.VInputTag(cms.InputTag("l1ExtraLayer2", "Tau"))
)
process.TauEmulIso = cms.EDAnalyzer(
        "RateTreeEmul",
        #src = cms.VInputTag(cms.InputTag("L1TCaloRCTToUpgradeConverter","tauJets"))
       # src = cms.VInputTag(cms.InputTag("CaloStage1FinalDigis", "rlxTaus"))
        src = cms.VInputTag(cms.InputTag("l1ExtraLayer2", "IsoTau"))
)



process.l1RCTParametersTest = cms.EDAnalyzer("L1RCTParametersTester")


process.p2 = cms.Path(
	process.TauEmul
	+process.TauEmulIso
	#+process.l1RCTParametersTest
)

process.schedule = cms.Schedule(
	process.p1,
	process.p2
)
