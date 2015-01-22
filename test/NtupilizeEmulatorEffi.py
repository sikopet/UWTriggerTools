import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.outputFile = "NtupilizeEmulatorEffi.root"
options.parseArguments()
process = cms.Process('EmulatorEffiNtupilizer')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')\

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Spit out filter efficiency at the end.
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#Max Events
process.maxEvents = cms.untracked.PSet(
	#input = cms.untracked.int32(options.maxEvents)
        input = cms.untracked.int32(3)
)

#Input
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles)
#        fileNames = cms.untracked.vstring("/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/C2D9064D-FB7F-E311-BFCC-003048679228.root")
)

#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string(options.outputFile)
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG')
process.GlobalTag.globaltag = cms.string('POSTLS162_V2::All')

process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')
process.load("L1Trigger.UWTriggerTools.recoObjects_cfi")

reco_object_step = process.recoObjects
process.p1 = cms.Path(
        reco_object_step+
	process.L1TCaloStage1_PPFromRaw
	+process.l1ExtraLayer2
)

process.TauEmulEffi = cms.EDAnalyzer(
	"MiniEffiTree",
        src = cms.VInputTag("recoTaus"),
        #src = cms.VInputTag("isoTaus"),
        genSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Tau")),
        #genSrc = cms.VInputTag(cms.InputTag("caloStage1Digis","rlxTaus")),
        srcJet = cms.VInputTag("ak5PFJets"),
        genSrcJet = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Central")),
        srcEg = cms.VInputTag("recoElecs"),
        genSrcEg = cms.VInputTag(cms.InputTag("l1ExtraLayer2","NonIsolated")),
        metSrc = cms.VInputTag("tcMet"),
	genMetSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","MET"))
)

#process.EGEmulEffi = cms.EDAnalyzer(
#	"MiniEffiTree",
#        src = cms.VInputTag("recoElecs"),
#        genSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","NonIsolated")),
#        srcJet = cms.VInputTag("ak5PFJets"),
#        genSrcJet = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Central")),
#        metSrc = cms.VInputTag("tcMet"),
#        genMetSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","MET"))
#)

process.p2 = cms.Path(
	process.TauEmulEffi
 	#+ process.EGEmulEffi
)

process.schedule = cms.Schedule(
	process.p1,
	process.p2
)
