import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.outputFile = "NtupilizeEmulatorEffiModifyParams.root"
options.register ('eActivityCut',   4, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Tau Veto HCAL Activity Threshold")
options.register ('hActivityCut',   4, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Tau Veto ECAL Activity Threshold")
options.register ('tauThresh',      7, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Tau Seed Threshold")
options.register ('tauNThresh',     0, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Tau Neighbor Seed Threshold")
options.register ('maxPtTauVeto',  64, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                  "Tau max Pt Tau Veto")
options.register ('tauIsoLUTFile',  "/afs/hep.wisc.edu/cms/aglevine/L1Emulator/CMSSW_7_4_0_pre8/src/L1Trigger/L1TCalorimeter/data/tauIsoLutTest0.15Iso.txt", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                  "Tau Isolation Cut")

options.parseArguments()
print '========Tau Parameter Configuration======='
print 'eActivityCut =   ',options.eActivityCut,' GeV'
print 'hActivityCut =   ',options.hActivityCut,' GeV'
print 'tauThresh    =   ',options.tauThresh,' GeV'
print 'tauNThresh   =   ',options.tauNThresh,' GeV'
print 'maxPtTauVeto =  ',options.maxPtTauVeto,' GeV'
print 'tauIsoLUTFile = ',options.tauIsoLUTFile
process = cms.Process('EmulatorEffiNtupilizer')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.Geometry.GeometryIdeal_cff')
# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Spit out filter efficiency at the end.
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#Max Events
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(options.maxEvents)
	#input = cms.untracked.int32(100)
)

#Input
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles)
        #fileNames = cms.untracked.vstring("/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/C2D9064D-FB7F-E311-BFCC-003048679228.root")
)

#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string(options.outputFile)
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG')
process.GlobalTag.globaltag = cms.string('POSTLS172_V2::All')
#########################
from L1Trigger.L1TCalorimeter.caloStage1Params_cfi import *
process.load("L1Trigger.L1TCalorimeter.caloStage1Params_cfi")
process.caloStage1Params.tauSeedThreshold         = cms.double(options.tauThresh)    #pre-RCT Calibration 7GeV
process.caloStage1Params.tauNeighbourThreshold    = cms.double(options.tauNThresh)   #pre-RCT Calibration 0GeV
process.caloStage1Params.tauMaxPtTauVeto          = cms.double(options.maxPtTauVeto) #pre-RCT Calibration 64GeV
process.caloStage1Params.tauIsoLUTFile            = cms.FileInPath(options.tauIsoLUTFile)  #pre-RCT Calibration 0.1

#process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')
# HCAL TP hack
process.load("L1Trigger.L1TCalorimeter.L1TRerunHCALTP_FromRaw_cff")

### Set RCT EG Activity Threshold and Hadronic Activity Threshold Here
#process.load("L1Trigger.L1TCalorimeter.caloStage1RCTLuts_cff")
process.load("L1Trigger.L1TCalorimeter.caloStage1RCTLuts_NewRCTCalib_cff")
process.RCTConfigProducers.hActivityCut = options.hActivityCut
process.RCTConfigProducers.eActivityCut = options.eActivityCut

### RCT To Digi Sequence
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

# RCT
# HCAL input would be from hcalDigis if hack not needed
process.load("L1Trigger.Configuration.SimL1Emulator_cff");
process.simRctDigis.ecalDigis = cms.VInputTag( cms.InputTag( 'ecalDigis:EcalTriggerPrimitives' ) )
process.simRctDigis.hcalDigis = cms.VInputTag( cms.InputTag( 'simHcalTriggerPrimitiveDigis' ) )

### stage 1 
process.load("L1Trigger.L1TCalorimeter.L1TCaloStage1_cff")

### L1Extra
process.load("L1Trigger.Configuration.L1Extra_cff")
process.l1ExtraLayer2 = L1Trigger.Configuration.L1Extra_cff.l1extraParticles.clone()
process.l1ExtraLayer2.isolatedEmSource    = cms.InputTag("simCaloStage1LegacyFormatDigis","isoEm")
process.l1ExtraLayer2.nonIsolatedEmSource = cms.InputTag("simCaloStage1LegacyFormatDigis","nonIsoEm")

process.l1ExtraLayer2.forwardJetSource = cms.InputTag("simCaloStage1LegacyFormatDigis","forJets")
process.l1ExtraLayer2.centralJetSource = cms.InputTag("simCaloStage1LegacyFormatDigis","cenJets")
process.l1ExtraLayer2.tauJetSource     = cms.InputTag("simCaloStage1LegacyFormatDigis","tauJets")
process.l1ExtraLayer2.isoTauJetSource  = cms.InputTag("simCaloStage1LegacyFormatDigis","isoTauJets")

process.l1ExtraLayer2.etTotalSource = cms.InputTag("simCaloStage1LegacyFormatDigis")
process.l1ExtraLayer2.etHadSource   = cms.InputTag("simCaloStage1LegacyFormatDigis")
process.l1ExtraLayer2.etMissSource  = cms.InputTag("simCaloStage1LegacyFormatDigis")
process.l1ExtraLayer2.htMissSource  = cms.InputTag("simCaloStage1LegacyFormatDigis")

process.l1ExtraLayer2.hfRingEtSumsSource    = cms.InputTag("simCaloStage1LegacyFormatDigis")
process.l1ExtraLayer2.hfRingBitCountsSource = cms.InputTag("simCaloStage1LegacyFormatDigis")

## update l1ExtraLayer2 muon tag
process.l1ExtraLayer2.muonSource = cms.InputTag("simGmtDigis")

#########################

# GT
                     # GT
from L1Trigger.Configuration.SimL1Emulator_cff import simGtDigis
process.simGtDigis = simGtDigis.clone()
process.simGtDigis.GmtInputTag = 'simGmtDigis'
process.simGtDigis.GctInputTag = 'caloStage1LegacyFormatDigis'
process.simGtDigis.TechnicalTriggersInputTags = cms.VInputTag( )

process.load("L1Trigger.UWTriggerTools.recoObjects_cfi")

#from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
#process.newRCTConfig = cms.ESSource("PoolDBESSource",
#    CondDBSetup,
#    connect = cms.string('frontier://FrontierPrep/CMS_COND_L1T'),
#    DumpStat=cms.untracked.bool(True),
#    toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string('L1RCTParametersRcd'),
#            tag = cms.string('L1RCTParametersRcd_L1TDevelCollisions_ExtendedScaleFactorsV1')
#        ) 
#    ) 
#)
#process.prefer("newRCTConfig")
reco_object_step = process.recoObjects
process.p1 = cms.Path(
        reco_object_step+
        process.L1TRerunHCALTP_FromRAW
        +process.ecalDigis
        +process.simRctDigis
        +process.L1TCaloStage1
        +process.simGtDigis
        +process.l1ExtraLayer2
)

process.TauEmulEffi = cms.EDAnalyzer(
	"MiniEffiTree",
        src = cms.VInputTag("recoTaus"),
        #src = cms.VInputTag("isoTaus"),
        L1Src = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Tau")),
        #L1Src = cms.VInputTag(cms.InputTag("caloStage1Digis","rlxTaus")),
        srcJet = cms.VInputTag("ak5PFJets"),
        L1SrcJet = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Central")),
        srcEg = cms.VInputTag("recoElecs"),
        L1SrcEg = cms.VInputTag(cms.InputTag("l1ExtraLayer2","NonIsolated")),
        metSrc = cms.VInputTag("tcMet"),
	L1MetSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","MET"))
)
process.TauEmulEffiIso =  cms.EDAnalyzer(
        "MiniEffiTree",
        src = cms.VInputTag("recoTaus"),
        #src = cms.VInputTag("isoTaus"),
        L1Src = cms.VInputTag(cms.InputTag("l1ExtraLayer2","IsoTau")),
        #L1Src = cms.VInputTag(cms.InputTag("caloStage1Digis","rlxTaus")),
        srcJet = cms.VInputTag("ak5PFJets"),
        L1SrcJet = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Central")),
        srcEg = cms.VInputTag("recoElecs"),
        L1SrcEg = cms.VInputTag(cms.InputTag("l1ExtraLayer2","NonIsolated")),
        metSrc = cms.VInputTag("tcMet"),
        L1MetSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","MET"))
)

process.l1RCTParametersTest = cms.EDAnalyzer("L1RCTParametersTester")


#process.EGEmulEffi = cms.EDAnalyzer(
#	"MiniEffiTree",
#        src = cms.VInputTag("recoElecs"),
#        L1Src = cms.VInputTag(cms.InputTag("l1ExtraLayer2","NonIsolated")),
#        srcJet = cms.VInputTag("ak5PFJets"),
#        L1SrcJet = cms.VInputTag(cms.InputTag("l1ExtraLayer2","Central")),
#        metSrc = cms.VInputTag("tcMet"),
#        L1MetSrc = cms.VInputTag(cms.InputTag("l1ExtraLayer2","MET"))
#)

process.p2 = cms.Path(
	process.TauEmulEffi
	+process.TauEmulEffiIso
	#+process.l1RCTParametersTest
 	#+ process.EGEmulEffi
)

process.schedule = cms.Schedule(
	process.p1,
	process.p2
)
