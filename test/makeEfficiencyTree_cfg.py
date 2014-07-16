#flake8: noqa
'''

Generate trees for measuring and comparing L1 and UCT efficiencies with
respect to RECO objects.

Usage:

    ./makeEfficiencyTree_cfg.py

Optional arguments:

    inputFiles=myFile.root outputFile=outputFile.root maxEvents=-1

Authors: L. Dodd, N. Woods, I. Ojalvo, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

# Get command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# Set useful defaults
#options.inputFiles = '/store/user/tapas/ETauSkim/skim_12_1_erV.root'
#options.inputFiles = '/store/user/ldodd/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6_Fall13dr-tsg_PU40bx25_POSTLS162_V2-v1/05a81b8d696d27a5c3c2ca036967addd/skim_20_1_63z.root'
#options.inputFiles = '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/CC8426A6-B27C-E311-8767-002618943934.root'
#options.inputFiles = '/store/user/ldodd/TT_Tune4C_13TeV-pythia8-tauola/TT_Tune4C_13TeV-pythia8-tauola-tsg_PU40bx25_POSTLS162_V2-v1/fb508503c16d6e4b02bc25104d11f7c2/skim_105_1_Equ.root'
#options.inputFiles = '/store/user/ldodd/TT_Tune4C_13TeV-pythia8-tauola/TT_Tune4C_13TeV-pythia8-tauola-tsg_PU40bx25_POSTLS162_V2-v1/fb508503c16d6e4b02bc25104d11f7c2/skim_104_1_uT6.root'
#options.outputFile = "uct_efficiency_tree.root"
options.register(
    'eicIsolationThreshold',
    3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "EIC Isolation threshold")
#options.register(
#    'hActivityCut',
#    5.0,
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.float,
#    "HCAL activity threshold")
options.register(
    'ecalCalib',
    'CALIB_V4',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'Can be CALIB_V1, CALIB_V3, or CALIB_V4')
options.register(
    'eicCardHcalOnly',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'If 1, turn off the ECAL for the stage1 EGTau path.')
options.register(
    'isMC',
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.')

options.parseArguments()

process = cms.Process("L1UCTEfficiency")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if 'CMSSW_6' in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = 'POSTLS162_V2::All'
    print "Using global tag for upgrade MC: %s" % process.GlobalTag.globaltag
    if not options.isMC:
        raise ValueError("There is no data in CMSSW 6, you must mean isMC=1")
else:
    if not options.isMC:
        # CMSSW 5 data
        process.GlobalTag.globaltag = 'GR_R_53_V21::All'
    else:
        # CMSSW 5 MC
        process.GlobalTag.globaltag = 'START53_V7B::All'
    process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    print "Using global tag for 52X data: %s" % process.GlobalTag.globaltag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    #input = cms.untracked.int32(500)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)

# Load emulation and RECO sequences
if not options.isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")

process.load("L1Trigger.UWTriggerTools.recoObjects_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")

if 'CMSSW_6' in os.environ['CMSSW_VERSION']:
   process.load("L1Trigger.UWTriggerTools.recoObjects_cfi")
else:
   process.load("L1Trigger.UWTriggerTools.recoObjects53X_cfi")


# Common branches to add to the ntuple
common_ntuple_branches = cms.PSet(
    index = cms.string("index"), # Index of reco object in the event
    nRecoObjects = cms.string("nTotalObjects"), # Number of reco objects in the event
    nPVs = cms.string("nPVs"), # number of reco'ed vertices in the event

    # Run, lumi, event number
    run = cms.string("id.run"),
    lumi = cms.string("id.luminosityBlock"),
    evt = cms.string("id.event"),

    recoPt = cms.string("reco.pt"),
    recoEta = cms.string("reco.eta"),
    recoPhi = cms.string("reco.phi"),

    # Whether there exists a L1/UCT object matched to reco
    l1Match = cms.string("l1Match"),
    l1gMatch = cms.string("l1gMatch"),

    l1Pt = cms.string("? l1Match ? l1.pt : 0"),
    l1Eta = cms.string("? l1Match ? l1.eta : 0"),
    l1Phi = cms.string("? l1Match ? l1.phi : 0"),
    l1Type = cms.string("? l1Match ? l1.type() : -1"),


    # TODO add L1extra eta/phi indices

    l1DPhi = cms.string("? l1Match ? deltaPhi(l1.phi, reco.phi) : -1"),
    l1DR = cms.string("? l1Match ? deltaR(l1.eta, l1.phi, reco.eta, reco.phi) : -1"),

    l1gPt = cms.string("? l1gMatch ? l1g.pt : 0"),
    l1gEta = cms.string("? l1gMatch ? l1g.eta : 0"),
    l1gPhi = cms.string("? l1gMatch ? l1g.phi : 0"),

    # For tuning isolation and PU subtraction
    l1gPUM0 = cms.string("? l1gMatch ? l1g.getFloat('puLevelPUM0', -4) : -2"),
    l1gPU = cms.string("? l1gMatch ? l1g.getFloat('puLevel', -4) : -2"),
    l1gPUUIC = cms.string("? l1gMatch ? l1g.getFloat('puLevelUIC', -4) : -2"),
    l1gRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedRegionEt', -4) : -2"),

    l1gEtaCode = cms.vstring("? l1gMatch ? l1g.getInt('rgnEta') : 0", "I"),
    l1gPhiCode = cms.vstring("? l1gMatch ? l1g.getInt('rgnPhi') : 0", "I"),

    l1gDPhi = cms.string("? l1gMatch ? deltaPhi(l1g.phi, reco.phi) : -1"),
    l1gDEta = cms.string("? l1gMatch ? l1g.eta - reco.eta : -10"),
    l1gDR = cms.string("? l1gMatch ? deltaR(l1g.eta, l1g.phi, reco.eta, reco.phi) : -1"),

)

jet_branches = cms.PSet(
    l1gNWRegion = cms.string("? l1gMatch ? l1g.getFloat('neighborNW_et', -4) : -2"),
    l1gNERegion = cms.string("? l1gMatch ? l1g.getFloat('neighborNE_et', -4) : -2"),
    l1gSERegion = cms.string("? l1gMatch ? l1g.getFloat('neighborSE_et', -4) : -2"),
    l1gSWRegion = cms.string("? l1gMatch ? l1g.getFloat('neighborSW_et', -4) : -2"),
    l1gWRegion = cms.string("? l1gMatch ? l1g.getFloat('neighborW_et', -4) : -2"),
    l1gERegion = cms.string("? l1gMatch ? l1g.getFloat('neighborE_et', -4) : -2"),
    l1gNRegion = cms.string("? l1gMatch ? l1g.getFloat('neighborN_et', -4) : -2"),
    l1gSRegion = cms.string("? l1gMatch ? l1g.getFloat('neighborS_et', -4) : -2"),
    l1gJetSeed = cms.string("? l1gMatch ? l1g.getFloat('jetseed_et', -4) : -2"),
)


# Specific to EG tau objects
egtau_branches = cms.PSet(
    l1gSecondRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedSecondRegionEt', -4) : -2"),
    l1gJetPt = cms.string("? l1gMatch ? l1g.getFloat('associatedJetPt', -4) : -2"),
    l1gEllIso = cms.string("? l1gMatch ? l1g.getInt('ellIsolation', -4) : -2"),
    l1gTauVeto = cms.string("? l1gMatch ? l1g.getInt('tauVeto', -4) : -2"),
    l1gTauVetoNeighbor = cms.string("? l1gMatch ? l1g.getInt('associated4x4Tau', -4) : -2"),
    l1gMIP = cms.string("? l1gMatch ? l1g.getInt('mipBit', -4) : -2"),
    l1gIsEle = cms.string("? l1gMatch ? l1g.getInt('isEle', -4) : -2"),
)


# Keep track of electron isolation values
electron_branches = cms.PSet(
    dr03TkSumPt  = cms.string("reco.dr03TkSumPt"),
    dr03EcalRecHitSumEt  = cms.string("reco.dr03EcalRecHitSumEt"),
    dr03HcalTowerSumEt  = cms.string("reco.dr03HcalTowerSumEt"),
    dr03CombinedEt  = cms.string("reco.dr03TkSumPt + reco.dr03EcalRecHitSumEt + reco.dr03HcalTowerSumEt"),
)

# Keep track of information about the ECAL/HCAL composition of taus
tau_branches = cms.PSet(
    emFraction = cms.string("reco.emFraction"),
    decayMode = cms.string("reco.decayMode"),
#    recoEcal = cms.string("ecalEnergy"),
#    recoHcal = cms.string("hcalEnergy"),
    # EK - as far as I can tell, this does not use the lead track at all
#    hcal = cms.string("reco.hcalTotOverPLead"),
)

process.isoTauEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    #recoSrc = cms.VInputTag("isoTaus"),
    recoSrc = cms.VInputTag("recoTaus"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Tau")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)

process.isoTauCheckEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    #recoSrc = cms.VInputTag("isoTaus"),
    recoSrc = cms.VInputTag("isoTaus"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Tau")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)

# Define the tree producers
process.rlxTauEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoTaus"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Tau")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)


# Define the tree producers
process.rlxTauPlusJetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoTaus"),
    l1Src = cms.VInputTag(
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Central"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        tau_branches,
    )
)

# Note that the input electron collection is not isolated, this needs to be done
# at the ntuple level.
process.isoEGEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoElecs"),
    l1Src = cms.VInputTag(cms.InputTag("l1extraParticles", "Isolated")),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        electron_branches,
    )
)

process.rlxEGEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoElecs"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "NonIsolated"),
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        electron_branches,
    )
)

# So we can compare relaxed UCT + ntuple isolation cuts versus stock L1 IsoEG
process.rlxUCTisoL1EGEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoElecs"),
    l1Src = cms.VInputTag(
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
        egtau_branches,
        electron_branches,
    )
)

# Package all of the lepton efficiencies into one sequence
process.leptonEfficiencies = cms.Sequence(

    process.isoTauEfficiency *
    process.isoTauCheckEfficiency *
    process.rlxTauEfficiency
    *process.rlxTauPlusJetEfficiency *

    process.isoEGEfficiency *
    process.rlxEGEfficiency *
    process.rlxUCTisoL1EGEfficiency
)

process.jetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "JetUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches, jet_branches
    )
)

process.corrjetEfficiency = cms.EDAnalyzer(
    "EfficiencyTree",
    recoSrc = cms.VInputTag("recoJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "CorrJetUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevelPUM0Unpacked"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
)

process.highPtPF = cms.EDFilter(
    "GenericPFCandidateSelector",
    src = cms.InputTag("particleFlow"),
    cut = cms.string("charge != 0 && pt > 10")
)

process.printPF = cms.EDAnalyzer(
    "CandInfoPrinter",
    src = cms.InputTag("highPtPF"),
    printHeader=cms.bool(True),
    pt = cms.string("pt"),
    eta = cms.string("eta"),
    pdgId = cms.string("pdgId")
)


process.printTaus = cms.EDAnalyzer(
    "CandInfoPrinter",
    src = cms.InputTag("isoTaus"),
    printHeader=cms.bool(True),
    pt = cms.string("pt"),
    ncands = cms.string("signalPFCands.size"),
    eta = cms.string("eta"),
    phi = cms.string("phi"),
    dm = cms.string("decayMode"),
)

process.printTaus2 = cms.EDAnalyzer(
    "CandInfoPrinter",
    src = cms.InputTag("recoTaus"),
    printHeader=cms.bool(True),
    pt = cms.string("pt"),
    ncands = cms.string("signalPFCands.size"),
    eta = cms.string("eta"),
    phi = cms.string("phi"),
    dm = cms.string("decayMode"),
)



process.printTPGs = cms.EDFilter(
    "TPGDebugger",
    ecalSrc = cms.InputTag("ecalDigis:EcalTriggerPrimitives"),
    #hcalSrc = cms.InputTag("hackHCALMIPs"),
    hcalSrc = cms.InputTag("simHcalTriggerPrimitiveDigis"),
    toPrint = cms.VPSet(
        # everything in all events
        cms.PSet(
            run = cms.uint32(0),
            minIEta = cms.int32(-1000),
            maxIEta = cms.int32(1000),
            maxIPhi = cms.int32(1000),
            minIPhi = cms.int32(-1000),
        )
    )
)

reco_object_step = process.recoObjects

# Determine if we want truth matching or not.
if options.isMC:
    reco_object_step = process.recoObjects_truthMatched
    process.rlxTauPlusJetEfficiency.recoSrc = cms.VInputTag("trueTaus")
    process.isoTauEfficiency.recoSrc = cms.VInputTag("trueTaus")
    process.rlxTauEfficiency.recoSrc = cms.VInputTag("trueTaus")
    process.isoTauEfficiency.recoSrc = cms.VInputTag("trueTaus")
    process.isoTauCheckEfficiency.recoSrc = cms.VInputTag("trueTaus")
    process.rlxTauEfficiency.recoSrc = cms.VInputTag("trueTaus")
    process.printTaus.src=cms.InputTag("trueTaus")

process.p1 = cms.Path(
    reco_object_step *
    process.emulationSequence
#    *process.printTaus 
    #process.highPtPF *
    #process.printPF *
    #process.printTPGs *
    #process.dump *
    #process.pionEfficiency *
)

print "Building Stage1 trees"
process.p1 += process.leptonEfficiencies

################################################################################
###  Semileptonic ttbar skim for sums ###########################################
################################################################################

process.oneMuon = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tightMuons"),
    minNumber = cms.uint32(1),
)

process.jetsPt30 = cms.EDFilter(
    "PFJetSelector",
    src = cms.InputTag("ak5PFJetsNOMuons"),
    filter = cms.bool(True),
    cut = cms.string("pt > 30")
)

process.atLeastThreeJets = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("jetsPt30"),
    minNumber = cms.uint32(3),
)

# Computing CaloMET significance
process.load("RecoMET.METProducers.CaloMETSignif_cfi")
process.metsignificance.noHF = True

# Computing RECO level Sum ET and Sum HT
process.load("L1Trigger.UWTriggerTools.PFSumET_cfi")

process.l1SumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    tree2015 =cms.bool(False),
    l1MHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticles", "MET"),
    # Evan said change l1METSigSrc to match recoMETSigSrc
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    #l1METSigSrc = cms.InputTag("metsignificance"),
    # fixme
    l1SHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticles", "MET"),
    recoMHTSrc = cms.InputTag("pfSumET", "mht"),
    recoMETSrc = cms.InputTag("metNoHF"), # calomet
    #recoMETSigSrc = cms.InputTag("metsignificance"), # does not work in 62X - why?
    recoMETSigSrc  = cms.InputTag("metNoHF"),
    recoSHTSrc = cms.InputTag("pfSumET", "sht"),
    recoSETSrc = cms.InputTag("pfSumET", "set"),
    recoPFMETSrc = cms.InputTag("pfMet"), # pfmet
    recoFile  = cms.bool(True)
)
process.uctSumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    tree2015 =cms.bool(True),
    l1MHTSrc = cms.InputTag("UCT2015Producer", "MHTUnpacked"),
    l1METSrc = cms.InputTag("UCT2015Producer", "METUnpacked"),
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    l1SHTSrc = cms.InputTag("UCT2015Producer", "SHTUnpacked"),
    l1SETSrc = cms.InputTag("UCT2015Producer", "SETUnpacked"),
    recoMETSrc = cms.InputTag("metNoHF"), # calomet
    recoPFMETSrc = cms.InputTag("pfMet"), # pfmet
    recoMHTSrc = cms.InputTag("pfSumET", "mht"), # calomet
#    recoMETSigSrc = cms.InputTag("metsignificance"), #  # does not work in 62X - why?
    recoMETSigSrc  = cms.InputTag("metNoHF"),
    recoSHTSrc = cms.InputTag("pfSumET", "sht"),
    recoSETSrc = cms.InputTag("pfSumET", "set"),
    recoFile  = cms.bool(True)

)

# Make a version of UCT without PU corrections.
process.UCT2015ProducerNoPU = process.UCT2015Producer.clone(
    puMultCorrect = False
)
process.uctSumsNoPUEfficiency = process.uctSumsEfficiency.clone(
    l1MHTSrc = cms.InputTag("UCT2015ProducerNoPU", "MHTUnpacked"),
    l1METSrc = cms.InputTag("UCT2015ProducerNoPU", "METUnpacked"),
    l1METSigSrc = cms.InputTag("UCT2015ProducerNoPU", "METSIGUnpacked"),
    l1SHTSrc = cms.InputTag("UCT2015ProducerNoPU", "SHTUnpacked"),
    l1SETSrc = cms.InputTag("UCT2015ProducerNoPU", "SETUnpacked"),
)

process.semileptonicTTBarPath = cms.Path(
    process.cleanJets *
    process.jetEfficiency *
    process.corrjetEfficiency*
    process.oneMuon *
    process.jetsPt30 *
    process.atLeastThreeJets *
    process.pfSumET *
    #process.metsignificance *
    process.l1SumsEfficiency *
    process.uctSumsEfficiency
    # w/o PU corrections
    #process.UCT2015ProducerNoPU *
    #process.uctSumsNoPUEfficiency
)
# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',    
    'keep *_*_*_L1UCTEfficiency',
   )
)

process.out = cms.EndPath(process.output)



process.schedule = cms.Schedule(
    process.p1,
    process.semileptonicTTBarPath 
       # process.out
)

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

#eic = options.eicIsolationThreshold
#print "Setting EIC threshold to %i" % eic
#process.RCTConfigProducers.eicIsolationThreshold = eic
#hActivity = options.hActivityCut
#print "Setting hActivity threshold to %f" % hActivity
#process.RCTConfigProducers.hActivityCut = hActivity


