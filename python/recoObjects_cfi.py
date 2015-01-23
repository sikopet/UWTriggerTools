#flake8: noqa
import FWCore.ParameterSet.Config as cms

'''

Select good RECO-level muons, electrons, taus, and jets for matching.

Authors: Laura Dodd, Nate Woods, Maria Cepeda, Evan Friis (UW Madison)

'''


#########################################################################
##-------- Find Isolated Muons and Electrons

tightMuons = cms.EDFilter(
    "MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string(
        "pt>10 && isGlobalMuon && isTrackerMuon && abs(eta)<2.5"
        " && globalTrack().normalizedChi2<10"
        " && globalTrack().hitPattern().numberOfValidMuonHits>0"
        " && globalTrack().hitPattern().numberOfValidPixelHits>0"
        " && numberOfMatchedStations>1"
        " && globalTrack().hitPattern().trackerLayersWithMeasurement>5"
        " && (pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.2"
    ),
    filter = cms.bool(False),
)

recoElecs = cms.EDFilter(
    "GsfElectronSelector",
    src = cms.InputTag( "gsfElectrons" ),
    cut = cms.string(
        " (et>20.0)"
        " && (!(-0.02<convDist<0.02 && -0.02<convDcot<0.02))"
        " && (  (isEB"
        " && abs(sigmaIetaIeta)<0.01 &&  abs(deltaPhiSuperClusterTrackAtVtx)<0.03 && abs(deltaEtaSuperClusterTrackAtVtx)<0.004 )"
        " || (isEE"
        " && abs(sigmaIetaIeta)<0.03 &&   abs(deltaPhiSuperClusterTrackAtVtx)<0.02 &&abs(deltaEtaSuperClusterTrackAtVtx)<0.005 )"
        ")"
    ),
    filter = cms.bool(False),
)

# A special, isolated collection of electrons for cleaning the jets
isoElecs = cms.EDFilter(
    "GsfElectronSelector",
    src = cms.InputTag('recoElecs'),
    cut = cms.string(
        "(dr03TkSumPt +  dr03EcalRecHitSumEt+  dr03HcalTowerSumEt)/pt  < 0.3 "
    ),
    filter = cms.bool(False),
)

# Apply jet energy corrections - disabled until we can fix the global tag issue.
#calibratedAK5PFJets = cms.EDProducer(
#    'PFJetCorrectionProducer',
#    src = cms.InputTag('ak5PFJets'),
#    correctors = cms.vstring("ak5PFL1FastL2L3Residual")
#)

##########################################################################
##-------- Remove electrons and muons from jet collection ----------------------
ak5PFJetsNOMuons = cms.EDProducer(
    "PFJetCleaner",
    #srcJets = cms.InputTag("calibratedAK5PFJets"),
    srcJets = cms.InputTag("ak5PFJets"),
    module_label = cms.string(""),
    idLevel = cms.int32(0), # No Jet Id Required
    etaMin  =  cms.double(0.0),
    etaMax  =  cms.double(100),
    ptMin   =  cms.double(0.0),
    srcObjects = cms.VInputTag(cms.InputTag("tightMuons")),
    deltaRMin = cms.double(0.5)
)

recoJets = cms.EDProducer(
    "PFJetCleaner",
    srcJets = cms.InputTag("ak5PFJetsNOMuons"),
    module_label = cms.string(""),
    idLevel = cms.int32(0),  # No Jet Id required
    etaMin  =  cms.double(0.0),
    etaMax  =  cms.double(100),
    ptMin   =  cms.double(0.0),
    srcObjects = cms.VInputTag(cms.InputTag("isoElecs")),
    deltaRMin = cms.double(0.5)
)

cleanJets=cms.Sequence(
    tightMuons*
    recoElecs*
    isoElecs*
    #calibratedAK5PFJets*
    ak5PFJetsNOMuons*
    recoJets
)

# Rerun the PFTau sequence
from Configuration.StandardSequences.GeometryIdeal_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from RecoTauTag.Configuration.RecoPFTauTag_cff import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import * 

# Select good taus
dmTaus = cms.EDFilter(
    "PFTauSelector",
    src = cms.InputTag("hpsPFTauProducer"),
    cut = cms.string("pt > 10 & abs(eta) < 2.5"),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
            selectionCut=cms.double(0.5)
        ),
    ),
    filter = cms.bool(False)
)

isoTaus = cms.EDFilter(
    "PFTauSelector",
    src = cms.InputTag("hpsPFTauProducer"),
    cut = cms.string("pt > 10 & abs(eta) < 2.5"),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
            selectionCut=cms.double(0.5)
        ),
        cms.PSet(
            #discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr"),
            discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
            selectionCut=cms.double(0.5)
        ),
    ),
    filter = cms.bool(False)
)

recoTaus = cms.EDFilter(
    "PFTauSelector",
    src = cms.InputTag("hpsPFTauProducer"),
    cut = cms.string("pt > 10 & abs(eta) < 2.5"),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
            selectionCut=cms.double(0.5)
        ),
        cms.PSet(
            discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
            selectionCut=cms.double(0.5)
        ),
        cms.PSet(
            discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3"),
            selectionCut=cms.double(0.5)
        ),
        cms.PSet(
            discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
            selectionCut=cms.double(0.5)
        ),
    ),
    filter = cms.bool(False)
)


recoObjects = cms.Sequence(
    cleanJets *
#    PFTau * 
    dmTaus *
    isoTaus*
    recoTaus
)

# For matching taus to MC truth

from RecoTauTag.TauTagTools.TauTruthProduction_cfi import tauGenJets, trueHadronicTaus, genParticles

from RecoTauTag.TauTagTools.RecoTauTruthMatching_cfi import recoTauTruthMatcher
recoTauTruthMatcher.src = "recoTaus"
recoTauTruthMatcher.matched = "trueHadronicTaus"

trueTaus = cms.EDFilter(
    "CandViewGenJetMatchRefSelector",
    filter = cms.bool(False),
    src = cms.InputTag("recoTaus"),
    matching = cms.InputTag("recoTauTruthMatcher")
)

#recoObjects_truthMatched = cms.Sequence(
#    recoObjects *
#    genParticles *
#    tauGenJets *
#    trueHadronicTaus *
#    recoTauTruthMatcher *
#    trueTaus
#)
