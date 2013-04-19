// -*- C++ -*-
//
// Package:    UCTStage1BProducer
// Class:      UCTStage1BProducer
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Thu Jun  7 13:29:52 CDT 2012
//


// system include files
#include <memory>
#include <math.h>
#include <vector>
#include <list>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "L1Trigger/UCT2015/interface/UCTCandidate.h"
#include "L1Trigger/UCT2015/interface/helpers.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std;
using namespace edm;

class UCTStage1BProducer : public edm::EDProducer {
public:

  // Concrete collection of L1Gobjects (with extra tuning information)
  typedef vector<UCTCandidate> UCTCandidateCollection;
  typedef std::auto_ptr<UCTCandidateCollection> UCTCandidateCollectionPtr;

  explicit UCTStage1BProducer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // Helper methods

  void makeEGs();
  void makeTaus();
  void makeJets();
  void puSubtraction();

  // Adds information about matched regions into a UCTCandidate
  void fillRegionInfo(UCTCandidate& toFill,
      const L1CaloRegionCollection& oldRegions,
      const L1CaloRegionCollection& ecalRegions) const;

  // ----------member data ---------------------------

  bool puCorrect;

  unsigned int egSeed;
  unsigned int tauSeed;

  double regionalHoECut;
  double egRelativeRgnIsolationCut;
  double egRelativeJetIsolationCut;
  double egRelativeEMRgnIsolationCut;
  double egRelativeEMJetIsolationCut;
  double tauRelativeRgnIsolationCut;
  double tauRelativeJetIsolationCut;
  double tauRelativeEMRgnIsolationCut;
  double tauRelativeEMJetIsolationCut;

  double egLSB_;
  double tauLSB_;
  double regionLSB_;

  double puLevelUIC;
  double puLevel;

  Handle<L1CaloRegionCollection> stage1Regions;
  Handle<L1CaloRegionCollection> emRegions;
  Handle<L1CaloEmCollection> tauCands;
  Handle<UCTCandidateCollection> emClusters;

  list<UCTCandidate> rlxEGList;
  list<UCTCandidate> isoEGList;

  list<UCTCandidate> rlxTauList;
  list<UCTCandidate> isoTauList;

  list<UCTCandidate> jetList;

};

void UCTStage1BProducer::fillRegionInfo(
    UCTCandidate& toFill,
    const L1CaloRegionCollection& stage1Regions,
    const L1CaloRegionCollection& emRegions) const {

  unsigned iphi = toFill.getInt("rgnPhi");
  unsigned ieta = toFill.getInt("rgnEta");

  // Region position, relative to iphi, ieta.
  typedef std::pair<int, int> RegionPosition;

  typedef std::map<RegionPosition, UCTRegion> RegionsInfo;

  // keep track of regions at all positions
  RegionsInfo regions;

  // First fill information about the classic regions which have E+H energy,
  // and MIP/tauVeto information.
  for(L1CaloRegionCollection::const_iterator region = stage1Regions.begin();
      region != stage1Regions.end(); region++) {
    int regionEta = region->gctEta();
    int regionPhi = region->gctPhi();
    int deltaEta = regionEta - ieta;
    int deltaPhi = deltaPhiWrapAtN(18, regionPhi, iphi);
    if (std::abs(deltaEta) < 2 && std::abs(deltaPhi) < 2) {
      UCTRegion uctRegion;
      uctRegion.etaPos = deltaEta;
      uctRegion.phiPos = deltaPhi;
      uctRegion.mip = region->mip();
      uctRegion.tauVeto = region->tauVeto();
      uctRegion.et = region->et()*regionLSB_;
      regions[std::make_pair(deltaEta, deltaPhi)] = uctRegion;
    }
  }

  // Now fill the ECAL only information from the CTP regions
  for(L1CaloRegionCollection::const_iterator region = emRegions.begin();
      region != emRegions.end(); region++) {
    int regionEta = region->gctEta();
    int regionPhi = region->gctPhi();
    int deltaEta = regionEta - ieta;
    int deltaPhi = deltaPhiWrapAtN(18, regionPhi, iphi);
    if (std::abs(deltaEta) < 2 && std::abs(deltaPhi) < 2) {

      regions[std::make_pair(deltaEta, deltaPhi)].ecalEt = region->et();
    }
  }

  std::vector<UCTRegion> flatRegions;
  for (RegionsInfo::const_iterator regionInfo = regions.begin();
      regionInfo != regions.end(); ++regionInfo) {
    flatRegions.push_back(regionInfo->second);
  }
  toFill.setRegions(flatRegions);
}

UCTStage1BProducer::UCTStage1BProducer(const edm::ParameterSet& iConfig) :
  puCorrect(iConfig.getParameter<bool>("puCorrect")),
  egSeed(iConfig.getParameter<unsigned int>("egSeed")),
  tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
  regionalHoECut(iConfig.getParameter<double>("regionalHoECut")),
  egRelativeRgnIsolationCut(iConfig.getParameter<double>("egRelativeRgnIsolationCut")),
  egRelativeJetIsolationCut(iConfig.getParameter<double>("egRelativeJetIsolationCut")),
  egRelativeEMRgnIsolationCut(iConfig.getParameter<double>("egRelativeEMRgnIsolationCut")),
  egRelativeEMJetIsolationCut(iConfig.getParameter<double>("egRelativeEMJetIsolationCut")),
  tauRelativeRgnIsolationCut(iConfig.getParameter<double>("tauRelativeRgnIsolationCut")),
  tauRelativeJetIsolationCut(iConfig.getParameter<double>("tauRelativeJetIsolationCut")),
  tauRelativeEMRgnIsolationCut(iConfig.getParameter<double>("tauRelativeEMRgnIsolationCut")),
  tauRelativeEMJetIsolationCut(iConfig.getParameter<double>("tauRelativeEMJetIsolationCut")),
  egLSB_(iConfig.getParameter<double>("egLSB")),
  tauLSB_(iConfig.getParameter<double>("tauLSB")),
  regionLSB_(iConfig.getParameter<double>("regionLSB"))
{
  // Also declare we produce unpacked collections (which have more info)
  produces<UCTCandidateCollection>( "RelaxedEGUnpacked" ) ;
  produces<UCTCandidateCollection>( "IsolatedEGUnpacked" ) ;
  produces<UCTCandidateCollection>( "RelaxedTauUnpacked" ) ;
  produces<UCTCandidateCollection>( "IsolatedTauUnpacked" ) ;
  produces<UCTCandidateCollection>( "JetUnpacked" ) ;
}


// ------------ method called for each event  ------------
void
UCTStage1BProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  iEvent.getByLabel("uctDigis", stage1Regions);
  iEvent.getByLabel("uctDigis", tauCands);
  iEvent.getByLabel("UCT2015EClusterProducer", "EClustersUnpacked", emClusters);
  iEvent.getByLabel("UCT2015EClusterProducer", "ERegions", emRegions);

  makeEGs();
  makeJets();
  makeTaus();

  UCTCandidateCollectionPtr unpackedRlxTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedIsoTaus(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedRlxEGs(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedIsoEGs(new UCTCandidateCollection);
  UCTCandidateCollectionPtr unpackedJets(new UCTCandidateCollection);

  for(list<UCTCandidate>::iterator rlxTau = rlxTauList.begin();
      rlxTau != rlxTauList.end();
      rlxTau++) {
    unpackedRlxTaus->push_back(*rlxTau);
  }
  for(list<UCTCandidate>::iterator isoTau = isoTauList.begin();
      isoTau != isoTauList.end();
      isoTau++) {
    unpackedIsoTaus->push_back(*isoTau);
  }
  for(list<UCTCandidate>::iterator rlxEG = rlxEGList.begin();
      rlxEG != rlxEGList.end();
      rlxEG++) {
    unpackedRlxEGs->push_back(*rlxEG);
  }
  for(list<UCTCandidate>::iterator isoEG = isoEGList.begin();
      isoEG != isoEGList.end();
      isoEG++) {
    unpackedIsoEGs->push_back(*isoEG);
  }
  for(list<UCTCandidate>::iterator jet = jetList.begin();
      jet != jetList.end();
      jet++) {
    unpackedJets->push_back(*jet);
  }

  iEvent.put(unpackedRlxTaus, "RelaxedTauUnpacked");
  iEvent.put(unpackedIsoTaus, "IsolatedTauUnpacked");
  iEvent.put(unpackedRlxEGs, "RelaxedEGUnpacked");
  iEvent.put(unpackedIsoEGs, "IsolatedEGUnpacked");
  iEvent.put(unpackedJets, "JetsUnpacked");

}

void UCTStage1BProducer::puSubtraction()
{
  puLevel = 0;
  puLevelUIC = 0;
  double r_puLevelUIC=0.0;

  int puCount = 0;
  double Rarea=0.0;

  for(L1CaloRegionCollection::const_iterator stage1Region = stage1Regions->begin();
      stage1Region != stage1Regions->end(); stage1Region++){
    double regionEt = stage1Region->et()*regionLSB_;
    double puETMax = 10;
    if(regionEt <= puETMax) {
      puLevel += regionEt;
      puCount++;
      r_puLevelUIC += regionEt;
      Rarea += getRegionArea(stage1Region->gctEta());
    }
  }
  puLevel = puLevel / puCount;
  puLevelUIC = r_puLevelUIC / Rarea;
}

void UCTStage1BProducer::makeEGs() {
  rlxEGList.clear();
  isoEGList.clear();
  for(vector<UCTCandidate>::const_iterator emCluster = emClusters->begin();
      emCluster != emClusters->end();
      emCluster++) {
    if(emCluster->et() > egSeed) {
      double emClusterEt = emCluster->et();
      unsigned emClusterRegionIPhi = emCluster->getInt("rgnPhi");
      unsigned emClusterRegionIEta = emCluster->getInt("rgnEta");

      double C = 0;
      double N = 0;
      double S = 0;
      double E = 0;
      double W = 0;
      double NE = 0;
      double SE = 0;
      double NW = 0;
      double SW = 0;

      int mipBitAtCenter = -1;
      int tauVetoBitAtCenter = -1;

      for(L1CaloRegionCollection::const_iterator region = stage1Regions->begin();
	  region != stage1Regions->end(); region++) {

        // in the physical scale.
        double regionEt = region->et()*regionLSB_;

	if((region->gctPhi() == emClusterRegionIPhi) &&
	   (region->gctEta() == emClusterRegionIEta)) {
	  C = regionEt;
          mipBitAtCenter = region->mip();
          tauVetoBitAtCenter = region->tauVeto();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare emCluster to the total E+H within the region to define rgnIsolation
      // Since the emCluster is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnIsolation
      double rgnIsolation = std::max(0., C - emClusterEt);
      double associatedRegionEt = C;
      double annulusRegions[] = {N, S, E, W, NE, NW, SE, SW};
      double associatedSecondRegionEt = *std::max_element(
          annulusRegions, annulusRegions+8);
      double associatedJetPt = C + N + S + E + W + NE + NW + SE + SW;
      double jetIsolation = associatedJetPt - emClusterEt;

      // Now compute EM only isolations
      C = 0;
      N = 0;
      S = 0;
      E = 0;
      W = 0;
      NE = 0;
      SE = 0;
      NW = 0;
      SW = 0;
      for(L1CaloRegionCollection::const_iterator region = emRegions->begin();
	  region != emRegions->end(); region++) {
        // The region ET is already in the physical scale, from the
        // ECluster producer.
        double regionEt = region->et();
	if((region->gctPhi() == emClusterRegionIPhi) &&
	   (region->gctEta() == emClusterRegionIEta)) {
	  C = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 0) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == 1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), emClusterRegionIPhi) == -1) &&
		(region->gctEta() == (emClusterRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare emCluster to the total E within the region to define rgnEMIsolation
      // Since the emCluster is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnEMIsolation
      double rgnEMIsolation = std::max(0., C - emClusterEt);
      double associatedRegionEtEM = C;
      double annulusRegionsEM[] = {N, S, E, W, NE, NW, SE, SW};
      double associatedSecondRegionEtEM = *std::max_element(
          annulusRegionsEM, annulusRegionsEM+8);
      double associatedJetPtEM = C + N + S + E + W + NE + NW + SE + SW;
      double jetEMIsolation = C + N + S + E + W + NE + NW + SE + SW - emClusterEt;
      // Every emCluster which passes HoECut is an egObject
      double h = associatedRegionEt - emClusterEt;
      double e = emClusterEt;
      if(true || (h / e) < regionalHoECut) {
        // make a copy of the candidate
        UCTCandidate egCand = *emCluster;
	// egObject that passes all isolations is an isolated egObject
	double relativeRgnEMIsolation = (double) rgnEMIsolation / (double) emCluster->pt();
	double relativeJetEMIsolation = (double) jetEMIsolation / (double) emCluster->pt();
	double relativeRgnIsolation = ((double) rgnIsolation) / ((double) emCluster->pt());
	double relativeJetIsolation = ((double) jetIsolation) / ((double) emCluster->pt());

        // Most of the data is already embedded into the ECluster.

        egCand.setFloat("associatedRegionEt", associatedRegionEt);
        egCand.setFloat("associatedSecondRegionEt", associatedSecondRegionEt);
        egCand.setFloat("associatedJetPt", associatedJetPt);

        egCand.setFloat("associatedRegionEtEM", associatedRegionEtEM);
        egCand.setFloat("associatedSecondRegionEtEM", associatedSecondRegionEtEM);
        egCand.setFloat("associatedJetPtEM", associatedJetPtEM);

        egCand.setFloat("effArea", getRegionArea(egCand.getInt("rgnEta")));
        egCand.setFloat("puLevel", puLevel);
        egCand.setFloat("puLevelUIC", puLevelUIC);

        egCand.setFloat("associatedRegionH", h);
        egCand.setFloat("associatedRegionE", e);

        egCand.setInt("mipBit", mipBitAtCenter);
        egCand.setInt("tauVeto", tauVetoBitAtCenter);

	rlxEGList.push_back(egCand);

	if(relativeRgnIsolation < egRelativeRgnIsolationCut &&
	   relativeJetIsolation < egRelativeJetIsolationCut &&
	   relativeRgnEMIsolation < egRelativeEMRgnIsolationCut &&
	   relativeJetEMIsolation < egRelativeEMJetIsolationCut) {
	  isoEGList.push_back(rlxEGList.back());
	}
      }
    }
  }
  rlxEGList.sort();
  isoEGList.sort();
  rlxEGList.reverse();
  isoEGList.reverse();
}

void UCTStage1BProducer::makeJets() {
  jetList.clear();
  for(L1CaloRegionCollection::const_iterator newRegion = stage1Regions->begin();
      newRegion != stage1Regions->end(); newRegion++) {
    double regionET = newRegion->et()*regionLSB_;
    if(regionET > tauSeed) {
      double neighborN_et = 0;
      double neighborS_et = 0;
      double neighborE_et = 0;
      double neighborW_et = 0;
      double neighborNE_et = 0;
      double neighborSW_et = 0;
      double neighborNW_et = 0;
      double neighborSE_et = 0;
      unsigned int nNeighbors = 0;
      for(L1CaloRegionCollection::const_iterator neighbor = stage1Regions->begin();
	  neighbor != stage1Regions->end(); neighbor++) {
	double neighborET = neighbor->et()*regionLSB_;
	if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	   (newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborN_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta()    ) == neighbor->gctEta()) {
	  neighborS_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborE_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborW_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborNE_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborSW_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
		(newRegion->gctEta() - 1) == neighbor->gctEta()) {
	  neighborNW_et = neighborET;
          nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
		(newRegion->gctEta() + 1) == neighbor->gctEta()) {
	  neighborSE_et = neighborET;
          nNeighbors++;
	  continue;
	}
      }
      if(regionET > neighborN_et &&
	 regionET > neighborNW_et &&
	 regionET > neighborW_et &&
	 regionET > neighborSW_et &&
	 regionET >= neighborNE_et &&
	 regionET >= neighborE_et &&
	 regionET >= neighborSE_et &&
	 regionET >= neighborS_et) {
	unsigned int jetET = regionET +
	  neighborN_et + neighborS_et + neighborE_et + neighborW_et +
	  neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;
	// Temporarily use the region granularity -- we will try to improve as above when code is debugged
	int jetPhi = newRegion->gctPhi();
	int jetEta = newRegion->gctEta();

        bool neighborCheck = (nNeighbors == 8);
        // On the eta edge we only expect 5 neighbors
        if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
          neighborCheck = true;

        if (!neighborCheck) {
          std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
          assert(false);
        }
        UCTCandidate theJet(jetET, convertRegionEta(jetEta), convertRegionPhi(jetPhi));
        theJet.setInt("rgnEta", jetEta);
        theJet.setInt("rgnPhi", jetPhi);
        // Embed the puLevel information in the jet object for later tuning
        theJet.setFloat("puLevel", puLevel);
        theJet.setFloat("puLevelUIC", puLevelUIC);
        // Store information about the "core" PT of the jet (central region)
        theJet.setFloat("associatedRegionEt", regionET);
        jetList.push_back(theJet);
      }
    }
  }
  jetList.sort();
  jetList.reverse();
}

void UCTStage1BProducer::makeTaus() {
  rlxTauList.clear();
  isoTauList.clear();
  for(L1CaloEmCollection::const_iterator tauCand = tauCands->begin();
      tauCand != tauCands->end(); tauCand++){
    double tauCandEt = tauCand->rank() * tauLSB_;
    if(tauCandEt > tauSeed) {
      unsigned tauCandRegionIPhi = tauCand->regionId().iphi();
      unsigned tauCandRegionIEta = tauCand->regionId().ieta();
      double C = 0;
      double N = 0;
      double S = 0;
      double E = 0;
      double W = 0;
      double NE = 0;
      double SE = 0;
      double NW = 0;
      double SW = 0;
      int tauVetoBitAtCenter = -1;
      int mipBitAtCenter = -1;
      for(L1CaloRegionCollection::const_iterator region = stage1Regions->begin();
	  region != stage1Regions->end(); region++) {

        double regionEt = region->et()*regionLSB_;


	if((region->gctPhi() == tauCandRegionIPhi) &&
	   (region->gctEta() == tauCandRegionIEta)) {
	  C = regionEt;
          tauVetoBitAtCenter = region->tauVeto();
          mipBitAtCenter = region->mip();
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare tauCand to the total E+H within the region to define rgnIsolation
      // Since the tauCand is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnIsolation
      double associatedJetPt = C + N + S + E + W + NE + NW + SE + SW;
      double associatedRegionEt = C;
      double annulusRegions[] = {N, S, E, W, NE, NW, SE, SW};
      double associatedSecondRegionEt = *std::max_element(
          annulusRegions, annulusRegions+8);

      double rgnIsolation = std::max(0., C - tauCandEt);
      double jetIsolation = C + N + S + E + W + NE + NW + SE + SW - tauCandEt;

      // Now compute EM only isolations
      C = 0;
      N = 0;
      S = 0;
      E = 0;
      W = 0;
      NE = 0;
      SE = 0;
      NW = 0;
      SW = 0;
      for(L1CaloRegionCollection::const_iterator region = emRegions->begin();
	  region != emRegions->end(); region++) {
        // EM regions are already in correct scale from ECluster producer.
        double regionEt = region->et();
	if((region->gctPhi() == tauCandRegionIPhi) &&
	   (region->gctEta() == tauCandRegionIEta)) {
	  C = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  N = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta    ))) {
	  S = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  E = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == 0) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  W = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  NE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == +1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  NW = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta + 1))) {
	  SE = regionEt;
	}
	else if((deltaPhiWrapAtN(18, region->gctPhi(), tauCandRegionIPhi) == -1) &&
		(region->gctEta() == (tauCandRegionIEta - 1))) {
	  SW = regionEt;
	}
      }
      // Compare tauCand to the total E within the region to define rgnEMIsolation
      // Since the tauCand is found including the neighbor EM energy,
      // whereas the regions do not include neighbors, one has to take
      // care to ignore negative rgnEMIsolation
      double associatedJetPtEM = C + N + S + E + W + NE + NW + SE + SW;
      double associatedRegionEtEM = C;
      double annulusRegionsEM[] = {N, S, E, W, NE, NW, SE, SW};
      double associatedSecondRegionEtEM = *std::max_element(
          annulusRegionsEM, annulusRegionsEM+8);

      double rgnEMIsolation = std::max(0., C - tauCandEt);
      unsigned jetEMIsolation = C + N + S + E + W + NE + NW + SE + SW - tauCandEt;
      // Use finer grain position resolution if possible using the emClusters
      // with default value being center of the region
      // phi: 0-17 becomes 0-71; eta: 0-13 becomes 0-55;
      double tauCandEta = convertRegionEta(tauCandRegionIEta);
      double tauCandPhi = convertRegionPhi(tauCandRegionIPhi);

      const UCTCandidate* matchedEmCluster = NULL;

      for(vector<UCTCandidate>::const_iterator emCluster = emClusters->begin();
	  emCluster != emClusters->end();
	  emCluster++) {
	if(std::abs(deltaPhiWrapAtN(18, emCluster->getInt("rgnPhi"), (int)tauCandRegionIPhi)) < 2
            && std::abs(emCluster->getInt("rgnEta") - (int)tauCandRegionIEta) < 2) {
          if (!matchedEmCluster || matchedEmCluster->pt() < emCluster->pt()) {
            tauCandEta = emCluster->eta();
            tauCandPhi = emCluster->phi();
            matchedEmCluster = &(*emCluster);
          }
	}
      }
      UCTCandidate theTau(tauCandEt, tauCandEta, tauCandPhi);
      // tauObject that passes all isolations is an isolated tauObject
      double relativeRgnEMIsolation = ((double) rgnEMIsolation) / ((double) tauCandEt);
      double relativeJetEMIsolation = ((double) jetEMIsolation) / ((double) tauCandEt);
      double relativeRgnIsolation = ((double) rgnIsolation) / ((double) tauCandEt);
      double relativeJetIsolation = ((double) jetIsolation) / ((double) tauCandEt);

      theTau.setFloat("associatedRegionEt", associatedRegionEt);
      theTau.setFloat("associatedSecondRegionEt", associatedSecondRegionEt);
      theTau.setFloat("associatedJetPt", associatedJetPt);

      theTau.setFloat("associatedRegionEtEM", associatedRegionEtEM);
      theTau.setFloat("associatedSecondRegionEtEM", associatedSecondRegionEtEM);
      theTau.setFloat("associatedJetPtEM", associatedJetPtEM);

      theTau.setFloat("emClusterEt", matchedEmCluster ? matchedEmCluster->et() : -1);
      theTau.setFloat("emClusterCenterEt", matchedEmCluster ? matchedEmCluster->getFloat("emClusterCenterEt") : -1);
      theTau.setInt("emClusterCenterFG", matchedEmCluster ? matchedEmCluster->getInt("emClusterCenterFG") : -1);
      theTau.setFloat("emClusterStripEt", matchedEmCluster ? matchedEmCluster->getFloat("emClusterStripEt") : -1);
      theTau.setFloat("puLevel", puLevel);
      theTau.setFloat("puLevelUIC", puLevelUIC);
      theTau.setFloat("puLevelEM", matchedEmCluster ? matchedEmCluster->getFloat("puLevelEM") : -1);
      theTau.setFloat("puLevelUICEM", matchedEmCluster ? matchedEmCluster->getFloat("puLevelUICEM") : -1);
      theTau.setFloat("effArea", getRegionArea(tauCandRegionIEta));
      theTau.setInt("rgnEta", tauCandRegionIEta);
      theTau.setInt("rgnPhi", tauCandRegionIPhi);

      theTau.setInt("mipBit", mipBitAtCenter);
      theTau.setInt("tauVeto", tauVetoBitAtCenter);

      fillRegionInfo(theTau, *stage1Regions, *emRegions);

      rlxTauList.push_back(theTau);

      if(relativeRgnIsolation < tauRelativeRgnIsolationCut &&
	 relativeJetIsolation < tauRelativeJetIsolationCut &&
	 relativeRgnEMIsolation < tauRelativeEMRgnIsolationCut &&
	 relativeJetEMIsolation < tauRelativeEMJetIsolationCut) {
	isoTauList.push_back(rlxTauList.back());
      }
    }
  }
  rlxTauList.sort();
  isoTauList.sort();
  rlxTauList.reverse();
  isoTauList.reverse();
}

//define this as a plug-in
DEFINE_FWK_MODULE(UCTStage1BProducer);
