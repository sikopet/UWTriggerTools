#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2TauAlgorithmImp.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class RateTreeEmul : public edm::EDAnalyzer {
  public:
    RateTreeEmul(const edm::ParameterSet& pset);
    virtual ~RateTreeEmul();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    VInputTag src_;
    TTree* tree;
    std::vector<Float_t>* pts_;
    std::vector<Float_t>* etas_;
    std::vector<Float_t>* phis_;
    std::vector<Int_t>* isoFlags_;
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;
    Int_t nPU_;
};

RateTreeEmul::RateTreeEmul(const edm::ParameterSet& pset) {
  // Initialize the ntuple builder
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Ntuple", "Ntuple");
  pts_ = new std::vector<Float_t>();
  etas_ = new std::vector<Float_t>();
  phis_ = new std::vector<Float_t>();
  isoFlags_ = new std::vector<Int_t>();

  tree->Branch("pt", "std::vector<float>", &pts_);
  tree->Branch("eta", "std::vector<float>", &etas_);
  tree->Branch("phi", "std::vector<float>", &phis_);
  tree->Branch("run", &run_, "run/i");
  //tree->Branch("isoFlags", "std::vector<int>", &isoFlags_);
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt", &event_, "evt/l");
  tree->Branch("nPU", &nPU_, "nPU/i");

  src_ = pset.getParameter<VInputTag>("src");
}

RateTreeEmul::~RateTreeEmul() {
  delete pts_;
  delete etas_;
  delete phis_;
  delete isoFlags_;
}


namespace {

  // Predicate to sort candidates by descending pt
  class CandPtSorter {
    public:
      bool operator()(const reco::Candidate* candA, const reco::Candidate* candB)
        const {
          return candA->pt() > candB->pt();
        }
  };

  // Turn a set of InputTags into a colleciton of candidate pointers.
  std::vector<const reco::Candidate*> getCollections(
      const edm::Event& evt, const VInputTag& collections) {
    std::vector<const reco::Candidate*> output;
    // Loop over collections
    for (size_t i = 0; i < collections.size(); ++i) {
      edm::Handle<edm::View<reco::Candidate> > handle;
      evt.getByLabel(collections[i], handle);
      // Loop over objects in current collection
      for (size_t j = 0; j < handle->size(); ++j) {
        const reco::Candidate& object = handle->at(j);
        output.push_back(&object);
      }
    }
    return output;
  }

}

void RateTreeEmul::analyze(const edm::Event& evt, const edm::EventSetup& es) {
  //std::cout << "getting objects" << std::endl;	
  // Get the objects.
  std::vector<const reco::Candidate*> objects = getCollections(
      evt, src_);
  //std::cout << "sorting objects" << std::endl;
  std::sort(objects.begin(), objects.end(), CandPtSorter());

  //std::cout << "clearing objects" << std::endl;
  // Clear previous event's objects
  pts_->clear();
  etas_->clear();
  phis_->clear();
  isoFlags_->clear();
  
  
  //std::cout << "setting up meta info" <<std::endl;
  // Setup meta info
  //  
  run_ = evt.id().run();
  lumi_ = evt.id().luminosityBlock();
  event_ = evt.id().event();
  
  //std::cout << "event: " << event_ << std::endl;
  
  //std::cout << "getting pileup summary" << std::endl;
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  evt.getByLabel("addPileupInfo", PupInfo);

  nPU_=-1;
  for(std::vector<PileupSummaryInfo>::const_iterator i = PupInfo->begin(); i!=PupInfo->end();++i) {
   int BX = i->getBunchCrossing();
   if(BX==0) {
     nPU_ =i->getTrueNumInteractions();
   }
   else {
     nPU_=i->getPU_NumInteractions();
   }
  }


 
  //std::cout << "pushing back" << std::endl;
  for (size_t i = 0; i < objects.size(); ++i) {
    pts_->push_back(objects[i]->pt());
    etas_->push_back(objects[i]->eta());
    phis_->push_back(objects[i]->phi());
    //std::cout <<"about to make dynamic cast"<<std::endl;
    //const l1t::Tau* tau = dynamic_cast <const l1t::Tau*>(objects[i]);
    //std::cout << "about to push back hwIso" <<std::endl;
    //std::cout << tau[i].hwIso() << std::endl;
    //isoFlags_->push_back(tau[i].hwIso());
  }


  tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RateTreeEmul);
