/*
 * =====================================================================================
 *
 *       Filename:  MiniEffiTree.cc
 *        Company:  UW Madison
 *
 * =====================================================================================
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include "TTree.h"

typedef std::vector<edm::InputTag> VInputTag;

class MiniEffiTree : public edm::EDAnalyzer {
 public:
  MiniEffiTree(const edm::ParameterSet& pset);
  virtual ~MiniEffiTree();
  void analyze(const edm::Event& evt, const edm::EventSetup& es);
 private:
  VInputTag src_;
  VInputTag L1Src_;
  VInputTag srcJet_;
  VInputTag L1SrcJet_;
  VInputTag srcEg_;
  VInputTag L1SrcEg_;
  VInputTag metSrc_;
  VInputTag L1MetSrc_;

  TTree* tree;
  std::vector<Float_t>* pts_;
  std::vector<Float_t>* etas_;
  std::vector<Float_t>* phis_;
  std::vector<Int_t>* L1matches_;

  std::vector<Float_t>* L1pts_;
  std::vector<Float_t>* L1etas_;
  std::vector<Float_t>* L1phis_;
  std::vector<Float_t>* L1charges_;

  std::vector<Float_t>* L1Matchedpts_;
  std::vector<Float_t>* L1Matchedetas_;
  std::vector<Float_t>* L1Matchedphis_;
  std::vector<Float_t>* L1Matchedcharges_;

  std::vector<Float_t>* L1dRToJet_;
  std::vector<Float_t>* L1MatchedJet_;

  std::vector<Float_t>* jetpts_;
  std::vector<Float_t>* jetetas_;
  std::vector<Float_t>* jetphis_;
  std::vector<Int_t>* L1Jetmatches_;

  std::vector<Float_t>* L1Jetpts_;
  std::vector<Float_t>* L1Jetetas_;
  std::vector<Float_t>* L1Jetphis_;

  std::vector<Float_t>* egpts_;
  std::vector<Float_t>* egetas_;
  std::vector<Float_t>* egphis_;
  std::vector<Int_t>* egTaumatches_;
  std::vector<Float_t>* dREgTau_;

  std::vector<Float_t>* L1Egpts_;
  std::vector<Float_t>* L1Egetas_;
  std::vector<Float_t>* L1Egphis_;
  std::vector<Int_t>* L1EgTaumatches_;
  std::vector<Float_t>* dRL1EgTau_;

  Float_t L1Met_Et_;
  Float_t L1Met_Phi_;
  Float_t met_Et_;
  Float_t met_Phi_;


  UInt_t run_;
  UInt_t lumi_;
  ULong64_t event_;
  Int_t nPU_;
 
  Int_t count_;
};

MiniEffiTree::MiniEffiTree(const edm::ParameterSet& pset) {
 // Initialize the ntuple builder
 edm::Service<TFileService> fs;
 tree = fs->make<TTree>("Ntuple", "Ntuple");
 pts_ = new std::vector<Float_t>();
 etas_ = new std::vector<Float_t>();
 phis_ = new std::vector<Float_t>();
 L1matches_ = new std::vector<Int_t>();

 L1pts_ = new std::vector<Float_t>();
 L1etas_ = new std::vector<Float_t>();
 L1phis_ = new std::vector<Float_t>();
 L1charges_ = new std::vector<Float_t>();

 L1Matchedpts_ = new std::vector<Float_t>();
 L1Matchedetas_ = new std::vector<Float_t>();
 L1Matchedphis_ = new std::vector<Float_t>();
 L1Matchedcharges_ = new std::vector<Float_t>();

 L1dRToJet_ = new std::vector<Float_t>();
 L1MatchedJet_  = new std::vector<Float_t>();


 jetpts_ = new std::vector<Float_t>();
 jetetas_ = new std::vector<Float_t>();
 jetphis_ = new std::vector<Float_t>();
 L1Jetmatches_ = new std::vector<Int_t>();

 L1Jetpts_ = new std::vector<Float_t>();
 L1Jetetas_ = new std::vector<Float_t>();
 L1Jetphis_ = new std::vector<Float_t>();

 egpts_ = new std::vector<Float_t>();
 egetas_ = new std::vector<Float_t>();
 egphis_ = new std::vector<Float_t>();
 egTaumatches_ = new std::vector<Int_t>();
 dREgTau_ = new std::vector<Float_t>();

 L1Egpts_ = new std::vector<Float_t>();
 L1Egetas_ = new std::vector<Float_t>();
 L1Egphis_ = new std::vector<Float_t>();
 L1EgTaumatches_ = new std::vector<Int_t>();
 dRL1EgTau_ = new std::vector<Float_t>();

 tree->Branch("pt", "std::vector<float>", &pts_);
 tree->Branch("eta", "std::vector<float>", &etas_);
 tree->Branch("phi", "std::vector<float>", &phis_);
 tree->Branch("L1matches", "std::vector<int>", &L1matches_);

 tree->Branch("L1pt", "std::vector<float>", &L1pts_);
 tree->Branch("L1eta", "std::vector<float>", &L1etas_);
 tree->Branch("L1phi", "std::vector<float>", &L1phis_);
 tree->Branch("L1charge", "std::vector<float>", &L1charges_);

 tree->Branch("L1Matchedpt", "std::vector<float>", &L1Matchedpts_);
 tree->Branch("L1Matchedeta", "std::vector<float>", &L1Matchedetas_);
 tree->Branch("L1Matchedphi", "std::vector<float>", &L1Matchedphis_);
 tree->Branch("L1Matchedcharge", "std::vector<float>", &L1Matchedcharges_);

 tree->Branch("L1dRToJet", "std::vector<float>", &L1dRToJet_);
 tree->Branch("L1MatchedJet", "std::vector<float>", &L1MatchedJet_);

 tree->Branch("jetpt", "std::vector<float>", &jetpts_);
 tree->Branch("jeteta", "std::vector<float>", &jetetas_);
 tree->Branch("jetphi", "std::vector<float>", &jetphis_);
 tree->Branch("L1Jetmatches", "std::vector<int>", &L1Jetmatches_);

 tree->Branch("L1Jetpt", "std::vector<float>", &L1Jetpts_);
 tree->Branch("L1Jeteta", "std::vector<float>", &L1Jetetas_);
 tree->Branch("L1Jetphi", "std::vector<float>", &L1Jetphis_);

 tree->Branch("egpt", "std::vector<float>", &egpts_);
 tree->Branch("egeta", "std::vector<float>", &egetas_);
 tree->Branch("egphi", "std::vector<float>", &egphis_);
 tree->Branch("egTaumatches", "std::vector<int>", &egTaumatches_);
 tree->Branch("dREgTau", "std::vector<float>", &dREgTau_);

 tree->Branch("L1Egpt", "std::vector<float>", &L1Egpts_);
 tree->Branch("L1Egeta", "std::vector<float>", &L1Egetas_);
 tree->Branch("L1Egphi", "std::vector<float>", &L1Egphis_);
 tree->Branch("L1EgTaumatches", "std::vector<int>", &L1EgTaumatches_);
 tree->Branch("dRL1EgTau" , "std::vector<float>", &dRL1EgTau_);

 tree->Branch("L1Met_Et", &L1Met_Et_, "L1Met_Et/F");
 tree->Branch("L1Met_Phi", &L1Met_Phi_, "L1Met_Phi/F");
 tree->Branch("met_Et", &met_Et_, "met_Et/F");
 tree->Branch("met_Phi", &met_Phi_, "met_Phi/F");

 tree->Branch("run", &run_, "run/i");
 tree->Branch("lumi", &lumi_, "lumi/i");
 tree->Branch("evt", &event_, "evt/l");
 tree->Branch("nPU", &nPU_, "nPU/i");

 src_ = pset.getParameter<VInputTag>("src");
 L1Src_ = pset.getParameter<VInputTag>("L1Src");
 srcJet_ = pset.getParameter<VInputTag>("srcJet");
 L1SrcJet_ = pset.getParameter<VInputTag>("L1SrcJet");
 srcEg_ = pset.getParameter<VInputTag>("srcEg");
 L1SrcEg_ = pset.getParameter<VInputTag>("L1SrcEg");
 metSrc_ = pset.getParameter<VInputTag>("metSrc");
 L1MetSrc_ = pset.getParameter<VInputTag>("L1MetSrc");

 count_=0;
}

MiniEffiTree::~MiniEffiTree() {
 delete pts_;
 delete etas_;
 delete phis_;
 delete L1matches_;
 delete L1pts_;
 delete L1phis_;
 delete L1etas_;
 delete L1charges_;
 delete L1Matchedpts_;
 delete L1Matchedphis_;
 delete L1Matchedetas_;
 delete L1Matchedcharges_;
 delete jetpts_;
 delete jetetas_;
 delete jetphis_;
 delete L1Jetmatches_;
 delete L1Jetpts_;
 delete L1Jetphis_;
 delete L1Jetetas_;

 delete egpts_;
 delete egetas_;
 delete egphis_;
 delete egTaumatches_;
 delete dREgTau_;
 delete L1Egpts_;
 delete L1Egphis_;
 delete L1Egetas_;
 delete L1EgTaumatches_;
 delete dRL1EgTau_;

 delete L1MatchedJet_;
 delete L1dRToJet_;
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

void MiniEffiTree::analyze(const edm::Event& evt, const edm::EventSetup& es) {

 count_++;

 // Get the objects.
 std::vector<const reco::Candidate*> objects = getCollections(
   evt, src_);
 std::vector<const reco::Candidate*> L1Objects = getCollections(
   evt, L1Src_);
 std::vector<const reco::Candidate*> jetObjects = getCollections(
   evt, srcJet_);
 std::vector<const reco::Candidate*> L1JetObjects = getCollections(
   evt, L1SrcJet_);
 std::vector<const reco::Candidate*> egObjects = getCollections(
   evt, srcEg_);
 std::vector<const reco::Candidate*> L1EgObjects = getCollections(
   evt, L1SrcEg_);
 std::vector<const reco::Candidate*> metObject = getCollections(
   evt, metSrc_);
 std::vector<const reco::Candidate*> L1MetObject = getCollections(
   evt, L1MetSrc_);

 std::sort(objects.begin(), objects.end(), CandPtSorter());
 std::sort(L1Objects.begin(), L1Objects.end(), CandPtSorter());




 // Clear previous event's objects
 pts_->clear();
 etas_->clear();
 phis_->clear();
 L1matches_->clear();
 L1pts_->clear();
 L1etas_->clear();
 L1phis_->clear();
 L1charges_->clear();
 L1Matchedpts_->clear();
 L1Matchedetas_->clear();
 L1Matchedphis_->clear();
 L1Matchedcharges_->clear();
 jetpts_->clear();
 jetetas_->clear();
 jetphis_->clear();
 L1Jetmatches_->clear();
 L1Jetpts_->clear();
 L1Jetetas_->clear();
 L1Jetphis_->clear();
 egpts_->clear();
 egetas_->clear();
 egphis_->clear();
 egTaumatches_->clear();
 L1Egpts_->clear();
 L1Egphis_->clear();
 L1Egetas_->clear();
 L1EgTaumatches_->clear();
 L1dRToJet_->clear();
 L1MatchedJet_->clear();

 // Setup meta info
 run_ = evt.id().run();
 lumi_ = evt.id().luminosityBlock();
 event_ = evt.id().event();

 edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
 evt.getByLabel("addPileupInfo", PupInfo);

 nPU_=-1;
 for(std::vector<PileupSummaryInfo>::const_iterator i = PupInfo->begin();
   i!=PupInfo->end();++i) {
  int BX = i->getBunchCrossing();
  if(BX==0) {
   nPU_ =i->getTrueNumInteractions();
   nPU_=i->getPU_NumInteractions();
  }
 }
 

// std::cout<<"----------------------------"<<count_<<std::endl;

// std::cout<<"Trigger Objects!"<<std::endl;
 for (size_t i = 0; i < objects.size(); ++i) {
//  std::cout<<objects[i]->pt()<<"   "<<objects[i]->eta()<<"   "<<objects[i]->phi()<<std::endl;
  pts_->push_back(objects[i]->pt());
  etas_->push_back(objects[i]->eta());
  phis_->push_back(objects[i]->phi());
 }

 //std::cout<<"L1 Objects!"<<std::endl;
 for (size_t i = 0; i < L1Objects.size(); ++i) {
  //std::cout<<L1Objects[i]->pt()<<"   "<<L1Objects[i]->eta()<<"   "<<L1Objects[i]->phi()<<std::endl;
  L1pts_->push_back(L1Objects[i]->pt());
  L1etas_->push_back(L1Objects[i]->eta());
  L1phis_->push_back(L1Objects[i]->phi());
  L1charges_->push_back(L1Objects[i]->charge());
 }


 for (size_t i = 0; i < egObjects.size(); ++i) {
//  std::cout<<objects[i]->pt()<<"   "<<objects[i]->eta()<<"   "<<objects[i]->phi()<<std::endl;
  egpts_->push_back(egObjects[i]->pt());
  egetas_->push_back(egObjects[i]->eta());
  egphis_->push_back(egObjects[i]->phi());
 }

// std::cout<<"Gen Objects!"<<std::endl;
 for (size_t i = 0; i < L1EgObjects.size(); ++i) {
//  std::cout<<L1Objects[i]->pt()<<"   "<<L1Objects[i]->eta()<<"   "<<L1Objects[i]->phi()<<std::endl;
  L1Egpts_->push_back(L1EgObjects[i]->pt());
  L1Egetas_->push_back(L1EgObjects[i]->eta());
  L1Egphis_->push_back(L1EgObjects[i]->phi());
 }

 for (size_t i =0; i < objects.size(); ++i){
  int matchj=-1;
  double minAngle=0.5;
  for (size_t j=0; j<L1Objects.size(); ++j){
   //double deltaEta=(etas_->at(j) - L1Objects[i]->eta() );
   //double deltaPhi=reco::deltaPhi(phis_->at(j),L1Objects[i]->phi()); 
   //double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   double deltaEta=(objects[i]->eta() - L1Objects[j]->eta() );
   double deltaPhi=reco::deltaPhi(objects[i]->phi(),L1Objects[j]->phi()); 
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi);
   std::cout << "pt, L1pt: " << objects[i]->pt() <<", "<<L1Objects[j]->pt() <<", "<<std::endl;

   std::cout <<"dEta, dPhi, dR: " << deltaEta << ", " << deltaPhi << ", " << dR << std::endl;
   if(dR<minAngle) {minAngle=dR; 
      matchj = j;
      std::cout << "dR < minAngle, matched!" << std::endl;}
  }
  for (size_t j=0; j<L1Objects.size(); ++j){
    if (j == matchj){
      L1matches_->push_back(1);
      std::cout << "pushing back 1 for match pt, L1pt: " << objects[i]->pt() <<", "<<L1Objects[j]->pt() <<std::endl;
      L1Matchedpts_->push_back(L1Objects[j]->pt());
      L1Matchedetas_->push_back(L1Objects[j]->eta());
      L1Matchedphis_->push_back(L1Objects[j]->phi());
      L1Matchedcharges_->push_back(L1Objects[j]->charge());
    }
//    else{
 //     L1matches_->push_back(-1);
  //    std::cout << "pushing back -1 for nonmatch pt, L1pt :" << objects[i]->pt() <<", "<<L1Objects[j]->pt() <<std::endl;
   // }
  }
 }

 for (size_t i =0; i < objects.size(); ++i){
  int match=-1;
  //double minAngle=0.5;
  double minAngle = 6.0;
  std::cout << "starting objects loop" << std::endl;
  std::cout << egObjects.size() << std::endl;
  for (size_t j=0; j<egObjects.size(); ++j){
   std::cout << "starting eg objects loop" << std::endl;
   double deltaEta=(objects[i]->eta() - egObjects[j]->eta() );
   double deltaPhi=reco::deltaPhi(objects[i]->phi(),egObjects[j]->phi()); 
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi);

   std::cout <<"EGdEta, EGdPhi, EGdR " << deltaEta << " " << deltaPhi << " " << dR << std::endl;
   if(dR<minAngle) {minAngle=dR; match=j;}
  }
  egTaumatches_->push_back(match);
  dREgTau_->push_back(minAngle);
  
 }



 for (size_t i =0; i < L1Objects.size(); ++i){
  int match=-1;
  //double minAngle=0.5;
  double minAngle = 6.0;
  for (size_t j=0; j<L1EgObjects.size(); ++j){
   double deltaEta=(L1Objects[i]->eta() - L1EgObjects[j]->eta() );
   double deltaPhi=reco::deltaPhi(L1Objects[i]->phi(),L1EgObjects[j]->phi()); 
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi);

   std::cout <<"L1EGdEta, L1EGdPhi, L1EGdR " << deltaEta << " " << deltaPhi << " " << dR << std::endl;
   if(dR<minAngle) {minAngle=dR; match=j;}
  }
  L1EgTaumatches_->push_back(match);
  dRL1EgTau_->push_back(minAngle);
 }



//  std::cout<<"-----> matched to"<<match<<std::endl;
/*
  double closestDR=10e6;
  double JET=0;
 for (size_t j = 0; j < L1JetObjects.size(); ++j) {
   double deltaEta=(L1etas_->at(i) - L1JetObjects[j]->eta() );
   double deltaPhi=reco::deltaPhi(L1phis_->at(i),L1JetObjects[j]->phi());
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   if(dR<0.05) {
                continue;
   }  
   if(dR<closestDR) {closestDR=dR;
                    JET=L1JetObjects[j]->pt(); }
 }
 if (closestDR < 1000){
   L1dRToJet_->push_back(closestDR);
 }
 else{
   L1dRToJet_->push_back(-5);
 }
 L1MatchedJet_->push_back(JET);
*/
 


// std::cout<<"jet Objects!"<<std::endl;
 for (size_t i = 0; i < jetObjects.size(); ++i) {
//  std::cout<<jetObjects[i]->pt()<<"   "<<jetObjects[i]->eta()<<"   "<<jetObjects[i]->phi()<<std::endl;
  jetpts_->push_back(jetObjects[i]->pt());
  jetetas_->push_back(jetObjects[i]->eta());
  jetphis_->push_back(jetObjects[i]->phi());
 }

// std::cout<<"Gen Jet Objects!"<<std::endl;
 for (size_t i = 0; i < L1JetObjects.size(); ++i) {

  bool isEle=false;
  for (size_t j=0; j<L1etas_->size(); j++){
   double deltaEta=(L1etas_->at(j) - L1JetObjects[i]->eta() );
   double deltaPhi=reco::deltaPhi(L1phis_->at(j),L1JetObjects[i]->phi());
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   if(dR<0.05) {isEle=true;}
  }
  if (isEle) continue;

  if(L1JetObjects[i]->pt()<10) continue;
//  std::cout<<L1JetObjects[i]->pt()<<"   "<<L1JetObjects[i]->eta()<<"   "<<L1JetObjects[i]->phi()<<std::endl;
  L1Jetpts_->push_back(L1JetObjects[i]->pt());
  L1Jetetas_->push_back(L1JetObjects[i]->eta());
  L1Jetphis_->push_back(L1JetObjects[i]->phi());
  int match=-1;
  double minAngle=2;
  for (size_t j=0; j<jetpts_->size(); j++){
   double deltaEta=(jetetas_->at(j) - L1JetObjects[i]->eta() );
   double deltaPhi=reco::deltaPhi(jetphis_->at(j),L1JetObjects[i]->phi());
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   if(dR<minAngle) {minAngle=dR; match=j;}
  }
  L1Jetmatches_->push_back(match);
//  std::cout<<"-----> matched to"<<match<<std::endl;
 }

// std::cout<<"MET ?"<<std::endl;
 met_Et_=-10, met_Phi_=-10,L1Met_Et_=-10, L1Met_Phi_=-10;

 if(metObject.size()==0) {std::cout<<"?????"<<std::endl;}
 else{ met_Et_=metObject[0]->pt();
       met_Phi_=metObject[0]->phi();
  //     std::cout<<"MET: "<<met_Et_<<"   "<<met_Phi_<<std::endl;
 }
 if(L1MetObject.size()==0) { std::cout<<"?????"<<std::endl; }
 else{ L1Met_Et_=L1MetObject[0]->pt();
       L1Met_Phi_=L1MetObject[0]->phi();
//       std::cout<<"GenMET: "<<L1Met_Et_<<"   "<<L1Met_Phi_<<std::endl;
 }

 tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniEffiTree);
