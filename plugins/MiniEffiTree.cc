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
  VInputTag genSrc_;
  VInputTag srcJet_;
  VInputTag genSrcJet_;
  VInputTag srcEg_;
  VInputTag genSrcEg_;
  VInputTag metSrc_;
  VInputTag genMetSrc_;

  TTree* tree;
  std::vector<Float_t>* pts_;
  std::vector<Float_t>* etas_;
  std::vector<Float_t>* phis_;
  std::vector<Int_t>* genmatches_;

  std::vector<Float_t>* genpts_;
  std::vector<Float_t>* genetas_;
  std::vector<Float_t>* genphis_;
  std::vector<Float_t>* gencharges_;
  std::vector<Float_t>* gendRToJet_;
  std::vector<Float_t>* genMatchedJet_;

  std::vector<Float_t>* jetpts_;
  std::vector<Float_t>* jetetas_;
  std::vector<Float_t>* jetphis_;
  std::vector<Int_t>* genJetmatches_;

  std::vector<Float_t>* genJetpts_;
  std::vector<Float_t>* genJetetas_;
  std::vector<Float_t>* genJetphis_;

  std::vector<Float_t>* egpts_;
  std::vector<Float_t>* egetas_;
  std::vector<Float_t>* egphis_;
  std::vector<Int_t>* egTaumatches_;
  std::vector<Float_t>* dREgTau_;

  std::vector<Float_t>* genEgpts_;
  std::vector<Float_t>* genEgetas_;
  std::vector<Float_t>* genEgphis_;
  std::vector<Int_t>* genEgTaumatches_;
  std::vector<Float_t>* dRgenEgTau_;

  Float_t genMet_Et_;
  Float_t genMet_Phi_;
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
 genmatches_ = new std::vector<Int_t>();

 genpts_ = new std::vector<Float_t>();
 genetas_ = new std::vector<Float_t>();
 genphis_ = new std::vector<Float_t>();
 gencharges_ = new std::vector<Float_t>();
 gendRToJet_ = new std::vector<Float_t>();
 genMatchedJet_  = new std::vector<Float_t>();


 jetpts_ = new std::vector<Float_t>();
 jetetas_ = new std::vector<Float_t>();
 jetphis_ = new std::vector<Float_t>();
 genJetmatches_ = new std::vector<Int_t>();

 genJetpts_ = new std::vector<Float_t>();
 genJetetas_ = new std::vector<Float_t>();
 genJetphis_ = new std::vector<Float_t>();

 egpts_ = new std::vector<Float_t>();
 egetas_ = new std::vector<Float_t>();
 egphis_ = new std::vector<Float_t>();
 egTaumatches_ = new std::vector<Int_t>();
 dREgTau_ = new std::vector<Float_t>();

 genEgpts_ = new std::vector<Float_t>();
 genEgetas_ = new std::vector<Float_t>();
 genEgphis_ = new std::vector<Float_t>();
 genEgTaumatches_ = new std::vector<Int_t>();
 dRgenEgTau_ = new std::vector<Float_t>();

 tree->Branch("pt", "std::vector<float>", &pts_);
 tree->Branch("eta", "std::vector<float>", &etas_);
 tree->Branch("phi", "std::vector<float>", &phis_);
 tree->Branch("genmatches", "std::vector<int>", &genmatches_);

 tree->Branch("genpt", "std::vector<float>", &genpts_);
 tree->Branch("geneta", "std::vector<float>", &genetas_);
 tree->Branch("genphi", "std::vector<float>", &genphis_);
 tree->Branch("gencharge", "std::vector<float>", &gencharges_);
 tree->Branch("gendRToJet", "std::vector<float>", &gendRToJet_);
 tree->Branch("genMatchedJet", "std::vector<float>", &genMatchedJet_);

 tree->Branch("jetpt", "std::vector<float>", &jetpts_);
 tree->Branch("jeteta", "std::vector<float>", &jetetas_);
 tree->Branch("jetphi", "std::vector<float>", &jetphis_);
 tree->Branch("genJetmatches", "std::vector<int>", &genJetmatches_);

 tree->Branch("genJetpt", "std::vector<float>", &genJetpts_);
 tree->Branch("genJeteta", "std::vector<float>", &genJetetas_);
 tree->Branch("genJetphi", "std::vector<float>", &genJetphis_);

 tree->Branch("egpt", "std::vector<float>", &egpts_);
 tree->Branch("egeta", "std::vector<float>", &egetas_);
 tree->Branch("egphi", "std::vector<float>", &egphis_);
 tree->Branch("egTaumatches", "std::vector<int>", &egTaumatches_);
 tree->Branch("dREgTau", "std::vector<float>", &dREgTau_);

 tree->Branch("genEgpt", "std::vector<float>", &genEgpts_);
 tree->Branch("genEgeta", "std::vector<float>", &genEgetas_);
 tree->Branch("genEgphi", "std::vector<float>", &genEgphis_);
 tree->Branch("genEgTaumatches", "std::vector<int>", &genEgTaumatches_);
 tree->Branch("dRgenEgTau" , "std::vector<float>", &dRgenEgTau_);

 tree->Branch("genMet_Et", &genMet_Et_, "genMet_Et/F");
 tree->Branch("genMet_Phi", &genMet_Phi_, "genMet_Phi/F");
 tree->Branch("met_Et", &met_Et_, "met_Et/F");
 tree->Branch("met_Phi", &met_Phi_, "met_Phi/F");

 tree->Branch("run", &run_, "run/i");
 tree->Branch("lumi", &lumi_, "lumi/i");
 tree->Branch("evt", &event_, "evt/l");
 tree->Branch("nPU", &nPU_, "nPU/i");

 src_ = pset.getParameter<VInputTag>("src");
 genSrc_ = pset.getParameter<VInputTag>("genSrc");
 srcJet_ = pset.getParameter<VInputTag>("srcJet");
 genSrcJet_ = pset.getParameter<VInputTag>("genSrcJet");
 srcEg_ = pset.getParameter<VInputTag>("srcEg");
 genSrcEg_ = pset.getParameter<VInputTag>("genSrcEg");
 metSrc_ = pset.getParameter<VInputTag>("metSrc");
 genMetSrc_ = pset.getParameter<VInputTag>("genMetSrc");

 count_=0;
}

MiniEffiTree::~MiniEffiTree() {
 delete pts_;
 delete etas_;
 delete phis_;
 delete genmatches_;
 delete genpts_;
 delete genphis_;
 delete genetas_;
 delete gencharges_;
 delete jetpts_;
 delete jetetas_;
 delete jetphis_;
 delete genJetmatches_;
 delete genJetpts_;
 delete genJetphis_;
 delete genJetetas_;

 delete egpts_;
 delete egetas_;
 delete egphis_;
 delete egTaumatches_;
 delete dREgTau_;
 delete genEgpts_;
 delete genEgphis_;
 delete genEgetas_;
 delete genEgTaumatches_;
 delete dRgenEgTau_;

 delete genMatchedJet_;
 delete gendRToJet_;
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
 std::vector<const reco::Candidate*> genObjects = getCollections(
   evt, genSrc_);
 std::vector<const reco::Candidate*> jetObjects = getCollections(
   evt, srcJet_);
 std::vector<const reco::Candidate*> genJetObjects = getCollections(
   evt, genSrcJet_);
 std::vector<const reco::Candidate*> egObjects = getCollections(
   evt, srcEg_);
 std::vector<const reco::Candidate*> genEgObjects = getCollections(
   evt, genSrcEg_);
 std::vector<const reco::Candidate*> metObject = getCollections(
   evt, metSrc_);
 std::vector<const reco::Candidate*> genMetObject = getCollections(
   evt, genMetSrc_);

 std::sort(objects.begin(), objects.end(), CandPtSorter());
 std::sort(genObjects.begin(), genObjects.end(), CandPtSorter());




 // Clear previous event's objects
 pts_->clear();
 etas_->clear();
 phis_->clear();
 genmatches_->clear();
 genpts_->clear();
 genetas_->clear();
 genphis_->clear();
 gencharges_->clear();
 jetpts_->clear();
 jetetas_->clear();
 jetphis_->clear();
 genJetmatches_->clear();
 genJetpts_->clear();
 genJetetas_->clear();
 genJetphis_->clear();
 egpts_->clear();
 egetas_->clear();
 egphis_->clear();
 egTaumatches_->clear();
 genEgpts_->clear();
 genEgphis_->clear();
 genEgetas_->clear();
 genEgTaumatches_->clear();
 gendRToJet_->clear();
 genMatchedJet_->clear();

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

// std::cout<<"Gen Objects!"<<std::endl;
 for (size_t i = 0; i < genObjects.size(); ++i) {
//  std::cout<<genObjects[i]->pt()<<"   "<<genObjects[i]->eta()<<"   "<<genObjects[i]->phi()<<std::endl;
  genpts_->push_back(genObjects[i]->pt());
  genetas_->push_back(genObjects[i]->eta());
  genphis_->push_back(genObjects[i]->phi());
  gencharges_->push_back(genObjects[i]->charge());
 }


 for (size_t i = 0; i < egObjects.size(); ++i) {
//  std::cout<<objects[i]->pt()<<"   "<<objects[i]->eta()<<"   "<<objects[i]->phi()<<std::endl;
  egpts_->push_back(egObjects[i]->pt());
  egetas_->push_back(egObjects[i]->eta());
  egphis_->push_back(egObjects[i]->phi());
 }

// std::cout<<"Gen Objects!"<<std::endl;
 for (size_t i = 0; i < genEgObjects.size(); ++i) {
//  std::cout<<genObjects[i]->pt()<<"   "<<genObjects[i]->eta()<<"   "<<genObjects[i]->phi()<<std::endl;
  genEgpts_->push_back(genEgObjects[i]->pt());
  genEgetas_->push_back(genEgObjects[i]->eta());
  genEgphis_->push_back(genEgObjects[i]->phi());
 }

 for (size_t i =0; i < objects.size(); ++i){
  int match=-1;
  double minAngle=0.5;
  for (size_t j=0; j<genObjects.size(); ++j){
   //double deltaEta=(etas_->at(j) - genObjects[i]->eta() );
   //double deltaPhi=reco::deltaPhi(phis_->at(j),genObjects[i]->phi()); 
   //double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   double deltaEta=(objects[i]->eta() - genObjects[j]->eta() );
   double deltaPhi=reco::deltaPhi(objects[i]->phi(),genObjects[j]->phi()); 
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi);

   std::cout <<"dEta, dPhi, dR " << deltaEta << " " << deltaPhi << " " << dR << std::endl;
   if(dR<minAngle) {minAngle=dR; match=j;}
  }
  genmatches_->push_back(match);
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



 for (size_t i =0; i < genObjects.size(); ++i){
  int match=-1;
  //double minAngle=0.5;
  double minAngle = 6.0;
  for (size_t j=0; j<genEgObjects.size(); ++j){
   double deltaEta=(genObjects[i]->eta() - genEgObjects[j]->eta() );
   double deltaPhi=reco::deltaPhi(genObjects[i]->phi(),genEgObjects[j]->phi()); 
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi);

   std::cout <<"genEGdEta, genEGdPhi, genEGdR " << deltaEta << " " << deltaPhi << " " << dR << std::endl;
   if(dR<minAngle) {minAngle=dR; match=j;}
  }
  genEgTaumatches_->push_back(match);
  dRgenEgTau_->push_back(minAngle);
 }



//  std::cout<<"-----> matched to"<<match<<std::endl;
/*
  double closestDR=10e6;
  double JET=0;
 for (size_t j = 0; j < genJetObjects.size(); ++j) {
   double deltaEta=(genetas_->at(i) - genJetObjects[j]->eta() );
   double deltaPhi=reco::deltaPhi(genphis_->at(i),genJetObjects[j]->phi());
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   if(dR<0.05) {
                continue;
   }  
   if(dR<closestDR) {closestDR=dR;
                    JET=genJetObjects[j]->pt(); }
 }
 if (closestDR < 1000){
   gendRToJet_->push_back(closestDR);
 }
 else{
   gendRToJet_->push_back(-5);
 }
 genMatchedJet_->push_back(JET);
*/
 


// std::cout<<"jet Objects!"<<std::endl;
 for (size_t i = 0; i < jetObjects.size(); ++i) {
//  std::cout<<jetObjects[i]->pt()<<"   "<<jetObjects[i]->eta()<<"   "<<jetObjects[i]->phi()<<std::endl;
  jetpts_->push_back(jetObjects[i]->pt());
  jetetas_->push_back(jetObjects[i]->eta());
  jetphis_->push_back(jetObjects[i]->phi());
 }

// std::cout<<"Gen Jet Objects!"<<std::endl;
 for (size_t i = 0; i < genJetObjects.size(); ++i) {

  bool isEle=false;
  for (size_t j=0; j<genetas_->size(); j++){
   double deltaEta=(genetas_->at(j) - genJetObjects[i]->eta() );
   double deltaPhi=reco::deltaPhi(genphis_->at(j),genJetObjects[i]->phi());
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   if(dR<0.05) {isEle=true;}
  }
  if (isEle) continue;

  if(genJetObjects[i]->pt()<10) continue;
//  std::cout<<genJetObjects[i]->pt()<<"   "<<genJetObjects[i]->eta()<<"   "<<genJetObjects[i]->phi()<<std::endl;
  genJetpts_->push_back(genJetObjects[i]->pt());
  genJetetas_->push_back(genJetObjects[i]->eta());
  genJetphis_->push_back(genJetObjects[i]->phi());
  int match=-1;
  double minAngle=2;
  for (size_t j=0; j<jetpts_->size(); j++){
   double deltaEta=(jetetas_->at(j) - genJetObjects[i]->eta() );
   double deltaPhi=reco::deltaPhi(jetphis_->at(j),genJetObjects[i]->phi());
   double dR=sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi)   ;
   if(dR<minAngle) {minAngle=dR; match=j;}
  }
  genJetmatches_->push_back(match);
//  std::cout<<"-----> matched to"<<match<<std::endl;
 }

// std::cout<<"MET ?"<<std::endl;
 met_Et_=-10, met_Phi_=-10,genMet_Et_=-10, genMet_Phi_=-10;

 if(metObject.size()==0) {std::cout<<"?????"<<std::endl;}
 else{ met_Et_=metObject[0]->pt();
       met_Phi_=metObject[0]->phi();
  //     std::cout<<"MET: "<<met_Et_<<"   "<<met_Phi_<<std::endl;
 }
 if(genMetObject.size()==0) { std::cout<<"?????"<<std::endl; }
 else{ genMet_Et_=genMetObject[0]->pt();
       genMet_Phi_=genMetObject[0]->phi();
//       std::cout<<"GenMET: "<<genMet_Et_<<"   "<<genMet_Phi_<<std::endl;
 }

 tree->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniEffiTree);
