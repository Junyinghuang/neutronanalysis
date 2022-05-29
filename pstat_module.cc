////////////////////////////////////////////////////////////////////////
// Class:       NeutronCaptureMC
// Plugin Type: analyzer (art v3_01_02)
// File:        NeutronCaptureMC_module.cc
//
// Generated at Tue Mar 26 11:44:15 2019 by Jingbo Wang using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TTree.h"
#include "art_root_io/TFileService.h" 
#include "TH1F.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcore/Geometry/Geometry.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

namespace test {
  class NeutronCaptureMC;
}

class test::NeutronCaptureMC : public art::EDAnalyzer {
public:
  explicit NeutronCaptureMC(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutronCaptureMC(NeutronCaptureMC const&) = delete;
  NeutronCaptureMC(NeutronCaptureMC&&) = delete;
  NeutronCaptureMC& operator=(NeutronCaptureMC const&) = delete;
  NeutronCaptureMC& operator=(NeutronCaptureMC&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const &p);

private:
  //analyzer;
  // Declare member data here.
  TTree *fTree; ///< My TTree
  int fEvent;
    std::string fTruthLabel;

  // Reco variables

  std::vector<float> en_n;
  std::vector<float> en_e;
  std::vector<float> en_mu;
  std::vector<int> pdg;
  std::vector<float> Emuv;
  std::vector<float> nnv;

  double fEkGen;

    
};


test::NeutronCaptureMC::NeutronCaptureMC(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
  // More initializers here.
{
	this->reconfigure(p); 
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::NeutronCaptureMC::beginJob() {
	//create output tree
    //art::ServiceHandle<geo::Geometry> geo_serv;
    
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("mytree", "My Tree");
    fTree->Branch("event", &fEvent, "event/I");
    fTree->Branch("en_e",&en_e);
    fTree->Branch("en_mu",&en_mu);
    fTree->Branch("en_n",&en_n);
    fTree->Branch("pdg",&pdg);
    fTree->Branch("Emuv",&Emuv);
    fTree->Branch("nnv",&nnv);
}


void test::NeutronCaptureMC::endJob() {
    
}

void test::NeutronCaptureMC::reconfigure(fhicl::ParameterSet const & p)
{
 fTruthLabel = p.get<std::string>("TruthLabel");
}

void test::NeutronCaptureMC::analyze(art::Event const& e)
{
    //art::ServiceHandle<geo::Geometry> geo_serv;
  // Implementation of required member function here
    en_n.clear();
    en_e.clear();
    en_mu.clear();
    pdg.clear();
    Emuv.clear();
    nnv.clear();
  fEvent = e.id().event();
    
    double Emu=0;
    double nn=0;
  
  if(!e.isRealData()) {
    // Get the list of MC particles from GEANT
    auto mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    if(mcParticles.isValid()) {
    	
    	for(auto &trueParticle : *mcParticles) {	      
    		fEkGen = (std::sqrt(trueParticle.P()*trueParticle.P() + trueParticle.Mass()*trueParticle.Mass()) - trueParticle.Mass()) * 1000; // MeVs
            pdg.push_back(trueParticle.PdgCode());
            if(trueParticle.Process() == "primary") {
                std::cout<<"PdgCode, Process, Total E(GeV), KineticE(MeV) = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<", "<<fEkGen<<std::endl;
            }
    	  if(trueParticle.Process() == "primary"&&trueParticle.PdgCode() ==2112) {
              en_n.push_back(fEkGen);
    	  }
            if(trueParticle.Process() == "primary"&&trueParticle.PdgCode() ==11) {
                en_e.push_back(fEkGen);
            }
            if(trueParticle.Process() == "primary"&&trueParticle.PdgCode() ==13) {
                en_mu.push_back(fEkGen);
                Emu=fEkGen;
            }
            if(trueParticle.Process() == "primary"&&trueParticle.PdgCode() ==-13) {
                en_mu.push_back(fEkGen);
            }
            if(trueParticle.PdgCode() ==2112) {
                en_n.push_back(fEkGen);
                nn++;
            }
    	}
        Emuv.push_back(Emu);
        nnv.push_back(nn);
        std::cout<<"muon energy: "<<Emu<<" number of neutron: "<<nn<<std::endl;
  }
  }
  // Analysis goes here...
  fTree->Fill();
}



DEFINE_ART_MODULE(test::NeutronCaptureMC)
	

