////////////////////////////////////////////////////////////////////////
// Class:       NeutronCaptureFinder
// Plugin Type: analyzer (art v3_01_02)
// File:        NeutronCaptureFinder_module.cc
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
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace test {
  class NeutronCaptureFinder;
}

class test::NeutronCaptureFinder : public art::EDAnalyzer {
public:
  explicit NeutronCaptureFinder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutronCaptureFinder(NeutronCaptureFinder const&) = delete;
  NeutronCaptureFinder(NeutronCaptureFinder&&) = delete;
  NeutronCaptureFinder& operator=(NeutronCaptureFinder const&) = delete;
  NeutronCaptureFinder& operator=(NeutronCaptureFinder&&) = delete;

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
  std::string fParticleLabel;
  std::string fTrackLabel; 
  std::string fHitLabel;
  std::string fSpacePointLabel;
  std::string fSliceLabel;
  art::InputTag fSimEdepTag;
  	
  // Truth variables
  float	fTrueEnergy;
  
  // Reco variables
  int fNPrimaries;
  float fLength;
  TH1F *fLengthHist;
  int fTrackHits;
  int fParticleHits;
  std::vector<float> fMcElectronX;
  std::vector<float> fMcElectronY; 
  std::vector<float> fMcElectronZ;
  std::vector<float> fRecoSpX;
  std::vector<float> fRecoSpY; 
  std::vector<float> fRecoSpZ;
  int fNAllHits = 0;
  int fNSlices = 0;
  double fEkGen;
  std::vector<int> fNhitsPerSlice;
  std::vector<int> fChargePerSlice;
  std::vector<double> fGammaEk_MC;  

};


test::NeutronCaptureFinder::NeutronCaptureFinder(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)  ,
    fNPrimaries(0)
  // More initializers here.
{
	this->reconfigure(p); 
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::NeutronCaptureFinder::beginJob() {
	//create output tree
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("mytree", "My Tree");
  fTree->Branch("event", &fEvent, "event/I"); 
  fTree->Branch("trueEnergy", &fTrueEnergy, "trueEnergy/F");
  fTree->Branch("nPrimaries", &fNPrimaries, "nPrimaries/I");
  fTree->Branch("nParticleHits", &fParticleHits, "nParticleHits/I");
  fTree->Branch("nAllHits", &fNAllHits, "nParticleHits/I");
  fTree->Branch("nSlices", &fNSlices, "nSlices/I");
  fTree->Branch("mcElectronX", &fMcElectronX);
  fTree->Branch("mcElectronY", &fMcElectronY);
  fTree->Branch("mcElectronZ", &fMcElectronZ);
  fTree->Branch("recoSpX", &fRecoSpX);
  fTree->Branch("recoSpY", &fRecoSpY);
  fTree->Branch("recoSpZ", &fRecoSpZ);
  fTree->Branch("fGammaEk_MC", &fGammaEk_MC);
  fLengthHist = tfs->make<TH1F>("hLength", ";Reconstructed  Track Length (cm)", 50, 0, 400);
}

void test::NeutronCaptureFinder::endJob() {
	
}

void test::NeutronCaptureFinder::reconfigure(fhicl::ParameterSet const & p)
{
 // Implementation of optional member function here.
 fTruthLabel = p.get<std::string>("TruthLabel");
 //fParticleLabel = p.get<std::string>("ParticleLabel");
 //fTrackLabel = p.get<std::string>("TrackLabel"); 
 fHitLabel = p.get<std::string>("HitLabel");
 fSpacePointLabel = p.get<std::string>("SpacePointLabel");
 fSliceLabel = p.get<std::string>("SliceLabel");
 fSimEdepTag = (art::InputTag(fTruthLabel, p.get<std::string>("SimEdepInstanceLabel")));
}

void test::NeutronCaptureFinder::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  fEvent = e.id().event();
  // Access the MC truth information
  fTrueEnergy = -999.;
  if(!e.isRealData()) {
    // Get the list of MC particles from GEANT
    auto mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    std::cout<<"NeutronCaptureFinder::analyze(): "<<"Number of MCParticles = "<<mcParticles->size()<<std::endl;
    if(mcParticles.isValid()) {
    	int NMCElectron = 0;
    	
    	for(auto &trueParticle : *mcParticles) {	      
    		fEkGen = (std::sqrt(trueParticle.P()*trueParticle.P() + trueParticle.Mass()*trueParticle.Mass()) - trueParticle.Mass()) * 1000; // MeVs
    	  //const simb::MCParticle trueParticle = mcParticles->at(t);
    	  if(trueParticle.Process() == "primary") {
    	    fTrueEnergy = trueParticle.E();
    	    std::cout<<"PdgCode, Process, KineticE = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<std::endl;
    	  }
    	  
    	  // true electron position
    	  if(trueParticle.PdgCode() ==11) {
    	  	NMCElectron ++;
    	  	fMcElectronX.push_back(trueParticle.Vx());
    	  	fMcElectronY.push_back(trueParticle.Vy());
    	  	fMcElectronZ.push_back(trueParticle.Vz());
    	  }
    	  // neutron capture gamma position
        if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "primary") {
        	std::cout<<"gamma = ("<< trueParticle.Vx()<<", "<<trueParticle.Vy()<<", "<<trueParticle.Vz()<<"), energy = "<<fEkGen<<std::endl;
        	fGammaEk_MC.push_back(fEkGen);
        }
    	}
    	std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of MCElectrons = "<<NMCElectron<<std::endl;   
      	
        
    } 
                             
    // Access the reconstructed information.    
    
    // get the list of Hits
    auto recoHits = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    // get the list of SpacePoints
    auto recoSpacePoints = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
    // get the list of Slices
    auto recoSlices = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
    // get the list of PFParticles
    //auto recoParticles = e.getValidHandle<std::vector<recob::PFParticle>>(fParticleLabel);
    // get the list of tracks 
    //auto recoTracks = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
    
    if(recoHits.isValid()) {
    	std::cout<<"NeutronCaptureFinder::analyze(): "<<"Number of recoHits = "<<recoHits->size()<<std::endl;    
      fNAllHits = recoHits->size();	
    }
    
    if(recoSpacePoints.isValid()) {
    	std::cout<<"NeutronCaptureFinder::analyze(): "<<"Number of recoSpacePoints = "<<recoSpacePoints->size()<<std::endl;         
      for(const auto &sp : *recoSpacePoints) {
        const double* xyz = sp.XYZ();
        fRecoSpX.push_back(xyz[0]);     
        fRecoSpY.push_back(xyz[1]);     
        fRecoSpZ.push_back(xyz[2]);     
      }

    }
    
    if(recoSlices.isValid()) {
    	std::cout<<"NeutronCaptureFinder::analyze(): "<<"Number of recoSlices = "<<recoSlices->size()<<std::endl;       
    	fNSlices = recoSlices->size();
      // find spacepoints for slice
      const art::FindManyP<recob::SpacePoint> findSpacePoints(recoSlices, e, fSpacePointLabel);
      // find hits for slice

    }
    
//    if(recoParticles.isValid()) {
//    	// Get the associations between the particles and tracks
//      const art::FindManyP<recob::Track> findTracks(recoParticles, e, fTrackLabel);
//      // find hits for track
//      const art::FindManyP<recob::Hit> findHits(recoTracks, e, fTrackLabel);
//      fLength = 0.0;
//      fParticleHits = 0;
//    	// Now let's have a look through these individual particles
//      for(unsigned int p=0; p<recoParticles->size(); p++) {
//        const recob::PFParticle particle = recoParticles->at(p);
//        if(particle.IsPrimary()) {
//          fNPrimaries ++;
//          // Get the association between particle and the tracks
//          const std::vector<art::Ptr<recob::Track>> pfpTracks = findTracks.at(p);
//          // Check that we have one track and that the particle is pion-like (PDG code 211 = pi+)
//          if(pfpTracks.size() == 1 && particle.PdgCode() == 211) {
//            art::Ptr<recob::Track> thisTrack = pfpTracks.at(0);
//            fLength = thisTrack->Length();
//            // Get the vector of hits associated to the track
//            const std::vector<art::Ptr<recob::Hit>> trackHits = findHits.at(thisTrack->ID()); //need to use the index of the track which is recob::Track::ID()
//            fTrackHits = trackHits.size();
//            fParticleHits += fTrackHits;
//          }
//        }
//      }
//      
//    }
  }
  // Analysis goes here...
  fTree->Fill(); 
  fLengthHist->Fill(fLength);
}



DEFINE_ART_MODULE(test::NeutronCaptureFinder)
	

