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
  /*std::string fTruthLabel;
  std::string fParticleLabel;
  std::string fTrackLabel; 
  std::string fHitLabel;
  std::string fSpacePointLabel;
  std::string fSliceLabel;
  art::InputTag fSimEdepTag;     */
    std::string fTruthLabel;
    std::vector<art::InputTag> fEDepTags;
  	
  // Truth variables
  float	fTrueEnergy;
    double Etot;
    double Etotar;
    int NElectrons;
    
  // Reco variables
  int fNPrimaries;
  float fLength;
  TH1F *fLengthHist;
    TH1D *gamen;
    TH1D *gamenall;
    TH1D *gamenoth;
    TH1D *Etoth;
    TH1D *gamenar;
    TH1D *gamenar8;
    TH1D *Etotarh;
    TH1D *Eenh;
    Double_t ecode;
    //geo::Geometry* geo_serv;
    const TGeoMaterial * test_material;
    //TGeoElement* elem;
    geo::Point_t test_point;
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
    std::vector<double> ec;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> xa;
    std::vector<double> ya;
    std::vector<double> za;
    std::vector<double> xar;
    std::vector<double> yar;
    std::vector<double> zar;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    std::vector<double> xc;
    std::vector<double> yc;
    std::vector<double> zc;
    std::vector<double> Etotv;
    std::vector<double> gammaenc;
    std::vector<double> gammaenpe;
    std::vector<double> gammaennc;
    std::vector<double> gammaenar;
    std::vector<double> gammaenoth;
    std::vector<double> gammaenothall;
    std::vector<double> neutronenpri;
    std::vector<double> neutronenall;
    std::vector<double> enumber;
    //std::vector<double> vcode;

    art::ServiceHandle<geo::Geometry> geo_serv;
    
};


test::NeutronCaptureMC::NeutronCaptureMC(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)  ,
    fNPrimaries(0)
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
    fTree->Branch("ec",&ec);
    fTree->Branch("x",&x);
    fTree->Branch("y",&y);
    fTree->Branch("z",&z);
    fTree->Branch("xa",&xa);
    fTree->Branch("ya",&ya);
    fTree->Branch("za",&za);
    fTree->Branch("xar",&xar);
    fTree->Branch("yar",&yar);
    fTree->Branch("zar",&zar);
    fTree->Branch("xs",&xs);
    fTree->Branch("ys",&ys);
    fTree->Branch("zs",&zs);
    fTree->Branch("xc",&xc);
    fTree->Branch("yc",&yc);
    fTree->Branch("zc",&zc);
    fTree->Branch("gammaennc",&gammaennc);
    fTree->Branch("gammaenc",&gammaenc);
    fTree->Branch("gammaenpe",&gammaenpe);
    fTree->Branch("gammaenar",&gammaenar);
    fTree->Branch("gammaenoth",&gammaenoth);
    fTree->Branch("gammaenothall",&gammaenothall);
    fTree->Branch("neutronenpri",&neutronenpri);
    fTree->Branch("neutronenall",&neutronenall);
    fTree->Branch("enumber",&enumber);
  /*fTree->Branch("trueEnergy", &fTrueEnergy, "trueEnergy/F");
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
  fTree->Branch("fGammaEk_MC", &fGammaEk_MC);*/
  fLengthHist = tfs->make<TH1F>("hLength", ";Reconstructed  Track Length (cm)", 50, 0, 400);
    gamen = tfs->make<TH1D>("gamma energy",";gamma energy (MeV);Counts",60,0,20);
    gamenall = tfs->make<TH1D>("gamma energy all",";gamma energy (MeV);Counts",60,0,20);
    Etoth = tfs->make<TH1D>("nCapture total energy",";nCapture total energy(MeV);Counts",60,0,20);
    gamenar = tfs->make<TH1D>("gamma energy Ar",";gamma energy (MeV);Counts",60,0,20);
    Etotarh = tfs->make<TH1D>("nCapture total energy Ar",";nCapture total energy(MeV);Counts",60,0,20);
    gamenar8 = tfs->make<TH1D>("gamma energy Ar 8",";gamma energy (MeV);Counts",80,0,8);
    Eenh = tfs->make<TH1D>("electron energy",";electron energy (MeV);Counts",700,0,0.7);
    gamenoth = tfs->make<TH1D>("gamma energy other",";gamma energy (MeV);Counts",90,0,30);
}

void test::NeutronCaptureMC::endJob() {
    /*int xs=x.size();
    int ys=y.size();
    int zs=z.size();
    std::cout<<"x"<<std::endl;
    for(int i=0;i<xs;i++){
        std::cout<<x[i]<<std::endl;
    }
    std::cout<<"y"<<std::endl;
    for(int i=0;i<ys;i++){
        std::cout<<y[i]<<std::endl;
    }
    std::cout<<"z"<<std::endl;
    for(int i=0;i<zs;i++){
        std::cout<<z[i]<<std::endl;
    }
    
    int xas=xa.size();
    int yas=ya.size();
    int zas=za.size();
    std::cout<<"xa"<<std::endl;
    for(int i=0;i<xas;i++){
        std::cout<<xa[i]<<std::endl;
    }
    std::cout<<"ya"<<std::endl;
    for(int i=0;i<yas;i++){
        std::cout<<ya[i]<<std::endl;
    }
    std::cout<<"za"<<std::endl;
    for(int i=0;i<zas;i++){
        std::cout<<za[i]<<std::endl;
    }
    
    std::cout<<"Etot"<<std::endl;
    int Es=Etotv.size();
    for(int i=0;i<Es;i++){
        std::cout<<i<<" "<<Etotv[i]<<std::endl;
    }*/
    /*std::cout<<"code"<<std::endl;
    int cs=vcode.size();
    for(int i=0;i<cs;i++){
        std::cout<<i<<" "<<vcode[i]<<std::endl;
    }*/
    
}

void test::NeutronCaptureMC::reconfigure(fhicl::ParameterSet const & p)
{
 // Implementation of optional member function here.
 fTruthLabel = p.get<std::string>("TruthLabel");
 //fParticleLabel = p.get<std::string>("ParticleLabel");
 //fTrackLabel = p.get<std::string>("TrackLabel");
 //fHitLabel = p.get<std::string>("HitLabel");
 //fSpacePointLabel = p.get<std::string>("SpacePointLabel");
 //fSliceLabel = p.get<std::string>("SliceLabel");
 //fSimEdepTag = (art::InputTag(fTruthLabel, p.get<std::string>("SimEdepInstanceLabel")));
    //fSimEdepTag = (art::InputTag(fTruthLabel, p.get<std::string>("EDepModuleLabels")));
    fEDepTags = p.get<std::vector<art::InputTag>>("EDepModuleLabels");
}

void test::NeutronCaptureMC::analyze(art::Event const& e)
{
    //art::ServiceHandle<geo::Geometry> geo_serv;
  // Implementation of required member function here.
    Etot=0;
    Etotar=0;
  fEvent = e.id().event();
    ec.clear();
    x.clear();
    y.clear();
    z.clear();
    xa.clear();
    ya.clear();
    za.clear();
    xar.clear();
    yar.clear();
    zar.clear();
    xs.clear();
    ys.clear();
    zs.clear();
    xc.clear();
    yc.clear();
    zc.clear();
    gammaenc.clear();
    gammaenpe.clear();
    gammaennc.clear();
    gammaenar.clear();
    gammaenoth.clear();
    gammaenothall.clear();
    neutronenpri.clear();
    neutronenall.clear();
    enumber.clear();
    NElectrons=0;
  // Access the MC truth information
  //fTrueEnergy = -999.;!
  if(!e.isRealData()) {
      
      std::vector<std::vector<sim::SimEnergyDeposit> const*> edep_vecs;
      for (auto label : fEDepTags) {
        auto const& edep_handle = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(label);
        edep_vecs.push_back(edep_handle);
      }
      for (auto const& edeps : edep_vecs) { //loop over modules
        for (auto const& edep : *edeps) {
            NElectrons+=edep.NumElectrons();
        }
      }
      std::cout<<"NElectrons: "<<NElectrons<<std::endl;
      
    // Get the list of MC particles from GEANT
    auto mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of MCParticles = "<<mcParticles->size()<<std::endl;
    if(mcParticles.isValid()) {
    	int NMCElectron = 0;
    	
    	for(auto &trueParticle : *mcParticles) {	      
    		fEkGen = (std::sqrt(trueParticle.P()*trueParticle.P() + trueParticle.Mass()*trueParticle.Mass()) - trueParticle.Mass()) * 1000; // MeVs
    	  //const simb::MCParticle trueParticle = mcParticles->at(t);
    	  if(trueParticle.Process() == "primary"&&trueParticle.PdgCode() ==2112) {
    	    //fTrueEnergy = trueParticle.E();!
    	    //std::cout<<"PdgCode, Process, Total E(GeV), KineticE(MeV) = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<", "<<fEkGen<<std::endl;
              neutronenpri.push_back(fEkGen);
    	  }
            
            if(trueParticle.Process() == "primary") {

              std::cout<<"PdgCode, Process, Total E(GeV), KineticE(MeV) = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<", "<<fEkGen<<std::endl;
            }
            
            if(trueParticle.PdgCode() ==2112) {
                neutronenall.push_back(fEkGen);
            }
    	  
    	  // true electron position
    	  if(trueParticle.PdgCode() ==11) {
    	  	NMCElectron ++;
    	  	fMcElectronX.push_back(trueParticle.Vx());
    	  	fMcElectronY.push_back(trueParticle.Vy());
    	  	fMcElectronZ.push_back(trueParticle.Vz());
              //std::cout<<"PdgCode, Process, KineticE = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<std::endl;
    	  }
    	  // neutron capture gamma position
        if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "primary") {
        //if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "nCapture") {
        	//std::cout<<"gamma = ("<< trueParticle.Vx()<<", "<<trueParticle.Vy()<<", "<<trueParticle.Vz()<<"), energy = "<<fEkGen<<std::endl;
        	fGammaEk_MC.push_back(fEkGen);
        }
            
            
            
            if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "nCapture"){
                gamen->Fill(fEkGen);
                gammaennc.push_back(fEkGen);
                Etot+=fEkGen;
            }
            
            if(trueParticle.PdgCode() ==11 && trueParticle.Process() == "primary"){
                Eenh->Fill(fEkGen);
            }
            
            /*if(trueParticle.PdgCode() ==22){
                gamenall->Fill(fEkGen);
                std::cout<<"test point"<<std::endl;
                test_point.SetCoordinates(trueParticle.Vx(), trueParticle.Vy(), trueParticle.Vz());
                std::cout<<"test mat"<<std::endl;
                test_material = geo_serv->Material(test_point);
                //std::cout<<"elem"<<std::endl;
                //elem=test_material->GetBaseElement();
                std::cout<<"ecode"<<std::endl;
                ecode=test_material->GetZ();
                std::cout<<"vcode"<<std::endl;
                vcode.push_back(ecode);
            }*/
            
            if(trueParticle.PdgCode() ==22){
                gamenall->Fill(fEkGen);
                //std::cout<<"PdgCode, Process, KineticE = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<std::endl;
            }
            
            
            
            test_point.SetCoordinates(trueParticle.Vx(), trueParticle.Vy(), trueParticle.Vz());
            test_material = geo_serv->Material(test_point);
            ecode=test_material->GetZ();
            ec.push_back(ecode);
            //std::cout<<"ecode "<<ecode<<std::endl;
            
            if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "nCapture" && ecode==18){
                gamenar->Fill(fEkGen);
                gamenar8->Fill(fEkGen);
                gammaenar.push_back(fEkGen);
                Etotar+=fEkGen;
                xar.push_back(trueParticle.Vx());
                yar.push_back(trueParticle.Vy());
                zar.push_back(trueParticle.Vz());
            }
            
            if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "nCapture" && ecode==3.5){
                xs.push_back(trueParticle.Vx());
                ys.push_back(trueParticle.Vy());
                zs.push_back(trueParticle.Vz());
                gammaenpe.push_back(fEkGen);
            }
            
            if(trueParticle.PdgCode() ==22 && trueParticle.Process() == "nCapture" && ecode==6){
                xc.push_back(trueParticle.Vx());
                yc.push_back(trueParticle.Vy());
                zc.push_back(trueParticle.Vz());
                gammaenc.push_back(fEkGen);
            }
            
            if(trueParticle.PdgCode() ==22&&trueParticle.Process() == "nCapture" && ecode!=18){
                x.push_back(trueParticle.Vx());
                y.push_back(trueParticle.Vy());
                z.push_back(trueParticle.Vz());
                gamenoth->Fill(fEkGen);
                gammaenoth.push_back(fEkGen);
            }
            if(trueParticle.PdgCode() ==22){
                xa.push_back(trueParticle.Vx());
                ya.push_back(trueParticle.Vy());
                za.push_back(trueParticle.Vz());
                gammaenothall.push_back(fEkGen);
            }
            
            
            
    	}
        Etotv.push_back(Etot);
        if(Etot!=0){Etoth->Fill(Etot);}
        if(Etotar!=0){Etotarh->Fill(Etotar);}
    	std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of MCElectrons = "<<NMCElectron<<std::endl;
        enumber.push_back(NMCElectron);
        
    }    
                             
    // Access the reconstructed information.    
    
//    // get the list of Hits
//    auto recoHits = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
//    // get the list of SpacePoints
//    auto recoSpacePoints = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
//    // get the list of Slices
//    auto recoSlices = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
//    // get the list of PFParticles
//    //auto recoParticles = e.getValidHandle<std::vector<recob::PFParticle>>(fParticleLabel);
//    // get the list of tracks 
//    //auto recoTracks = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
//    
//    if(recoHits.isValid()) {
//    	std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of recoHits = "<<recoHits->size()<<std::endl;    
//      fNAllHits = recoHits->size();	
//    }
//    
//    if(recoSpacePoints.isValid()) {
//    	std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of recoSpacePoints = "<<recoSpacePoints->size()<<std::endl;         
//      for(const auto &sp : *recoSpacePoints) {
//        const double* xyz = sp.XYZ();
//        fRecoSpX.push_back(xyz[0]);     
//        fRecoSpY.push_back(xyz[1]);     
//        fRecoSpZ.push_back(xyz[2]);     
//      }
//
//    }
//    
//    if(recoSlices.isValid()) {
//    	std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of recoSlices = "<<recoSlices->size()<<std::endl;       
//    	fNSlices = recoSlices->size();
//      // find spacepoints for slice
//      const art::FindManyP<recob::SpacePoint> findSpacePoints(recoSlices, e, fSpacePointLabel);
//      // find hits for slice
//
//    }
    
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



DEFINE_ART_MODULE(test::NeutronCaptureMC)
	

