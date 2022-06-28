#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "cetlib/pow.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1F.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "larcore/Geometry/Geometry.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include <map>
#include <list>
#include <memory>

namespace test {
  class dataextract;
}

class test::dataextract : public art::EDAnalyzer {
public:
  explicit dataextract(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
    dataextract(dataextract const&) = delete;
    dataextract(dataextract&&) = delete;
    dataextract& operator=(dataextract const&) = delete;
    dataextract& operator=(dataextract&&) = delete;

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
  std::string fSliceLabel;
  art::InputTag fSimEdepTag;     */
    std::string fTruthLabel;
    //std::string fSCLabel;
    std::string fSpacePointLabel;
    std::string fHitLabel;
    std::string fSimChannelProducerLabel;
    std::string fSimChannelProducerInstance;
    art::InputTag fSimChannelProducerTag;
    //std::vector<art::InputTag> fEDepTags;
    //std::vector<int> SCChannelID;
    std::vector<int> HitChannelID;
    std::vector<int> HitPeakTime;
    std::vector<int> HitTrackID;
    std::vector<float> HitInt;
    std::vector<float> HitADC;
    std::vector<float> HitSigmaPT;
    std::vector<float> HitRMS;
    std::vector<float> HitPA;
    std::vector<float> HitSigmaPA;
    std::vector<float> HitSigmaInt;
    std::vector<float> HitGOF;
    std::vector<int> HitAncestorPDGnum;
    std::vector<int> HitAncestor;
    std::vector<int> HitAncestorPDG;
    std::vector<int> spx;
    std::vector<int> spy;
    std::vector<int> spz;
    std::vector<int> spAncestorPDG;
    std::vector<float> spInt;
    std::vector<float> spADC;
    std::map<int,int> getmother;
    std::map<int,int> getpdg;
    std::map<int,int> getAncestorpdg;
    std::map<int,float> getHitInt;
    std::map<int,float> getHitADC;
    std::map<int,int>getPDGnum;
      
  // Truth variables
  float    fTrueEnergy;
  // Reco variables
  int fNPrimaries;
    //geo::Geometry* geo_serv;
    const TGeoMaterial * test_material;
    //TGeoElement* elem;
    geo::Point_t test_point;
  double fEkGen;
    //std::vector<double> vcode;

    art::ServiceHandle<geo::Geometry> geo_serv;
    
};


test::dataextract::dataextract(fhicl::ParameterSet const& p)
  : EDAnalyzer(p),
    fNPrimaries(0)
  // More initializers here.
{
    this->reconfigure(p);
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::dataextract::beginJob() {
    //create output tree
    //art::ServiceHandle<geo::Geometry> geo_serv;
    
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("mytree", "My Tree");
    fTree->Branch("event", &fEvent, "event/I");
    fTree->Branch("HitChannelID", &HitChannelID);
    fTree->Branch("HitTrackID", &HitTrackID);
    fTree->Branch("HitInt", &HitInt);
    fTree->Branch("HitADC", &HitADC);
    fTree->Branch("HitSigmaPT", &HitSigmaPT);
    fTree->Branch("HitRMS", &HitRMS);
    fTree->Branch("HitPA", &HitPA);
    fTree->Branch("HitSigmaPA", &HitSigmaPA);
    fTree->Branch("HitSigmaInt", &HitSigmaInt);
    fTree->Branch("HitGOF", &HitGOF);
    fTree->Branch("HitPeakTime", &HitPeakTime);
    fTree->Branch("HitAncestor", &HitAncestor);
    fTree->Branch("HitAncestorPDG", &HitAncestorPDG);
    //fTree->Branch("SCChannelID", &SCChannelID);
    fTree->Branch("spAncestorPDG", &spAncestorPDG);
    fTree->Branch("HitAncestorPDGnum", &HitAncestorPDGnum);
    fTree->Branch("spInt", &spInt);
    fTree->Branch("spADC", &spADC);
    fTree->Branch("spx", &spx);
    fTree->Branch("spy", &spy);
    fTree->Branch("spz", &spz);
    
    /*getPDGnum.insert(std::pair<int,int>(-1,0));
    getPDGnum.insert(std::pair<int,int>(2112,1));
    getPDGnum.insert(std::pair<int,int>(13,2));
    getPDGnum.insert(std::pair<int,int>(-13,2));
    getPDGnum.insert(std::pair<int,int>(11,3));
    getPDGnum.insert(std::pair<int,int>(-11,3));*/
}

void test::dataextract::endJob() {

}

void test::dataextract::reconfigure(fhicl::ParameterSet const & p)
{
 // Implementation of optional member function here.
 fTruthLabel = p.get<std::string>("TruthLabel");
 //fParticleLabel = p.get<std::string>("ParticleLabel");
 //fTrackLabel = p.get<std::string>("TrackLabel");
 fHitLabel = p.get<std::string>("HitLabel");
 //fSimEdepTag = (art::InputTag(fTruthLabel, p.get<std::string>("SimEdepInstanceLabel")));
    //fSimEdepTag = (art::InputTag(fTruthLabel, p.get<std::string>("EDepModuleLabels")));
    //fEDepTags = p.get<std::vector<art::InputTag>>("EDepModuleLabels");
    //fSCLabel =p.get<std::string>("SimChannelLabel");
    fSimChannelProducerLabel=p.get<std::string>("SimChannelLabel");
    fSimChannelProducerInstance=p.get<std::string>("SimChannelInstance");
    fSpacePointLabel=p.get<std::string>("SpacePointLabel");
}




void test::dataextract::analyze(art::Event const& e)
{
    //art::ServiceHandle<geo::Geometry> geo_serv;
  // Implementation of required member function here.
    
  fEvent = e.id().event();
    
    HitChannelID.clear();
    HitPeakTime.clear();
    HitTrackID.clear();
    HitAncestor.clear();
    HitAncestorPDG.clear();
    spAncestorPDG.clear();
    spInt.clear();
    HitInt.clear();
    HitADC.clear();
    spADC.clear();
    spx.clear();
    spy.clear();
    spz.clear();
    getmother.clear();
    getpdg.clear();
    getAncestorpdg.clear();
    getHitInt.clear();
    getHitADC.clear();
    HitSigmaPT.clear();
    HitRMS.clear();
    HitPA.clear();
    HitSigmaPA.clear();
    HitSigmaInt.clear();
    HitGOF.clear();
    HitAncestorPDGnum.clear();
    getPDGnum.clear();
    getPDGnum.insert(std::pair<int,int>(-1,0));
    getPDGnum.insert(std::pair<int,int>(2112,1));
    getPDGnum.insert(std::pair<int,int>(13,2));
    getPDGnum.insert(std::pair<int,int>(-13,2));
    getPDGnum.insert(std::pair<int,int>(11,3));
    getPDGnum.insert(std::pair<int,int>(-11,3));
    //SCChannelID.clear();
  // Access the MC truth information
  //fTrueEnergy = -999.;!
    fSimChannelProducerTag = art::InputTag(fSimChannelProducerLabel, fSimChannelProducerInstance);
    auto scs = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelProducerTag);
    //std::cout<<scs->size()<<std::endl;
    /*if(scs.isValid()) {
        for(auto &sc : *scs) {
            SCChannelID.push_back(sc.Channel());
        }
    }*/
    
    art::Handle< std::vector<recob::SpacePoint> > pointsListHandle; // to get information about the Spacepoints
    std::vector<art::Ptr<recob::SpacePoint> > pointsList;
    if(e.getByLabel(fSpacePointLabel, pointsListHandle)){
      art::fill_ptr_vector(pointsList, pointsListHandle);
    }
    
    art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if(e.getByLabel(fHitLabel, hitListHandle)){
        art::fill_ptr_vector(hitlist, hitListHandle);
    }
    
    art::FindManyP<recob::Hit> hitsFromSps(pointsListHandle, e, fSpacePointLabel);
    
    //auto hits = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    
    auto mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    //std::cout<<hits->size()<<std::endl;
    for(auto &trueParticle : *mcParticles) {
        getmother.insert(std::pair<int,int>(trueParticle.TrackId(),trueParticle.Mother()));
        getpdg.insert(std::pair<int,int>(trueParticle.TrackId(),trueParticle.PdgCode()));
    }
    
    
    //if(hits.isValid()) {
        //for(auto &hit : *hits) {
        for(size_t i=0;i< hitlist.size();i++) {
            if (hitlist[i]->WireID().Plane != 2) continue;
            HitChannelID.push_back(hitlist[i]->Channel());
            HitPeakTime.push_back((int)hitlist[i]->PeakTime());
            HitInt.push_back(hitlist[i]->Integral());
            getHitInt.insert(std::pair<int,float>(i,hitlist[i]->Integral()));
            HitADC.push_back(hitlist[i]->SummedADC());
            getHitADC.insert(std::pair<int,float>(i,hitlist[i]->SummedADC()));
            HitSigmaPT.push_back(hitlist[i]->SigmaPeakTime());
            HitRMS.push_back(hitlist[i]->RMS());
            HitPA.push_back(hitlist[i]->PeakAmplitude());
            HitSigmaPA.push_back(hitlist[i]->SigmaPeakAmplitude());
            HitSigmaInt.push_back(hitlist[i]->SigmaIntegral());
            HitGOF.push_back(hitlist[i]->GoodnessOfFit());
            auto hitChannelNumber = hitlist[i]->Channel();
            auto hitpt = (int)hitlist[i]->PeakTime();
            //std::cout<<hitChannelNumber<<std::endl;
            int hitid=-1;
            int mother=-1;
            int mothertemp=-1;
            //int ancestor=-1;
            int pdg=-1;
            int pdgnum=-1;
            
            for(auto &sc : *scs) {
                auto simChannelNumber = sc.Channel();
                if (simChannelNumber != hitChannelNumber) continue;
                auto const& trackInfo=sc.TrackIDEs(hitpt, hitpt);
                if(trackInfo.size()!=0){
                    HitTrackID.push_back(trackInfo[0].trackID);
                    hitid=trackInfo[0].trackID;
                    int infosize=trackInfo.size();
                    std::cout<<"chan: "<<simChannelNumber<<" i: "<<hitpt<<" infosize: "<<infosize<<std::endl;
                    for(int j=0;j<infosize;j++){
                        std::cout<<"j: "<<j<<std::endl;
                        std::cout<<"numElectrons: "<<trackInfo[i].numElectrons<<std::endl;
                        std::cout<<"energy: "<<trackInfo[i].energy<<std::endl;
                        int sid=trackInfo[i].trackID;
                        std::cout<<"track id: "<<sid<<std::endl;
                }
                }
                else{
                    HitTrackID.push_back(-1);
                    HitAncestor.push_back(-1);
                    HitAncestorPDG.push_back(-1);
                    getAncestorpdg.insert(std::pair<int,int>(i,-1));
                }
                /*std::cout<<"tracksize: "<<trackInfo.size()<<std::endl;
                for(int j = 0; j < (int) trackInfo.size(); ++j){
                    std::cout<<"trackid: "<<trackInfo[j].trackID<<std::endl;
                }*/
            }
            
            if(hitid!=-1){
                for(auto &trueParticle : *mcParticles) {
                    auto mcid=trueParticle.TrackId();
                    if (mcid != hitid) continue;
                    mother = trueParticle.Mother();
                    mothertemp=hitid;
                    while (mother != 0)
                    {
                        mothertemp=mother;
                        mother=getmother[mother];
                    }
                    pdg=getpdg[mothertemp];
                    HitAncestor.push_back(mothertemp);
                    HitAncestorPDG.push_back(pdg);
                    getAncestorpdg.insert(std::pair<int,int>(i,pdg));
                    break;
                }
            }
            pdgnum=getPDGnum[pdg];
            std::cout<<"PDG: "<<pdg<<" PDGnum: "<<pdgnum<<std::endl;
            HitAncestorPDGnum.push_back(pdgnum);
        }
    
    for (int i = 0; i < (int) pointsList.size(); i++){
        auto & spsHit = hitsFromSps.at(i);
        const double* xyz = pointsList[i]->XYZ();
        for (int j = 0; j < (int) spsHit.size(); j++){
            if (spsHit[j]->WireID().Plane != 2) continue;
            //std::cout<<spsHit[j]->WireID().Plane<<std::endl;
            for(size_t k=0;k< hitlist.size();k++) {
                if(spsHit[j]->Channel()!=hitlist[k]->Channel()||spsHit[j]->PeakTime()!=hitlist[k]->PeakTime()) continue;
                spAncestorPDG.push_back(getAncestorpdg[k]);
                spInt.push_back(getHitInt[k]);
                spADC.push_back(getHitADC[k]);
                spx.push_back(xyz[0]);
                spy.push_back(xyz[1]);
                spz.push_back(xyz[2]);
                break;
            }
        }
    }
    //}
    // Get the list of MC particles from GEANT
    //auto mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    /*std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of MCParticles = "<<mcParticles->size()<<std::endl;
    if(mcParticles.isValid()) {
        //int NMCElectron = 0;
        
        for(auto &trueParticle : *mcParticles) {
            fEkGen = (std::sqrt(trueParticle.P()*trueParticle.P() + trueParticle.Mass()*trueParticle.Mass()) - trueParticle.Mass()) * 1000; // MeVs
          //const simb::MCParticle trueParticle = mcParticles->at(t);
          if(trueParticle.Process() == "primary"&&trueParticle.PdgCode() ==2112) {
            //fTrueEnergy = trueParticle.E();!
            //std::cout<<"PdgCode, Process, Total E(GeV), KineticE(MeV) = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<", "<<fEkGen<<" Mother: "<<trueParticlestd::endl;
          }
            
            
              std::cout<<"PdgCode, Process, Total E(GeV), KineticE(MeV) = "<<trueParticle.PdgCode()<<", "<<trueParticle.Process()<<", "<<trueParticle.E()<<", "<<fEkGen<<" Mother: "<<trueParticle.Mother()<<" ID: "<<trueParticle.TrackId()<<" NumberDaughters: "<<trueParticle.NumberDaughters()<<std::endl;
            
            
          // neutron capture gamma position
            }

        }*/
        
                             
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
//        std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of recoHits = "<<recoHits->size()<<std::endl;
//      fNAllHits = recoHits->size();
//    }
//
//    if(recoSpacePoints.isValid()) {
//        std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of recoSpacePoints = "<<recoSpacePoints->size()<<std::endl;
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
//        std::cout<<"NeutronCaptureMC::analyze(): "<<"Number of recoSlices = "<<recoSlices->size()<<std::endl;
//        fNSlices = recoSlices->size();
//      // find spacepoints for slice
//      const art::FindManyP<recob::SpacePoint> findSpacePoints(recoSlices, e, fSpacePointLabel);
//      // find hits for slice
//
//    }
    
//    if(recoParticles.isValid()) {
//        // Get the associations between the particles and tracks
//      const art::FindManyP<recob::Track> findTracks(recoParticles, e, fTrackLabel);
//      // find hits for track
//      const art::FindManyP<recob::Hit> findHits(recoTracks, e, fTrackLabel);
//      fLength = 0.0;
//      fParticleHits = 0;
//        // Now let's have a look through these individual particles
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
  // Analysis goes here...
  fTree->Fill();
}


DEFINE_ART_MODULE(test::dataextract)
    

