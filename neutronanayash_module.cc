#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"


//#include "larexamples/Algorithms/RemoveIsolatedSpacePoints/SpacePointIsolationAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
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

#include <memory> // std::make_unique()

using namespace std;

struct hitStruct {
    float PT;
    float cID;
    //int PlnNum;

    bool operator==(const hitStruct &o) const {
        return PT == o.PT && cID == o.cID;
    }

    bool operator<(const hitStruct &o)  const {
        return PT < o.PT || (PT == o.PT && cID < o.cID);
    }
};

/*
struct sptStruct {
    float posY;
    float posZ;

    bool operator==(const sptStruct &o) const {
        return posY == o.posY && posZ == o.posZ;
    }

    bool operator<(const sptStruct &o)  const {
        return posY < o.posY || (posY == o.posY && posZ < o.posZ);
    }
};
*/

struct sptStruct {
    float posX;
    float posY;
    float posZ;

    bool operator==(const sptStruct &other) const{ 
        return (posX == other.posX
        && posY == other.posY
        && posZ == other.posZ);
    }
};

namespace std{

    template <>
    struct hash<sptStruct>
    {
        std::size_t operator()(const sptStruct& k) const
        {
            using std::size_t;
            using std::hash;
            //using std::float;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return ((hash<float>()(k.posX)
            ^ (hash<float>()(k.posY) << 1)) >> 1)
            ^ (hash<float>()(k.posZ) << 1);
        }
    };
}

namespace dune{

    class neutronanayash : public art::EDAnalyzer {
    public:

        explicit neutronanayash(fhicl::ParameterSet const& pset);
        virtual ~neutronanayash();

        void beginJob();
        void endJob();
        void beginRun(const art::Run& run);
        void analyze(const art::Event& evt);
        void reset();
    
    private:

        const detinfo::DetectorProperties* fDetprop;

        float GetEhitMeV(const recob::Hit & hit) const;

        double fElectronsToGeV;

        std::ofstream file;

        geo::GeometryCore const * fGeometry;

        TTree* fEventTree;
        TTree* fEventHitsTree;
        Int_t    run;                  
        Int_t    subrun;               
        Int_t    event;
        Double_t evttime;

        calo::CalorimetryAlg fCalorimetryAlg;
        std::string fHitsModuleLabel;
        std::string fTrackModuleLabel;

        //Spacepoints
        std::string fSpacePointTag;

        //shower
        std::string fShowerTag;

        //3dcluster
        std::string f3dclusterTag;

        std::string f2dclusterTag;
        
        std::string fCalorimetryModuleLabel;

        

        bool  fSaveCaloInfo;
        bool  fSaveTrackInfo;

        //// Pandora
        std::vector<int> HitChannelID;
        std::vector<float> HitPeakTime;
        std::vector<int> HitClusterTag; //0 = Cluster Hit ; 1 = not a cluster hit
        //std::vector<bool> HitVetoTag; //0 = Bad Hit ; 1 = good hit
        std::vector<int> HitPlaneNo;
        std::vector<float> HitEnergy;

        // Trajcluster
        std::vector<int> TrajHitTag; //0 = Cluster Hit ; 1 = not a cluster hit

        ////// DBScan

        std::vector<float> HitSPx;
        std::vector<float> HitSPy;
        std::vector<float> HitSPz;
        //std::vector<bool> SptTrkTag; // 1 = Track point; 0 = not a track point
        //std::vector<int> sptSize;
        std::vector<int> SptClusterSize;
        std::vector<int> SptClusterTag; //0 = bad Spt ; 1 = good Spt
        std::vector<int> SptHitChannelID;
        std::vector<float> SptHitPeakTime;
        std::vector<int> SptHitTag; //0 = bad Hit ; 1 = good hit
        std::vector<int> SptHitPlaneNo;
        std::vector<float> SptHitEnergy;

        std::vector<float> TrajPandoraSPx;
        std::vector<float> TrajPandoraSPy;
        std::vector<float> TrajPandoraSPz;
        std::vector<int> TrajSptTag; //0 = bad Spt ; 1 = good Spt
        std::vector<int> PandoraSptTag; //0 = bad Spt ; 1 = good Spt



    };

    //========================================================================
neutronanayash::neutronanayash(fhicl::ParameterSet const& pset) :
        EDAnalyzer(pset),
        fCalorimetryAlg           (pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
        fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
        fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
        fSpacePointTag            (pset.get<std::string>("SpacePointLabel","")    ),
        fShowerTag                (pset.get<std::string>("ShowerLabel","")    ),
        f3dclusterTag             (pset.get<std::string>("Cluster3dLabel","")    ),
        f2dclusterTag             (pset.get<std::string>("Cluster2dLabel","")    ),
        fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
        fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
        fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false))
    {
        if (fSaveTrackInfo == false) fSaveCaloInfo = false;
        // get a pointer to the geometry service provider
        /*
        fGeometry = &*(art::ServiceHandle<geo::Geometry>());
        fDetprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        file.open("dump.txt");
        */
    }
 
    //========================================================================
neutronanayash::~neutronanayash(){
    }
    //========================================================================

    void neutronanayash::beginJob(){
        std::cout<<"job begin..."<<std::endl;

        art::ServiceHandle<sim::LArG4Parameters> larParameters;
        fElectronsToGeV = 1./larParameters->GeVToElectrons();
        std::cout<<"fElectronsToGeV = "<<fElectronsToGeV<<std::endl;

        art::ServiceHandle<art::TFileService> tfs;
        fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
        fEventTree->Branch("event", &event,"event/I");
        fEventTree->Branch("evttime",&evttime,"evttime/D");
        fEventTree->Branch("run", &run,"run/I");
        fEventTree->Branch("subrun", &subrun,"surbrun/I");

        fEventTree->Branch("HitChannelID", &HitChannelID);
        fEventTree->Branch("HitPeakTime", &HitPeakTime);
        fEventTree->Branch("HitClusterTag", &HitClusterTag);
        //fEventTree->Branch("HitVetoTag", &HitVetoTag);
        fEventTree->Branch("HitPlaneNo", &HitPlaneNo);
        fEventTree->Branch("HitEnergy", &HitEnergy);
        fEventTree->Branch("TrajHitTag", &TrajHitTag);

        fEventTree->Branch("HitSPx", &HitSPx);
        fEventTree->Branch("HitSPy", &HitSPy);
        fEventTree->Branch("HitSPz", &HitSPz);
        //fEventTree->Branch("SptTrkTag", &SptTrkTag);
        //fEventTree->Branch("sptSize", &sptSize);
        fEventTree->Branch("SptClusterTag", &SptClusterTag);
        fEventTree->Branch("SptClusterSize", &SptClusterSize);
        fEventTree->Branch("SptHitChannelID", &SptHitChannelID);
        fEventTree->Branch("SptHitPeakTime", &SptHitPeakTime);
        fEventTree->Branch("SptHitTag", &SptHitTag);
        //fEventTree->Branch("HitVetoTag", &HitVetoTag);
        fEventTree->Branch("SptHitPlaneNo", &SptHitPlaneNo);
        fEventTree->Branch("SptHitEnergy", &SptHitEnergy);

        fEventTree->Branch("TrajPandoraSPx", &TrajPandoraSPx);
        fEventTree->Branch("TrajPandoraSPy", &TrajPandoraSPy);
        fEventTree->Branch("TrajPandoraSPz", &TrajPandoraSPz);
        fEventTree->Branch("TrajSptTag", &TrajSptTag);
        fEventTree->Branch("PandoraSptTag", &PandoraSptTag);
    }

    //========================================================================
    void neutronanayash::endJob(){

    }

    //========================================================================
    void neutronanayash::beginRun(const art::Run&){
        mf::LogInfo("neutronana")<<"begin run..."<<std::endl;
    }
    //========================================================================

    void neutronanayash::analyze( const art::Event& evt){
        art::Handle< std::vector<recob::Track> > trackListHandle;
        std::vector<art::Ptr<recob::Track> > tracklist;
        if(evt.getByLabel("pandoraTrack",trackListHandle)){
          art::fill_ptr_vector(tracklist, trackListHandle);
        }
        else return;

        art::Handle< std::vector<recob::SpacePoint> > pointsListHandle; // to get information about the Spacepoints
        std::vector<art::Ptr<recob::SpacePoint> > pointsList;
        if(evt.getByLabel(fSpacePointTag, pointsListHandle)){
          art::fill_ptr_vector(pointsList, pointsListHandle);
        }

/*        
        art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
        std::vector<art::Ptr<recob::PFParticle> > pfplist;
        if(evt.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
*/      
        
        art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
        std::vector<art::Ptr<recob::Hit>> hitlist;
        if(evt.getByLabel(fHitsModuleLabel, hitListHandle)){
            art::fill_ptr_vector(hitlist, hitListHandle);
        }

        art::Handle< std::vector<recob::Slice> > sliceListHandle; // to get information about the slices
        std::vector<art::Ptr<recob::Slice>> sliceList;
        if(evt.getByLabel(f3dclusterTag, sliceListHandle)){
            art::fill_ptr_vector(sliceList, sliceListHandle);
        }

        art::Handle< std::vector<recob::Cluster> > trajListHandle; // to get information about the clusters
        std::vector<art::Ptr<recob::Cluster>> trajList;
        if(evt.getByLabel(f2dclusterTag, trajListHandle)){
            art::fill_ptr_vector(trajList, trajListHandle);
        }

        //Associations
        art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
        art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
        // art::FindManyP<recob::Track, recob::TrackHitMeta> thass(hitListHandle, evt, fTrackModuleLabel); //to associate hit 
        art::FindManyP<recob::Track> thass(hitListHandle, evt, fTrackModuleLabel); //to associate hit just trying
        art::FindManyP<recob::Shower> hitsToShower(hitListHandle, evt, fShowerTag); //to associate hit to showers

        art::FindManyP<recob::SpacePoint> spsFromHits(hitListHandle, evt, fHitsModuleLabel); //to associate hit with space point
        art::FindManyP<recob::Hit> hitsFromSps(pointsListHandle, evt, fSpacePointTag); //to associate space point to hit
        art::FindManyP<recob::SpacePoint> spsFromSlice(sliceListHandle, evt, f3dclusterTag); //to associate slices to spacepoint
        //art::FindManyP<recob::Cluster> clusterFromHit(hitListHandle, evt, f2dclusterTag); //to associate clusters to hit
        art::FindManyP<recob::Hit> hitsFromTraj(trajListHandle, evt, f2dclusterTag); //to associate clusters to hit

        art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
//        art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
        
//        art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
        art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");

        run = evt.run();
        subrun = evt.subRun();
        event = evt.id().event();
        art::Timestamp ts = evt.time();
        TTimeStamp tts(ts.timeHigh(), ts.timeLow());
        evttime=tts.AsDouble();

        size_t NTracks = tracklist.size();
        cout << "Total number of tracks : " << NTracks << endl;

        cout << "Total number of hits : " << hitlist.size() << endl;
        //int l = 0; //To count total non-track hits

        HitChannelID.clear();
        HitPeakTime.clear();
        HitClusterTag.clear();
        //HitVetoTag.clear();
        HitPlaneNo.clear();
        HitEnergy.clear();
        TrajHitTag.clear();

        HitSPx.clear();
        HitSPy.clear();
        HitSPz.clear();
        SptClusterTag.clear();
        SptClusterSize.clear();

        SptHitChannelID.clear();
        SptHitPeakTime.clear();
        SptHitPlaneNo.clear();
        SptHitEnergy.clear();
        SptHitTag.clear();
        //SptTrkTag.clear();
        //sptSize.clear();

        TrajPandoraSPx.clear();
        TrajPandoraSPy.clear();
        TrajPandoraSPz.clear();
        TrajSptTag.clear();
        PandoraSptTag.clear();

        /////////////////////////// Storing Space Points /////////////////////////////////////

        std::unordered_map<sptStruct, int> SptClusterMap;

        for (int i = 0; i < (int) sliceList.size(); ++i) {

            auto & sliceSps = spsFromSlice.at(i);

            for (int j = 0; j < (int) sliceSps.size(); ++j){
                const double* xyz = sliceSps[j]->XYZ();

                sptStruct spt;
                spt.posX = xyz[0];
                spt.posY = xyz[1];
                spt.posZ = xyz[2];

                SptClusterMap.emplace(spt, sliceSps.size());
            }
        }

        for (int i = 0; i < (int) pointsList.size(); ++i){

            auto & spsHit = hitsFromSps.at(i);
            const double* xyz = pointsList[i]->XYZ();

            int TAG = 0;

            HitSPx.push_back(xyz[0]);
            HitSPy.push_back(xyz[1]);
            HitSPz.push_back(xyz[2]);

            sptStruct spt;
            spt.posX = xyz[0];
            spt.posY = xyz[1];
            spt.posZ = xyz[2];

            std::unordered_map<sptStruct, int>::iterator itr = SptClusterMap.find(spt);
            if (itr != SptClusterMap.end())
            {
                SptClusterTag.push_back(0); //Cluster space point
                SptClusterSize.push_back(itr->second);
                TAG = 0; //Bad Hit
            } else {
                SptClusterTag.push_back(1); //Non-Cluster space point
                SptClusterSize.push_back(0);
                TAG = 1; //Good Hit
            }

            for (int j = 0; j < (int) spsHit.size(); ++j){
                SptHitChannelID.push_back(spsHit[j]->Channel());
                SptHitPeakTime.push_back(spsHit[j]->PeakTime());
                SptHitPlaneNo.push_back(spsHit[j]->WireID().Plane);
                SptHitEnergy.push_back(GetEhitMeV(*spsHit[j]));
                SptHitTag.push_back(TAG);
            }
        }

/*
        for (int i = 0; i < (int) sliceList.size(); ++i){

            auto & sliceSps = spsFromSlice.at(i);

            for (int j = 0; j < (int) sliceSps.size(); ++j){
                const double* xyz = sliceSps[j]->XYZ();
                HitSPx.push_back(xyz[0]);
                HitSPy.push_back(xyz[1]);
                HitSPz.push_back(xyz[2]);
                SptClusterTag.push_back(sliceSps.size());
            }
        }
*/

        ///////////////////////////////////////////////////////////////////// Trajcluster
        std::map<hitStruct, int> TrajHitMap;

        for(int i = 0; i < (int) trajList.size(); ++i){

            auto & TrajHit = hitsFromTraj.at(i);

            for (int j = 0; j < (int) TrajHit.size(); ++j)
            {   
                hitStruct vhit;
                vhit.PT = TrajHit[j]->PeakTime();
                vhit.cID = TrajHit[j]->Channel();

                TrajHitMap.emplace(vhit, TrajHit.size());
            }
        }
        

        ////////////////////////////////////////////////////// Storing Hits ////////////////////////////////////////////////

        //std::map<hitStruct, bool> HitTrkMap;

        for(size_t j=0;j< hitlist.size();j++){ //loop over all the hits

            auto & Tracks = thass.at(hitlist[j].key());
            auto & Shower = hitsToShower.at(hitlist[j].key());

            auto & sps = spsFromHits.at(hitlist[j].key());

            int TrajTag = 0;
            int PandoraTag = 0;

            hitStruct vhit;
            vhit.PT = hitlist[j]->PeakTime();
            vhit.cID = hitlist[j]->Channel();

            HitChannelID.push_back(hitlist[j]->Channel());
            HitPeakTime.push_back(hitlist[j]->PeakTime());
            HitEnergy.push_back(GetEhitMeV(*hitlist[j]));
            HitPlaneNo.push_back(hitlist[j]->WireID().Plane);

            //filling the map
            if(Tracks.empty() && Shower.empty()){
                HitClusterTag.push_back(1); // Not a Track or a Shower hit
                PandoraTag = 1;
                //HitTrkMap.emplace(vhit, 1); // Not a Track or a Shower hit
            } else {    

                HitClusterTag.push_back(0); // A Track or a Shower hit
                PandoraTag = 0;
                //HitTrkMap.emplace(vhit, 0); // A Track or a Shower hit
            }

            std::map<hitStruct, int>::iterator itr = TrajHitMap.find(vhit);
            if (itr != TrajHitMap.end())
            {
                TrajHitTag.push_back(0); //Cluster Hit
                TrajTag = 0;
            } else {
                TrajHitTag.push_back(1); //Non-Cluster Hit
                TrajTag = 1;
            }

            for (int i = 0; i < (int) sps.size(); ++i){

                const double* xyz = sps[i]->XYZ();

                TrajPandoraSPx.push_back(xyz[0]);
                TrajPandoraSPy.push_back(xyz[1]);
                TrajPandoraSPz.push_back(xyz[2]);

                TrajSptTag.push_back(TrajTag);
                PandoraSptTag.push_back(PandoraTag);
            }
        }
        
/*
        for(size_t j=0;j< hitlist.size();j++){ //loop over all the hits

            bool VetoTag = 1; //Not a track/shower hit

            std::map<hitStruct, bool>::iterator itr; // = HitTrkMap.find(vhit);

            //loop over track and shower hits
            for(itr = HitTrkMap.begin(); itr != HitTrkMap.end(); ++itr){
                if (itr->second == 0)
                {
                    hitStruct vhit = itr->first;
                    if(vhit.cID - 21 < hitlist[j]->Channel() && vhit.cID + 21 > hitlist[j]->Channel() && vhit.PT - 250 < hitlist[j]->PeakTime() && vhit.PT + 250 > hitlist[j]->PeakTime()){
                        HitVetoTag.push_back(0); //vetoed because too close to track or shower hit
                        VetoTag = 0;
                        break;
                    }
                }
            }

            if(VetoTag == 1){
                HitVetoTag.push_back(1);
            }
        }

        for (int i = 0; i < (int) pointsList.size(); ++i)
        {
            auto & spsHit = hitsFromSps.at(i);
            const double* xyz = pointsList[i]->XYZ();

            bool TrkPoint = 1;

            for (int j = 0; j < (int) spsHit.size(); ++j)
            {   
                hitStruct vhit;
                vhit.PT = spsHit[j]->PeakTime();
                vhit.cID = spsHit[j]->Channel();
                //vhit.PlnNum = spsHit[j]->WireID().Plane;

                std::map<hitStruct, bool>::iterator itr = HitTrkMap.find(vhit);
                if (itr != HitTrkMap.end() && itr->second == 0)
                {
                    TrkPoint = 0;
                    break;
                }
            }

            SptTrkTag.push_back(TrkPoint);
            HitSPx.push_back(xyz[0]);
            HitSPy.push_back(xyz[1]);
            HitSPz.push_back(xyz[2]);
            sptSize.push_back(spsHit.size());
        }

        */

        fEventTree->Fill();
        
    }

    float neutronanayash::GetEhitMeV(const recob::Hit & hit) const
    {
        float dqadc = hit.Integral();
        if (!std::isnormal(dqadc) || (dqadc < 0)) dqadc = 0.0;

        unsigned short plane = hit.WireID().Plane;
        // double tdrift = hit.PeakTime();
        float dqel = fCalorimetryAlg.ElectronsFromADCArea(dqadc, plane);

        float correllifetime = 1; // fCalorimetryAlg.LifetimeCorrection(tdrift, fT0);
        float dq = dqel * correllifetime * fElectronsToGeV * 1000;
        //double dq = dqel * correllifetime;
        if (!std::isnormal(dq) || (dq < 0)) dq = 0.0;   

        return dq; 
    } 

    DEFINE_ART_MODULE(neutronanayash)
}
