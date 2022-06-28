// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Data type includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larsim/IonizationScintillation/ISCalcSeparate.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"

// ROOT includes.
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
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

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

namespace label_gen{

  class labelgen : public art::EDAnalyzer{
  public:
 
    explicit labelgen(fhicl::ParameterSet const& pset);
    virtual ~labelgen();

    void beginJob();
    void beginRun(const art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt); 

  private:

    // Parameters in .fcl file
    std::string fRawDigitLabel;
    std::string fTPCInput;
    std::string fTPCInstance;
    std::string fTruthLabel;
    std::string fSimChannelProducerLabel;
    std::string fSimChannelProducerInstance;
    art::InputTag fSimChannelProducerTag;

    // Branch variables for tree
    TTree *fTree;
    unsigned int fEvent;
    unsigned int fRun;
    unsigned int fSubRun;


    // TPC
    unsigned int fNUCh;
    unsigned int fNVCh;
    unsigned int fNZCh;

    // find channel boundaries for each view
    unsigned int fUChanMin;
    unsigned int fUChanMax;
    unsigned int fVChanMin;
    unsigned int fVChanMax;
    unsigned int fZChanMin;
    unsigned int fZChanMax;
    unsigned int fNticks;

    unsigned int fNofAPA;
    unsigned int fChansPerAPA;


    //unsigned int fMinT, fMaxT, fMaxTimeRange; // unused

    std::vector<TH2S*> fTimeChanU;
    std::vector<TH2S*> fTimeChanV;
    std::vector<TH2S*> fTimeChanZ;
    std::vector<double> Z0;
    std::vector<double> Z1;
    std::vector<double> Z2;
    std::vector<double> Z3;
    std::vector<double> Z4;
    std::vector<double> Z5;
    int cn;
      int vi;
      
      std::vector<double> Z0l;
      std::vector<double> Z1l;
      std::vector<double> Z2l;
      std::vector<double> Z3l;
      std::vector<double> Z4l;
      std::vector<double> Z5l;
      int scid;
      int cnl;
      int vil;
      int mother;
      int mothertemp;
      int mothertemp2;
      int pdg;
      double fEkGen;
      int pdgnum;
      int Ek;
      
      std::map<int,int> getmother;
      std::map<int,int> getpdg;
      std::map<int,int> getPDGnum;
      std::map<int,double> getE;

    // define nADC counts for uncompressed vs compressed
    unsigned int nADC_uncompPed;

    geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());


 }; // RawEventDisplay

  //-----------------------------------------------------------------------

  labelgen::labelgen(fhicl::ParameterSet const& parameterSet): EDAnalyzer(parameterSet) {
    this->reconfigure(parameterSet);

  }

  //-----------------------------------------------------------------------


  labelgen::~labelgen() {}
   

  //-----------------------------------------------------------------------
  void labelgen::beginJob() {
    // place to define the histograms

    art::ServiceHandle<art::TFileService> tfs;
    //Histogram names and titles
    
      fTree = tfs->make<TTree>("mytree", "My Tree");
      fTree->Branch("event", &fEvent, "event/I");
      fTree->Branch("Z0", &Z0);
      fTree->Branch("Z1", &Z1);
      fTree->Branch("Z2", &Z2);
      fTree->Branch("Z3", &Z3);
      fTree->Branch("Z4", &Z5);
      fTree->Branch("Z5", &Z5);
      fTree->Branch("Z0l", &Z0l);
      fTree->Branch("Z1l", &Z1l);
      fTree->Branch("Z2l", &Z2l);
      fTree->Branch("Z3l", &Z3l);
      fTree->Branch("Z4l", &Z5l);
      fTree->Branch("Z5l", &Z5l);
      
      //getPDGnum.insert(std::pair<int,int>(2112,1));
      getPDGnum.insert(std::pair<int,int>(11,4));
      getPDGnum.insert(std::pair<int,int>(-11,4));
      getPDGnum.insert(std::pair<int,int>(13,5));
      getPDGnum.insert(std::pair<int,int>(-13,5));
      
    std::stringstream  name, title;

    unsigned int UChMin;
    unsigned int UChMax;
    unsigned int VChMin;
    unsigned int VChMax;
    unsigned int ZChMin;
    unsigned int ZChMax;
    TH2S* TempHisto;


    // Accquiring geometry data
    fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2;
    fChansPerAPA = fGeom->Nchannels()/fNofAPA;

    // taken from dune35t module a way to organise the channel mapping:
    // loop through channels in the first APA to find the channel boundaries for each view
    // will adjust for desired APA after
    fUChanMin = 0;
    fZChanMax = fChansPerAPA - 1;
    for ( unsigned int c = fUChanMin + 1; c < fZChanMax; c++ ){
      if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
        fVChanMin = c;
        fUChanMax = c - 1;
      }
      if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
        fZChanMin = c;
        fVChanMax = c-1;
      }
    }
    

    fNUCh=fUChanMax-fUChanMin+1;
    fNVCh=fVChanMax-fVChanMin+1;
    fNZCh=fZChanMax-fZChanMin+1;

      std::cout<<"fNofAPA: "<<fNofAPA<<" fChansPerAPA: "<<fChansPerAPA<<std::endl;
      std::cout<<"fUChanMin: "<<fUChanMin<<" fUChanMax: "<<fUChanMax<<" fVChanMin: "<<fVChanMin<<" fVChanMax: "<<fVChanMax<<" fZChanMin: "<<fZChanMin<<" fZChanMax: "<<fZChanMax<<std::endl;
      std::cout<<"fNUCh: "<<fNUCh<<" fNVCh: "<<fNVCh<<" fNZCh: "<<fNZCh<<std::endl;
    
    unsigned int minT = 0;
    unsigned int maxT = 0;
    minT = 0;
    maxT = fNticks;
    unsigned int binT = (maxT-minT);

    for(unsigned int i=0;i<fNofAPA;i++){
      UChMin=fUChanMin + i*fChansPerAPA;
      UChMax=fUChanMax + i*fChansPerAPA;
      VChMin=fVChanMin + i*fChansPerAPA;
      VChMax=fVChanMax + i*fChansPerAPA;
      ZChMin=fZChanMin + i*fChansPerAPA;
      ZChMax=fZChanMax + i*fChansPerAPA;

      // construct the histograms; TH2 constructors: ("Name", "Title", NxBin, xMin, xMax, NyBin, yMax, yMin)
      name.str("");
      name << "fTimeChanU";
      name <<  i;
      title.str("");
      title << "Time vs Channel(Plane U, APA";
      title << i<<")";
      TempHisto = tfs->make<TH2S>(name.str().c_str(),title.str().c_str(), UChMax - UChMin + 1, UChMin, UChMax, binT, minT, maxT);
      fTimeChanU.push_back(TempHisto);

      name.str("");
      name << "fTimeChanV";
      name << i;
      title.str("");
      title << "Time vs Channel(Plane V, APA";
      title << i<<")";
      TempHisto = tfs->make<TH2S>(name.str().c_str(),title.str().c_str(), VChMax - VChMin + 1, VChMin, VChMax, binT, minT, maxT);
      fTimeChanV.push_back(TempHisto);

      name.str("");
      name << "fTimeChanZ";
      name << i;
      title.str("");
      title << "Time vs Channel(Plane Z, APA";
      title <<i<<")";
      TempHisto = tfs->make<TH2S>(name.str().c_str(),title.str().c_str(), ZChMax - ZChMin + 1, ZChMin, ZChMax, binT, minT, maxT);
      fTimeChanZ.push_back(TempHisto);


      fTimeChanU[i]->SetStats(0);
      fTimeChanV[i]->SetStats(0);    
      fTimeChanZ[i]->SetStats(0);    


      fTimeChanU[i]->GetXaxis()->SetTitle("Channel"); fTimeChanU[i]->GetYaxis()->SetTitle("TDC");
      fTimeChanV[i]->GetXaxis()->SetTitle("Channel"); fTimeChanV[i]->GetYaxis()->SetTitle("TDC");
      fTimeChanZ[i]->GetXaxis()->SetTitle("Channel"); fTimeChanZ[i]->GetYaxis()->SetTitle("TDC");
    }


  }

  //-----------------------------------------------------------------------

  void labelgen::beginRun(const art::Run& run) {
    // place to read databases or run independent info
  }

  //-----------------------------------------------------------------------

  void labelgen::reconfigure(fhicl::ParameterSet const& p){

    // reconfigure without recompiling
    // read in the parameters from the .fcl file
    // allows for interactive changes of the parameter values

    fTPCInput       = p.get< std::string >("TPCInputModule");
    fTPCInstance    = p.get< std::string >("TPCInstanceName");
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    fNticks         = detProp.NumberTimeSamples();
    fTruthLabel = p.get<std::string>("TruthLabel");
    fSimChannelProducerLabel=p.get<std::string>("SimChannelLabel");
    fSimChannelProducerInstance=p.get<std::string>("SimChannelInstance");
    return;
  }



  //-----------------------------------------------------------------------

  void labelgen::analyze(const art::Event& event) {

    // called once per event

    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    std::cout << "EventNumber = " << fEvent << std::endl;
      
    Z0.clear();
    Z1.clear();
    Z2.clear();
    Z3.clear();
    Z4.clear();
    Z5.clear();
      Z0l.clear();
      Z1l.clear();
      Z2l.clear();
      Z3l.clear();
      Z4l.clear();
      Z5l.clear();
      getE.clear();
      
      getmother.clear();
      getpdg.clear();
      
      for(int i=0;i<5760000;i++){
          Z0.push_back(0);
          Z1.push_back(0);
          Z2.push_back(0);
          Z3.push_back(0);
          Z4.push_back(0);
          Z5.push_back(0);
          Z0l.push_back(0);
          Z1l.push_back(0);
          Z2l.push_back(0);
          Z3l.push_back(0);
          Z4l.push_back(0);
          Z5l.push_back(0);
      }

    // Get the objects holding raw information: RawDigit for TPC data
    art::Handle< std::vector<raw::RawDigit> > RawTPC;
    event.getByLabel(fTPCInput, fTPCInstance, RawTPC);


    // Fill pointer vectors - more useful form for the raw data
    // a more usable form
    std::vector< art::Ptr<raw::RawDigit> > RawDigits;
    art::fill_ptr_vector(RawDigits, RawTPC);



    // Loop over all RawDigits (entire channels)                                                                                                        
    for(auto const & dptr : RawDigits) {
      const raw::RawDigit & digit = *dptr;
      
      // Get the channel number for this digit
      uint32_t chan = digit.Channel();
      // number of samples in uncompressed ADC
      int nSamples = digit.Samples();
      unsigned int apa = std::floor( chan/fChansPerAPA );	  
      int pedestal = (int)digit.GetPedestal();
      
      std::vector<short> uncompressed(nSamples);
      // with pedestal	  
      raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());
      // subtract pedestals
      std::vector<short> uncompPed(nSamples);
      for (int i=0; i<nSamples; i++) uncompPed.at(i)=uncompressed.at(i)-pedestal;
      
      // number of ADC uncompressed without pedestal
      nADC_uncompPed=uncompPed.size();	  
      

      //Induction Plane	  
      if( fGeom->View(chan) == geo::kU){	
				for(unsigned int l=0;l<nADC_uncompPed;l++) {
	  			if(uncompPed.at(l)!=0){
	    			fTimeChanU[apa]->Fill(chan,l, uncompPed.at(l));
                    //std::cout<<"U APA: "<<apa<<" chan: "<<chan<<" l: "<<l<<" uncompPed.at(l): " <<uncompPed.at(l)<<" chan-apa*fChansPerAPA: "<<chan-apa*fChansPerAPA<<std::endl;
	  			}
				}
      }// end of U View
      
      //Induction Plane	  
      if( fGeom->View(chan) == geo::kV){	
				for(unsigned int l=0;l<nADC_uncompPed;l++) {
	  			if(uncompPed.at(l)!=0){
	    			fTimeChanV[apa]->Fill(chan,l, uncompPed.at(l));
                    //std::cout<<"V APA: "<<apa<<" chan: "<<chan<<" l: "<<l<<" uncompPed.at(l): " <<uncompPed.at(l)<<" chan-apa*fChansPerAPA-fNUCh: "<<chan-apa*fChansPerAPA-fNUCh<<std::endl;
	  			}
				}
      }// end of V View

      if ( fGeom->View(chan) == geo::kZ){
				for(unsigned int l=0;l<nADC_uncompPed;l++) {
	  			if(uncompPed.at(l)!=0){
	    			fTimeChanZ[apa]->Fill(chan,l, uncompPed.at(l));
                    cn=chan-apa*fChansPerAPA-fNUCh-fNVCh;
                    vi=6000*cn+l;
                    if(apa==0){
                        Z0[vi]=uncompPed.at(l);
                    }
                    else if(apa==1){
                        Z1[vi]=uncompPed.at(l);
                    }
                    else if(apa==2){
                        Z2[vi]=uncompPed.at(l);
                    }
                    else if(apa==3){
                        Z3[vi]=uncompPed.at(l);
                    }
                    else if(apa==4){
                        Z4[vi]=uncompPed.at(l);
                    }
                    else if(apa==5){
                        Z5[vi]=uncompPed.at(l);
                    }

                //std::cout<<"Z APA: "<<apa<<" chan: "<<chan<<" l: "<<l<<" uncompPed.at(l): " <<uncompPed.at(l)<<" chan-apa*fChansPerAPA-fNUCh-fNVCh: "<<chan-apa*fChansPerAPA-fNUCh-fNVCh<<std::endl;
	  			}
				}	
      }
    } // RawDigits
      
    fSimChannelProducerTag = art::InputTag(fSimChannelProducerLabel, fSimChannelProducerInstance);
    auto scs = event.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelProducerTag);
    auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    for(auto &trueParticle : *mcParticles) {
        fEkGen = (std::sqrt(trueParticle.P()*trueParticle.P() + trueParticle.Mass()*trueParticle.Mass()) - trueParticle.Mass()) * 1000; // MeVs
        getmother.insert(std::pair<int,int>(trueParticle.TrackId(),trueParticle.Mother()));
        getpdg.insert(std::pair<int,int>(trueParticle.TrackId(),trueParticle.PdgCode()));
        getE.insert(std::pair<int,double>(trueParticle.TrackId(),fEkGen));
    }
    for(auto &sc : *scs){
        auto simChannelNumber = sc.Channel();
        if ( fGeom->View(simChannelNumber) != geo::kZ) continue;
        unsigned int apal = std::floor( simChannelNumber/fChansPerAPA );
        std::cout<<"apal: "<<apal<<std::endl;
        //std::cout<<"simChannelNumber: "<<simChannelNumber<<std::endl;
        for(int i=0;i<6000;i++){
            auto const& trackInfo=sc.TrackIDEs(i, i);
            //int infosize=trackInfo.size();
            if((int)trackInfo.size()!=0){
                //std::cout<<"chan: "<<simChannelNumber<<" i: "<<i<<" infosize: "<<(int)trackInfo.size()<<std::endl;
                scid=trackInfo[0].trackID;
                cnl=simChannelNumber-apal*fChansPerAPA-fNUCh-fNVCh;
                vil=6000*cnl+i;
                for(auto &trueParticle : *mcParticles) {
                    auto mcid=trueParticle.TrackId();
                    if (mcid != scid) continue;
                    //fEkGen = (std::sqrt(trueParticle.P()*trueParticle.P() + trueParticle.Mass()*trueParticle.Mass()) - trueParticle.Mass()) * 1000; // MeVs
                    mother=trueParticle.Mother();
                    mothertemp=mcid;
                    while (mother != 0)
                    {
                        mothertemp2=mothertemp;
                        mothertemp=mother;
                        mother=getmother[mother];
                    }
                    pdg=getpdg[mothertemp];
                    if(pdg==2112){
                        Ek=getE[mothertemp2];
                        if(Ek<4){
                            pdgnum=1;
                        }
                        else if(Ek>=4 && Ek<5){
                            pdgnum=2;
                        }
                        else{
                            pdgnum=3;
                        }
                    }
                    else{
                        //pdg=getpdg[mothertemp];
                        pdgnum=getPDGnum[pdg];
                    }
                    if(apal==0){
                        Z0l[vil]=pdgnum;
                    }
                    else if(apal==1){
                        Z1l[vil]=pdgnum;
                    }
                    else if(apal==2){
                        Z2l[vil]=pdgnum;
                    }
                    else if(apal==3){
                        Z3l[vil]=pdgnum;
                    }
                    else if(apal==4){
                        Z4l[vil]=pdgnum;
                    }
                    else if(apal==5){
                        Z5l[vil]=pdgnum;
                    }
                    
                    break;
                }
                
            }
        }
    }
      
      
      
      fTree->Fill();
    return;
  }
  
}
DEFINE_ART_MODULE(label_gen::labelgen)
