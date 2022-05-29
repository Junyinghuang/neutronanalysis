// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

// C++ includes
#include <cmath>
#include <map>

namespace
{
  //double GetEdepHitsMeV(const std::vector< recob::Hit > & hits) const;

  bool TDCIDETimeCompare(const sim::TDCIDE&, const sim::TDCIDE&);

}

namespace lar {
  namespace example{
    
    class MyEnergyAnalysis : public art::EDAnalyzer {
      
    public:
      
      explicit MyEnergyAnalysis(fhicl::ParameterSet const & p);

      virtual void beginJob() override;

      virtual void beginRun(const art::Run& run) override;

      virtual void analyze(const art::Event& event) override;

    private:
      
      double GetEdepHitsMeV(const std::vector<recob::Hit> & hits);
      double GetElectrons(const std::vector<recob::Hit> & hits);

      calo::CalorimetryAlg fCalorimetryAlg;
      std::string fHitProducerLabel;
      std::string fClusterProducerLabel;
      std::string fSimChannelProducerLabel;
      std::string fSimChannelProducerInstance;

      art::InputTag fSimChannelProducerTag;

      TTree* fReconstructionNtuple;

      int fEvent;
      int fRun;
      int fSubRun;
      double fEnergyFromHits;
      double fnElectrons;
      double fAllChannelEnergy;
      double fT0 = 0;


      geo::GeometryCore const* fGeometryService;
      double fElectronsToGeV;
      int fTriggerOffset;

    };

    MyEnergyAnalysis::MyEnergyAnalysis(fhicl::ParameterSet const & p):
      EDAnalyzer(p),
      fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg")),
      fHitProducerLabel(p.get<std::string>("HitLabel")),
      fClusterProducerLabel(p.get<std::string>("ClusterLabel")),
      fSimChannelProducerLabel(p.get<std::string>("SimChannelLabel")),
      fSimChannelProducerInstance(p.get<std::string>("SimChannelInstance"))
    {
      
      fGeometryService = lar::providerFrom<geo::Geometry>();

    }

    void MyEnergyAnalysis::beginJob()
    {
      art::ServiceHandle<art::TFileService const> tfs;
      
      fReconstructionNtuple =
	tfs->make<TTree>("MyEnergyAnalysisReconstruction", "MyEnergyAnalysisReconstruction");

      fReconstructionNtuple->Branch("Event", &fEvent, "Event/I");
      fReconstructionNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
      fReconstructionNtuple->Branch("Run", &fRun, "Run/I");
      fReconstructionNtuple->Branch("EnergyFromHits", &fEnergyFromHits, "EnegyFromHits/D");
      fReconstructionNtuple->Branch("AllChannelEnergy", &fAllChannelEnergy, "AllChannelEnergy/D");
      fReconstructionNtuple->Branch("nElectrons", &fnElectrons, "nElectrons/D");
    }

    void MyEnergyAnalysis::beginRun(const art::Run& /*run*/)
    {
      art::ServiceHandle<sim::LArG4Parameters const> larParameters;
      fElectronsToGeV = 1. /larParameters->GeVToElectrons();
    }

    void MyEnergyAnalysis::analyze(const art::Event& event)
    {

      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      double e = 0;

      fSimChannelProducerTag = art::InputTag(fSimChannelProducerLabel, fSimChannelProducerInstance);

      const auto hitsHandle = event.getValidHandle< std::vector<recob::Hit> >(fHitProducerLabel);
      std::cout<<"Found "<<hitsHandle->size()<<" hits in total"<<std::endl;

      auto channelHandle = event.getValidHandle< std::vector<sim::SimChannel> >(fSimChannelProducerTag);
      auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      //auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();

      fTriggerOffset = trigger_offset(clock_data);

      const auto& hitListHandle = *event.getValidHandle< std::vector<recob::Hit> >(fHitProducerLabel);
      fEnergyFromHits = GetEdepHitsMeV(hitListHandle);
      fnElectrons = GetElectrons(hitListHandle);


      for(auto const& hit : (*hitsHandle)){
        
	auto hitChannelNumber = hit.Channel();
        
	typedef sim::SimChannel::StoredTDC_t TDC_t;
	TDC_t start_tdc = clock_data.TPCTick2TDC(hit.StartTick());
        TDC_t end_tdc = clock_data.TPCTick2TDC(hit.EndTick());
        TDC_t hitStart_tdc = clock_data.TPCTick2TDC(hit.PeakTime() - 3. * hit.SigmaPeakTime());
        TDC_t hitEnd_tdc = clock_data.TPCTick2TDC(hit.PeakTime() + 3. * hit.SigmaPeakTime());

	start_tdc = std::max(start_tdc, hitStart_tdc);
	end_tdc = std::max(end_tdc, hitEnd_tdc);

	for(auto const& channel : (*channelHandle)){
	  
	  auto simChannelNumber = channel.Channel();
	  if (simChannelNumber != hitChannelNumber) continue;

	  auto const& timeSlices = channel.TDCIDEMap();

	  sim::TDCIDE startTime;
	  sim::TDCIDE endTime;
	  startTime.first = start_tdc;
	  endTime.first = end_tdc;

	  auto const startPointer = 
	    std::lower_bound(timeSlices.begin(), timeSlices.end(), startTime, TDCIDETimeCompare);

	  auto const endPointer =
            std::upper_bound(startPointer, timeSlices.end(), endTime, TDCIDETimeCompare);

          if (startPointer == timeSlices.end() || startPointer == endPointer) continue;
          MF_LOG_DEBUG("AnalysisExample")
            << "Time slice start = " << (*startPointer).first << std::endl;

          for (auto slicePointer = startPointer; slicePointer != endPointer; slicePointer++) {
            auto const timeSlice = *slicePointer;
            auto time = timeSlice.first;

            MF_LOG_DEBUG("AnalysisExample")
              << "Hit index = " << hit.LocalIndex() << " channel number = " << hitChannelNumber
              << " start TDC tick = " << hit.StartTick() << " end TDC tick = " << hit.EndTick()
              << " peak TDC tick = " << hit.PeakTime()
              << " sigma peak time = " << hit.SigmaPeakTime()
              << " adjusted start TDC tick = " << clock_data.TPCTick2TDC(hit.StartTick())
              << " adjusted end TDC tick = " << clock_data.TPCTick2TDC(hit.EndTick())
              << " adjusted peak TDC tick = " << clock_data.TPCTick2TDC(hit.PeakTime())
              << " adjusted start_tdc = " << start_tdc << " adjusted end_tdc = " << end_tdc
              << " time = " << time << std::endl;

            // Loop over the energy deposits.
            auto const& energyDeposits = timeSlice.second;
            for (auto const& energyDeposit : energyDeposits) {
	      e += energyDeposit.numElectrons * fElectronsToGeV;
	    }

	  }

	}

      }


      fAllChannelEnergy = e*1000.;
      std::cout<<"fEnergyFromHits="<<fEnergyFromHits<<std::endl;

      std::cout<<"fAllChannelEnergy="<<fAllChannelEnergy<<std::endl;

      fReconstructionNtuple->Fill();


    } 

    double MyEnergyAnalysis::GetEdepHitsMeV(const std::vector< recob::Hit > & hits)
    {
      if (!hits.size()) return 0.0;
      unsigned short fBestView = 2;
      double dqsum = 0.0;
      for (size_t h = 0; h < hits.size(); ++h)
	{
	  double dqadc = hits[h].Integral();
	  if (!std::isnormal(dqadc) || (dqadc < 0)) continue;
	  
	  unsigned short plane = hits[h].WireID().Plane;
	  if (plane == fBestView)
	    {
	      //double tdrift = hits[h].PeakTime();
	      double dqel = fCalorimetryAlg.ElectronsFromADCArea(dqadc, plane);
	      
	      //double correllifetime = fCalorimetryAlg.LifetimeCorrection(tdrift, fT0);
	      double dq = dqel * fElectronsToGeV * 1000;
	      
	      if (!std::isnormal(dq) || (dq < 0)) continue;
	      
	      dqsum += dq;
	    } 
	}

      return dqsum; 
    }

    double MyEnergyAnalysis::GetElectrons(const std::vector< recob::Hit > & hits)
    {
      if (!hits.size()) return 0.0;
      unsigned short fBestView = 2;
      double dqelsum = 0.0;
      for (size_t h = 0; h < hits.size(); ++h)
	{
	  double dqadc = hits[h].Integral();
	  if (!std::isnormal(dqadc) || (dqadc < 0)) continue;
	  
	  unsigned short plane = hits[h].WireID().Plane;
	  if (plane == fBestView)
	    {
	      //double tdrift = hits[h].PeakTime();
	      double dqel = fCalorimetryAlg.ElectronsFromADCArea(dqadc, plane);
	      
	      //double correllifetime = fCalorimetryAlg.LifetimeCorrection(tdrift, fT0);
	      double dq = dqel * fElectronsToGeV * 1000;
	      
	      if (!std::isnormal(dq) || (dq < 0)) continue;
	      
	      dqelsum += dqel;
	    } 
	}

      return dqelsum; 
    }
    
    DEFINE_ART_MODULE(MyEnergyAnalysis)
      

  }
}
namespace{
 

  bool
  TDCIDETimeCompare(const sim::TDCIDE& lhs, const sim::TDCIDE& rhs)
  {
    return lhs.first < rhs.first;
  }

}





