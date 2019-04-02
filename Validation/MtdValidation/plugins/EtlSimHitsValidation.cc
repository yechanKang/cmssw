// -*- C++ -*-
//
// Package:    Validation/MtdValidation
// Class:      EtlSimHitsValidation
//
/**\class EtlSimHitsValidation EtlSimHitsValidation.cc Validation/MtdValidation/plugins/EtlSimHitsValidation.cc

 Description: ETL SIM hits validation

 Implementation:
     [Notes on implementation]
*/

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"


struct MTDHit {
  float energy;
  float time;
  float x;
  float y;
  float z;
};


class EtlSimHitsValidation : public DQMEDAnalyzer {

public:
  explicit EtlSimHitsValidation(const edm::ParameterSet&);
  ~EtlSimHitsValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  void bookHistograms(DQMStore::IBooker &,
		      edm::Run const&,
		      edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ------------ member data ------------

  const MTDGeometry* geom_;

  const std::string folder_;
  const std::string g4InfoLabel_;
  const std::string etlHitsCollection_;

  const float hitMinEnergy_;

  int eventCount_;

  edm::EDGetTokenT<edm::PSimHitContainer> etlSimHitsToken_;

  // --- histograms declaration

  MonitorElement* meNhits_[2];
  MonitorElement* meNtrkPerCell_[2];

  MonitorElement* meHitEnergy_[2];
  MonitorElement* meHitTime_[2];

  MonitorElement* meHitXlocal_[2];
  MonitorElement* meHitYlocal_[2];
  MonitorElement* meHitZlocal_[2];

  MonitorElement* meOccupancy_[2];

  MonitorElement* meHitX_[2];
  MonitorElement* meHitY_[2];
  MonitorElement* meHitZ_[2];
  MonitorElement* meHitPhi_[2];
  MonitorElement* meHitEta_[2];

  MonitorElement* meHitTvsE_[2];
  MonitorElement* meHitEvsPhi_[2];
  MonitorElement* meHitEvsEta_[2];
  MonitorElement* meHitTvsPhi_[2];
  MonitorElement* meHitTvsEta_[2];

};


// ------------ constructor and destructor --------------
EtlSimHitsValidation::EtlSimHitsValidation(const edm::ParameterSet& iConfig):
  geom_(nullptr),
  folder_(iConfig.getParameter<std::string>("folder")),
  g4InfoLabel_(iConfig.getParameter<std::string>("moduleLabelG4")),
  etlHitsCollection_(iConfig.getParameter<std::string>("etlSimHitsCollection")),
  hitMinEnergy_( iConfig.getParameter<double>("hitMinimumEnergy") ),
  eventCount_(0) {
  
  etlSimHitsToken_ = consumes <edm::PSimHitContainer> (edm::InputTag(std::string(g4InfoLabel_),
								     std::string(etlHitsCollection_)));

}

EtlSimHitsValidation::~EtlSimHitsValidation() {
}


// ------------ method called for each event  ------------
void EtlSimHitsValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;

  edm::LogInfo("EventInfo") << " Run = " << iEvent.id().run() << " Event = " << iEvent.id().event();

  edm::ESHandle<MTDGeometry> geometryHandle;
  iSetup.get<MTDDigiGeometryRecord>().get(geometryHandle);
  geom_ = geometryHandle.product();

  edm::Handle<edm::PSimHitContainer> etlSimHitsHandle;
  iEvent.getByToken(etlSimHitsToken_, etlSimHitsHandle);

  if( ! etlSimHitsHandle.isValid() ) {
    edm::LogWarning("DataNotFound") << "No ETL SIM hits found";
    return;
  }

  eventCount_++;
  
  std::unordered_map<uint32_t, MTDHit> m_etlHits[2];
  std::unordered_map<uint32_t, std::set<int> > m_etlTrkPerCell[2];

  // --- Loop over the BLT SIM hits
  for (auto const& simHit: *etlSimHitsHandle) {

    // --- Use only hits compatible with the in-time bunch-crossing
    if ( simHit.tof() < 0 || simHit.tof() > 25. ) continue;

    ETLDetId id = simHit.detUnitId();

    int idet = (id.zside()+1)/2;

    m_etlTrkPerCell[idet][id.rawId()].insert(simHit.trackId());

    auto simHitIt = m_etlHits[idet].emplace(id.rawId(),MTDHit()).first;

    // --- Accumulate the energy (in MeV) of SIM hits in the same detector cell
    (simHitIt->second).energy += 1000.*simHit.energyLoss();

    // --- Get the time of the first SIM hit in the cell
    if( (simHitIt->second).time==0 || simHit.tof()<(simHitIt->second).time ) {

      (simHitIt->second).time = simHit.tof();

      auto hit_pos = simHit.entryPoint();
      (simHitIt->second).x = hit_pos.x();
      (simHitIt->second).y = hit_pos.y();
      (simHitIt->second).z = hit_pos.z();

    }

  } // simHit loop


  // ==============================================================================
  //  Histogram filling
  // ==============================================================================

  for (int idet=0; idet<2; ++idet){

    meNhits_[idet]->Fill(m_etlHits[idet].size());

    for (auto const& hit: m_etlTrkPerCell[idet])
      meNtrkPerCell_[idet]->Fill((hit.second).size());


    for (auto const& hit: m_etlHits[idet]) {

      if ( (hit.second).energy < hitMinEnergy_ ) continue;

      // --- Get the SIM hit global position
      ETLDetId detId(hit.first); 

      DetId geoId = detId.geographicalId();
      const MTDGeomDet* thedet = geom_->idToDet(geoId);
      if( thedet == nullptr )
	throw cms::Exception("EtlSimHitsValidation") << "GeographicalID: " << std::hex << geoId.rawId()
						     << " (" << detId.rawId()<< ") is invalid!" << std::dec
						     << std::endl;

      Local3DPoint local_point(0.1*(hit.second).x,0.1*(hit.second).y,0.1*(hit.second).z);
      const auto& global_point = thedet->toGlobal(local_point);

      // --- Fill the histograms

      meHitEnergy_[idet]->Fill((hit.second).energy);
      meHitTime_[idet]->Fill((hit.second).time);

      meHitXlocal_[idet]->Fill((hit.second).x);
      meHitYlocal_[idet]->Fill((hit.second).y);
      meHitZlocal_[idet]->Fill((hit.second).z);

      meOccupancy_[idet]->Fill(global_point.x(),global_point.y());

      meHitX_[idet]->Fill(global_point.x()); 
      meHitY_[idet]->Fill(global_point.y()); 
      meHitZ_[idet]->Fill(global_point.z()); 
      meHitPhi_[idet]->Fill(global_point.phi()); 
      meHitEta_[idet]->Fill(global_point.eta()); 

      meHitTvsE_[idet]->Fill((hit.second).energy,(hit.second).time);  
      meHitEvsPhi_[idet]->Fill(global_point.phi(),(hit.second).energy);
      meHitEvsEta_[idet]->Fill(global_point.eta(),(hit.second).energy);
      meHitTvsPhi_[idet]->Fill(global_point.phi(),(hit.second).time);
      meHitTvsEta_[idet]->Fill(global_point.eta(),(hit.second).time);


    } // hit loop

  } // idet loop

}


// ------------ method for histogram booking ------------
void EtlSimHitsValidation::bookHistograms(DQMStore::IBooker & ibook,
                               edm::Run const& run,
                               edm::EventSetup const & iSetup) {

  ibook.setCurrentFolder(folder_);

  // --- histograms booking

  meNhits_[1]     = ibook.book1D("EtlNhitsZpos", "Number of ETL cells with SIM hits (+Z);N_{ETL cells}", 250, 0., 5000.);
  meNhits_[0]     = ibook.book1D("EtlNhitsZneg", "Number of ETL cells with SIM hits (-Z);N_{ETL cells}", 250, 0., 5000.);
  meNtrkPerCell_[1] = ibook.book1D("EtlNtrkPerCellZpos", "Number of tracks per ETL sensor (+Z);N_{trk}", 10, 0., 10.);
  meNtrkPerCell_[0] = ibook.book1D("EtlNtrkPerCellZneg", "Number of tracks per ETL sensor (-Z);N_{trk}", 10, 0., 10.);

  meHitEnergy_[1] = ibook.book1D("EtlHitEnergyZpos", "ETL SIM hits energy (+Z);E_{SIM} [MeV]", 200, 0., 2.);
  meHitEnergy_[0] = ibook.book1D("EtlHitEnergyZneg", "ETL SIM hits energy (-Z);E_{SIM} [MeV]", 200, 0., 2.);
  meHitTime_[1]   = ibook.book1D("EtlHitTimeZpos", "ETL SIM hits ToA (+Z);ToA_{SIM} [ns]", 250, 0., 25.);
  meHitTime_[0]   = ibook.book1D("EtlHitTimeZneg", "ETL SIM hits ToA (-Z);ToA_{SIM} [ns]", 250, 0., 25.);

  meHitXlocal_[1] = ibook.book1D("EtlHitXlocalZpos", "ETL SIM local X (+Z);X_{SIM}^{LOC} [mm]", 100, -25., 25.);
  meHitXlocal_[0] = ibook.book1D("EtlHitXlocalZneg", "ETL SIM local X (-Z);X_{SIM}^{LOC} [mm]", 100, -25., 25.);
  meHitYlocal_[1] = ibook.book1D("EtlHitYlocalZpos", "ETL SIM local Y (+Z);Y_{SIM}^{LOC} [mm]", 200, -50., 50.);
  meHitYlocal_[0] = ibook.book1D("EtlHitYlocalZneg", "ETL SIM local Y (-Z);Y_{SIM}^{LOC} [mm]", 200, -50., 50.);
  meHitZlocal_[1] = ibook.book1D("EtlHitZlocalZpos", "ETL SIM local Z (+Z);Z_{SIM}^{LOC} [mm]", 80, -0.2, 0.2);
  meHitZlocal_[0] = ibook.book1D("EtlHitZlocalZneg", "ETL SIM local Z (-Z);Z_{SIM}^{LOC} [mm]", 80, -0.2, 0.2);

  meOccupancy_[1] = ibook.book2D("EtlOccupancyZpos", "ETL SIM hits occupancy (+Z);X_{SIM} [cm];Y_{SIM} [cm]",
				  135, -135., 135.,  135, -135., 135.);
  meOccupancy_[0] = ibook.book2D("EtlOccupancyZneg", "ETL SIM hits occupancy (-Z);X_{SIM} [cm];Y_{SIM} [cm]",
				  135, -135., 135.,  135, -135., 135.);

  meHitX_[1]      = ibook.book1D("EtlHitXZpos", "ETL SIM hits X (+Z);X_{SIM} [cm]", 135, -135., 135.);
  meHitX_[0]      = ibook.book1D("EtlHitXZneg", "ETL SIM hits X (-Z);X_{SIM} [cm]", 135, -135., 135.);
  meHitY_[1]      = ibook.book1D("EtlHitYZpos", "ETL SIM hits Y (+Z);Y_{SIM} [cm]", 135, -135., 135.);
  meHitY_[0]      = ibook.book1D("EtlHitYZneg", "ETL SIM hits Y (-Z);Y_{SIM} [cm]", 135, -135., 135.);
  meHitZ_[1]      = ibook.book1D("EtlHitZZpos", "ETL SIM hits Z (+Z);Z_{SIM} [cm]", 100,  303.,  304.5);
  meHitZ_[0]      = ibook.book1D("EtlHitZZneg", "ETL SIM hits Z (-Z);Z_{SIM} [cm]", 100, -304.5, -303.);

  meHitPhi_[1]    = ibook.book1D("EtlHitPhiZpos", "ETL SIM hits #phi (+Z);#phi_{SIM} [rad]", 315, -3.15, 3.15);
  meHitPhi_[0]    = ibook.book1D("EtlHitPhiZneg", "ETL SIM hits #phi (-Z);#phi_{SIM} [rad]", 315, -3.15, 3.15);
  meHitEta_[1]    = ibook.book1D("EtlHitEtaZpos", "ETL SIM hits #eta (+Z);#eta_{SIM}", 200,  1.55,  3.05);
  meHitEta_[0]    = ibook.book1D("EtlHitEtaZneg", "ETL SIM hits #eta (-Z);#eta_{SIM}", 200, -3.05, -1.55);

  meHitTvsE_[1]    = ibook.bookProfile("EtlHitTvsEZpos", "ETL SIM time vs energy (+Z);E_{SIM} [MeV];T_{SIM} [ns]",
				       100, 0., 2., 0., 100.);
  meHitTvsE_[0]    = ibook.bookProfile("EtlHitTvsEZneg", "ETL SIM time vs energy (-Z);E_{SIM} [MeV];T_{SIM} [ns]",
				       100, 0., 2., 0., 100.);
  meHitEvsPhi_[1]  = ibook.bookProfile("EtlHitEvsPhiZpos", "ETL SIM energy vs #phi (+Z);#phi_{SIM} [rad];E_{SIM} [MeV]",
				       100, -3.15, 3.15, 0., 100.);
  meHitEvsPhi_[0]  = ibook.bookProfile("EtlHitEvsPhiZneg", "ETL SIM energy vs #phi (-Z);#phi_{SIM} [rad];E_{SIM} [MeV]",
				       100, -3.15, 3.15, 0., 100.);
  meHitEvsEta_[1]  = ibook.bookProfile("EtlHitEvsEtaZpos","ETL SIM energy vs #eta (+Z);#eta_{SIM};E_{SIM} [MeV]",
				       200, 1.55, 3.05, 0., 100.);
  meHitEvsEta_[0]  = ibook.bookProfile("EtlHitEvsEtaZneg","ETL SIM energy vs #eta (-Z);#eta_{SIM};E_{SIM} [MeV]",
				       200, -3.05, -1.55, 0., 100.);
  meHitTvsPhi_[1]  = ibook.bookProfile("EtlHitTvsPhiZpos", "ETL SIM time vs #phi (+Z);#phi_{SIM} [rad];T_{SIM} [ns]",
				       100, -3.15, 3.15, 0., 100.);
  meHitTvsPhi_[0]  = ibook.bookProfile("EtlHitTvsPhiZneg", "ETL SIM time vs #phi (-Z);#phi_{SIM} [rad];T_{SIM} [ns]",
				       100, -3.15, 3.15, 0., 100.);
  meHitTvsEta_[1] = ibook.bookProfile("EtlHitTvsEtaZpos","ETL SIM time vs #eta (+Z);#eta_{SIM};T_{SIM} [ns]",
				       200, 1.55, 3.05, 0., 100.);
  meHitTvsEta_[0] = ibook.bookProfile("EtlHitTvsEtaZpos","ETL SIM time vs #eta (-Z);#eta_{SIM};T_{SIM} [ns]",
				       200, -3.05, -1.55, 0., 100.);

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EtlSimHitsValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/ETL/SimHits");
  desc.add<std::string>("moduleLabelG4","g4SimHits");
  desc.add<std::string>("etlSimHitsCollection","FastTimerHitsEndcap");
  desc.add<double>("hitMinimumEnergy",0.1); // [MeV]

  descriptions.add("etlSimHits", desc);

}

DEFINE_FWK_MODULE(EtlSimHitsValidation);
