// -*- C++ -*-
//
// Package:    Validation/MtdValidation
// Class:      EtlRecHitsValidation
//
/**\class EtlRecHitsValidation EtlRecHitsValidation.cc Validation/MtdValidation/plugins/EtlRecHitsValidation.cc

 Description: ETL RECO hits validation

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

#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"


class EtlRecHitsValidation : public DQMEDAnalyzer {

public:
  explicit EtlRecHitsValidation(const edm::ParameterSet&);
  ~EtlRecHitsValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  void bookHistograms(DQMStore::IBooker &,
		      edm::Run const&,
		      edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ------------ member data ------------

  const MTDGeometry* geom_;

  const std::string folder_;
  const std::string infoLabel_;
  const std::string etlRecHitsCollection_;

  int eventCount_;

  edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;

  // --- histograms declaration

  MonitorElement* meNhits_[2];

  MonitorElement* meHitEnergy_[2];
  MonitorElement* meHitTime_[2];

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
EtlRecHitsValidation::EtlRecHitsValidation(const edm::ParameterSet& iConfig):
  geom_(nullptr),
  folder_(iConfig.getParameter<std::string>("folder")),
  infoLabel_(iConfig.getParameter<std::string>("moduleLabel")),
  etlRecHitsCollection_(iConfig.getParameter<std::string>("etlRecHitsCollection")),
  eventCount_(0) {

  etlRecHitsToken_ = consumes <FTLRecHitCollection> (edm::InputTag(std::string(infoLabel_),
								   std::string(etlRecHitsCollection_)));

}

EtlRecHitsValidation::~EtlRecHitsValidation() {
}


// ------------ method called for each event  ------------
void EtlRecHitsValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;

  edm::LogInfo("EventInfo") << " Run = " << iEvent.id().run() << " Event = " << iEvent.id().event();

  edm::ESHandle<MTDGeometry> geometryHandle;
  iSetup.get<MTDDigiGeometryRecord>().get(geometryHandle);
  geom_ = geometryHandle.product();

  edm::Handle<FTLRecHitCollection> etlRecHitsHandle;
  iEvent.getByToken(etlRecHitsToken_, etlRecHitsHandle);

  if( ! etlRecHitsHandle.isValid() ) {
    edm::LogWarning("DataNotFound") << "No ETL RecHits found";
    return;
  }

  eventCount_++;

  // --- Loop over the ELT RECO hits

  unsigned int n_reco_etl[2] = {0,0};

  for (const auto& recHit: *etlRecHitsHandle) {

    ETLDetId detId = recHit.id();

    DetId geoId = detId.geographicalId();
    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    if( thedet == nullptr )
      throw cms::Exception("EtlRecHitsValidation") << "GeographicalID: " << std::hex << geoId.rawId()
						   << " (" << detId.rawId()<< ") is invalid!" << std::dec
						   << std::endl;

    const auto& global_point = thedet->toGlobal(Local3DPoint(0.,0.,0.));

    // --- Fill the histograms

    int idet = (detId.zside()+1)/2;

    meHitEnergy_[idet]->Fill(recHit.energy());
    meHitTime_[idet]->Fill(recHit.time());

    meOccupancy_[idet]->Fill(global_point.x(),global_point.y());

    meHitX_[idet]->Fill(global_point.x());
    meHitY_[idet]->Fill(global_point.y());
    meHitZ_[idet]->Fill(global_point.z());
    meHitPhi_[idet]->Fill(global_point.phi());
    meHitEta_[idet]->Fill(global_point.eta());

    meHitTvsE_[idet]->Fill(recHit.energy(),recHit.time());
    meHitEvsPhi_[idet]->Fill(global_point.phi(),recHit.energy());
    meHitEvsEta_[idet]->Fill(global_point.eta(),recHit.energy());
    meHitTvsPhi_[idet]->Fill(global_point.phi(),recHit.time());
    meHitTvsEta_[idet]->Fill(global_point.eta(),recHit.time());

    n_reco_etl[idet]++;

  } // recHit loop

  meNhits_[0]->Fill(n_reco_etl[0]);
  meNhits_[1]->Fill(n_reco_etl[1]);

}


// ------------ method for histogram booking ------------
void EtlRecHitsValidation::bookHistograms(DQMStore::IBooker & ibook,
                               edm::Run const& run,
                               edm::EventSetup const & iSetup) {

  ibook.setCurrentFolder(folder_);

  // --- histograms booking

  meNhits_[0]     = ibook.book1D("EtlNhitsZneg", "Number of ETL RECO hits (-Z);N_{RECO}", 250, 0., 5000.);
  meNhits_[1]     = ibook.book1D("EtlNhitsZpos", "Number of ETL RECO hits (+Z);N_{RECO}", 250, 0., 5000.);

  meHitEnergy_[0] = ibook.book1D("EtlHitEnergyZneg", "ETL RECO hits energy (-Z);E_{RECO} [MeV]", 150, 0., 3.);
  meHitEnergy_[1] = ibook.book1D("EtlHitEnergyZpos", "ETL RECO hits energy (+Z);E_{RECO} [MeV]", 150, 0., 3.);

  meHitTime_[0] = ibook.book1D("EtlHitTimeZneg", "ETL RECO hits ToA (-Z);ToA_{RECO} [ns]", 250, 0., 25.);
  meHitTime_[1] = ibook.book1D("EtlHitTimeZpos", "ETL RECO hits ToA (+Z);ToA_{RECO} [ns]", 250, 0., 25.);

  meOccupancy_[0] = ibook.book2D("EtlOccupancyZneg","ETL RECO hits occupancy (-Z);X_{RECO} [cm];Y_{RECO} [cm]",
				 59, -130., 130.,  59, -130., 130.);
  meOccupancy_[1] = ibook.book2D("EtlOccupancyZpos","ETL DIGI hits occupancy (+Z);X_{RECO} [cm];Y_{RECO} [cm]",
				 59, -130., 130.,  59, -130., 130.);

  meHitX_[1]      = ibook.book1D("EtlHitXZpos", "ETL RECO hits X (+Z);X_{RECO} [cm]", 135, -135., 135.);
  meHitX_[0]      = ibook.book1D("EtlHitXZneg", "ETL RECO hits X (-Z);X_{RECO} [cm]", 135, -135., 135.);
  meHitY_[1]      = ibook.book1D("EtlHitYZpos", "ETL RECO hits Y (+Z);Y_{RECO} [cm]", 135, -135., 135.);
  meHitY_[0]      = ibook.book1D("EtlHitYZneg", "ETL RECO hits Y (-Z);Y_{RECO} [cm]", 135, -135., 135.);
  meHitZ_[1]      = ibook.book1D("EtlHitZZpos", "ETL RECO hits Z (+Z);Z_{RECO} [cm]", 100,  303.,  304.5);
  meHitZ_[0]      = ibook.book1D("EtlHitZZneg", "ETL RECO hits Z (-Z);Z_{RECO} [cm]", 100, -304.5, -303.);

  meHitPhi_[1]    = ibook.book1D("EtlHitPhiZpos", "ETL RECO hits #phi (+Z);#phi_{RECO} [rad]", 315, -3.15, 3.15);
  meHitPhi_[0]    = ibook.book1D("EtlHitPhiZneg", "ETL RECO hits #phi (-Z);#phi_{RECO} [rad]", 315, -3.15, 3.15);
  meHitEta_[1]    = ibook.book1D("EtlHitEtaZpos", "ETL RECO hits #eta (+Z);#eta_{RECO}", 200,  1.55,  3.05);
  meHitEta_[0]    = ibook.book1D("EtlHitEtaZneg", "ETL RECO hits #eta (-Z);#eta_{RECO}", 200, -3.05, -1.55);

  meHitTvsE_[1]    = ibook.bookProfile("EtlHitTvsEZpos", "ETL RECO time vs energy (+Z);E_{RECO} [MeV];ToA_{RECO} [ns]",
				       100, 0., 2., 0., 100.);
  meHitTvsE_[0]    = ibook.bookProfile("EtlHitTvsEZneg", "ETL RECO time vs energy (-Z);E_{RECO} [MeV];ToA_{RECO} [ns]",
				       100, 0., 2., 0., 100.);
  meHitEvsPhi_[1]  = ibook.bookProfile("EtlHitEvsPhiZpos", "ETL RECO energy vs #phi (+Z);#phi_{RECO} [rad];E_{RECO} [MeV]",
				       100, -3.15, 3.15, 0., 100.);
  meHitEvsPhi_[0]  = ibook.bookProfile("EtlHitEvsPhiZneg", "ETL RECO energy vs #phi (-Z);#phi_{RECO} [rad];E_{RECO} [MeV]",
				       100, -3.15, 3.15, 0., 100.);
  meHitEvsEta_[1]  = ibook.bookProfile("EtlHitEvsEtaZpos","ETL RECO energy vs #eta (+Z);#eta_{RECO};E_{RECO} [MeV]",
				       200, 1.55, 3.05, 0., 100.);
  meHitEvsEta_[0]  = ibook.bookProfile("EtlHitEvsEtaZneg","ETL RECO energy vs #eta (-Z);#eta_{RECO};E_{RECO} [MeV]",
				       200, -3.05, -1.55, 0., 100.);
  meHitTvsPhi_[1]  = ibook.bookProfile("EtlHitTvsPhiZpos", "ETL RECO time vs #phi (+Z);#phi_{RECO} [rad];ToA_{RECO} [ns]",
				       100, -3.15, 3.15, 0., 100.);
  meHitTvsPhi_[0]  = ibook.bookProfile("EtlHitTvsPhiZneg", "ETL RECO time vs #phi (-Z);#phi_{RECO} [rad];ToA_{RECO} [ns]",
				       100, -3.15, 3.15, 0., 100.);
  meHitTvsEta_[1] = ibook.bookProfile("EtlHitTvsEtaZpos","ETL RECO time vs #eta (+Z);#eta_{RECO};ToA_{RECO} [ns]",
				       200, 1.55, 3.05, 0., 100.);
  meHitTvsEta_[0] = ibook.bookProfile("EtlHitTvsEtaZpos","ETL RECO time vs #eta (-Z);#eta_{RECO};ToA_{RECO} [ns]",
				       200, -3.05, -1.55, 0., 100.);

}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EtlRecHitsValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/ETL/RecHits");
  desc.add<std::string>("moduleLabel","mtdRecHits");
  desc.add<std::string>("etlRecHitsCollection","FTLEndcap");

  descriptions.add("etlRecHits", desc);

}

DEFINE_FWK_MODULE(EtlRecHitsValidation);
