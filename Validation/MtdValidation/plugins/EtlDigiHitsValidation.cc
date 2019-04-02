// -*- C++ -*-
//
// Package:    Validation/MtdValidation
// Class:      EtlDigiHitsValidation
//
/**\class EtlDigiHitsValidation EtlDigiHitsValidation.cc Validation/MtdValidation/plugins/EtlDigiHitsValidation.cc

 Description: ETL DIGI hits validation

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

#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"

#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"


class EtlDigiHitsValidation : public DQMEDAnalyzer {

public:
  explicit EtlDigiHitsValidation(const edm::ParameterSet&);
  ~EtlDigiHitsValidation() override;

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
  const std::string etlDigisCollection_;

  int eventCount_;

  edm::EDGetTokenT<ETLDigiCollection> etlDigiHitsToken_;

  // --- histograms declaration

  MonitorElement* meNhits_[2];

  MonitorElement* meHitCharge_[2];
  MonitorElement* meHitTime_[2];

  MonitorElement* meOccupancy_[2];

  MonitorElement* meHitX_[2];
  MonitorElement* meHitY_[2];
  MonitorElement* meHitZ_[2];
  MonitorElement* meHitPhi_[2];
  MonitorElement* meHitEta_[2];

  MonitorElement* meHitTvsQ_[2];
  MonitorElement* meHitQvsPhi_[2];
  MonitorElement* meHitQvsEta_[2];
  MonitorElement* meHitTvsPhi_[2];
  MonitorElement* meHitTvsEta_[2];

};


// ------------ constructor and destructor --------------
EtlDigiHitsValidation::EtlDigiHitsValidation(const edm::ParameterSet& iConfig):
  geom_(nullptr),
  folder_(iConfig.getParameter<std::string>("folder")),
  infoLabel_(iConfig.getParameter<std::string>("moduleLabel")),
  etlDigisCollection_(iConfig.getParameter<std::string>("etlDigiHitsCollection")),
  eventCount_(0) {

  etlDigiHitsToken_ = consumes <ETLDigiCollection> (edm::InputTag(std::string(infoLabel_),
								  std::string(etlDigisCollection_)));

}

EtlDigiHitsValidation::~EtlDigiHitsValidation() {
}


// ------------ method called for each event  ------------
void EtlDigiHitsValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;

  edm::LogInfo("EventInfo") << " Run = " << iEvent.id().run() << " Event = " << iEvent.id().event();

  edm::ESHandle<MTDGeometry> geometryHandle;
  iSetup.get<MTDDigiGeometryRecord>().get(geometryHandle);
  geom_ = geometryHandle.product();

  edm::Handle<ETLDigiCollection> etlDigiHitsHandle;
  iEvent.getByToken(etlDigiHitsToken_, etlDigiHitsHandle);

  if( ! etlDigiHitsHandle.isValid() ) {
    edm::LogWarning("DataNotFound") << "No ETL DIGI hits found";
    return;
  }

  eventCount_++;
  
  // --- Loop over the ELT DIGI hits

  unsigned int n_digi_etl[2] = {0,0};

  for (const auto& dataFrame: *etlDigiHitsHandle) {

    // --- Get the on-time sample
    int isample = 2;

    const auto& sample = dataFrame.sample(isample);

    ETLDetId detId = dataFrame.id();

    DetId geoId = detId.geographicalId();

    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    if( thedet == nullptr )
      throw cms::Exception("EtlDigiHitsValidation") << "GeographicalID: " << std::hex << geoId.rawId()
						    << " (" << detId.rawId()<< ") is invalid!" << std::dec
						    << std::endl;
    
    Local3DPoint local_point(0.,0.,0.);
    const auto& global_point = thedet->toGlobal(local_point);

    // --- Fill the histograms

    int idet = (detId.zside()+1)/2;

    meHitCharge_[idet]->Fill(sample.data());
    meHitTime_[idet]->Fill(sample.toa());
    meOccupancy_[idet]->Fill(global_point.x(),global_point.y());

    
    meHitX_[idet]->Fill(global_point.x());
    meHitY_[idet]->Fill(global_point.y());
    meHitZ_[idet]->Fill(global_point.z());
    meHitPhi_[idet]->Fill(global_point.phi());
    meHitEta_[idet]->Fill(global_point.eta());

    meHitTvsQ_[idet]->Fill(sample.data(),sample.toa());
    meHitQvsPhi_[idet]->Fill(global_point.phi(),sample.data());
    meHitQvsEta_[idet]->Fill(global_point.eta(),sample.data());
    meHitTvsPhi_[idet]->Fill(global_point.phi(),sample.toa());
    meHitTvsEta_[idet]->Fill(global_point.eta(),sample.toa());

    n_digi_etl[idet]++;

  } // dataFrame loop

  meNhits_[0]->Fill(n_digi_etl[0]);
  meNhits_[1]->Fill(n_digi_etl[1]);

}


// ------------ method for histogram booking ------------
void EtlDigiHitsValidation::bookHistograms(DQMStore::IBooker & ibook,
                               edm::Run const& run,
                               edm::EventSetup const & iSetup) {

  ibook.setCurrentFolder(folder_);

  // --- histograms booking

  meNhits_[0]     = ibook.book1D("EtlNhitsZneg", "Number of ETL DIGI hits (-Z);N_{DIGI}", 250, 0., 5000.);
  meNhits_[1]     = ibook.book1D("EtlNhitsZpos", "Number of ETL DIGI hits (+Z);N_{DIGI}", 250, 0., 5000.);

  meHitCharge_[0] = ibook.book1D("EtlHitChargeZneg", "ETL DIGI hits charge (-Z);Q_{DIGI} [ADC counts]",
				 256, 0., 256.);
  meHitCharge_[1] = ibook.book1D("EtlHitChargeZpos", "ETL DIGI hits charge (+Z);Q_{DIGI} [ADC counts]",
				 256, 0., 256.);
  meHitTime_[0]   = ibook.book1D("EtlHitTimeZneg", "ETL DIGI hits ToA (-Z);ToA_{DIGI} [TDC counts]", 
				 1000, 0., 2000.);
  meHitTime_[1]   = ibook.book1D("EtlHitTimeZpos", "ETL DIGI hits ToA (+Z);ToA_{DIGI} [TDC counts]", 
				 1000, 0., 2000.);

  meOccupancy_[0] = ibook.book2D("EtlOccupancyZneg","ETL DIGI hits occupancy (-Z);X_{DIGI} [cm];Y_{DIGI} [cm]",
				 135, -135., 135.,  135, -135., 135.);
  meOccupancy_[1] = ibook.book2D("EtlOccupancyZpos","ETL DIGI hits occupancy (+Z);X_{DIGI} [cm];Y_{DIGI} [cm]",
				 135, -135., 135.,  135, -135., 135.);

  meHitX_[0]      = ibook.book1D("EtlHitXZneg", "ETL DIGI hits X (-Z);X_{DIGI} [cm]", 135, -135., 135.);
  meHitX_[1]      = ibook.book1D("EtlHitXZpos", "ETL DIGI hits X (+Z);X_{DIGI} [cm]", 135, -135., 135.);
  meHitY_[0]      = ibook.book1D("EtlHitYZneg", "ETL DIGI hits Y (-Z);Y_{DIGI} [cm]", 135, -135., 135.);
  meHitY_[1]      = ibook.book1D("EtlHitYZpos", "ETL DIGI hits Y (+Z);Y_{DIGI} [cm]", 135, -135., 135.);
  meHitZ_[0]      = ibook.book1D("EtlHitZZneg", "ETL DIGI hits Z (-Z);Z_{DIGI} [cm]", 100, -304.5, -303.);
  meHitZ_[1]      = ibook.book1D("EtlHitZZpos", "ETL DIGI hits Z (+Z);Z_{DIGI} [cm]", 100,  303.,  304.5);

  meHitPhi_[0]    = ibook.book1D("EtlHitPhiZneg", "ETL DIGI hits #phi (-Z);#phi_{DIGI} [rad]", 315, -3.15, 3.15);
  meHitPhi_[1]    = ibook.book1D("EtlHitPhiZpos", "ETL DIGI hits #phi (+Z);#phi_{DIGI} [rad]", 315, -3.15, 3.15);
  meHitEta_[0]    = ibook.book1D("EtlHitEtaZneg", "ETL DIGI hits #eta (-Z);#eta_{DIGI}", 200,  1.55,  3.05);
  meHitEta_[1]    = ibook.book1D("EtlHitEtaZpos", "ETL DIGI hits #eta (+Z);#eta_{DIGI}", 200, -3.05, -1.55);


  meHitTvsQ_[0]   = ibook.bookProfile("EtlHitTvsQZneg", "ETL DIGI ToA vs charge (-Z);Q_{DIGI} [ADC counts];ToA_{DIGI} [TDC counts]",
				       256, 0., 256., 0., 1024.);
  meHitTvsQ_[1]   = ibook.bookProfile("EtlHitTvsQZpos", "ETL DIGI ToA vs charge (+Z);Q_{DIGI} [ADC counts];ToA_{DIGI} [TDC counts]",
				       256, 0., 256., 0., 1024.);
  meHitQvsPhi_[0] = ibook.bookProfile("EtlHitQvsPhiZneg", "ETL DIGI charge vs #phi (-Z);#phi_{DIGI} [rad];Q_{DIGI} [ADC counts]",
				      100, -3.15, 3.15, 0., 1024.);
  meHitQvsPhi_[1] = ibook.bookProfile("EtlHitQvsPhiZpos", "ETL DIGI charge vs #phi (+Z);#phi_{DIGI} [rad];Q_{DIGI} [ADC counts]",
				      100, -3.15, 3.15, 0., 1024.);
  meHitQvsEta_[0] = ibook.bookProfile("EtlHitQvsEtaZneg","ETL DIGI charge vs #eta (-Z);#eta_{DIGI};Q_{DIGI} [ADC counts]",
				      100, -3.05, -1.55, 0., 1024.);
  meHitQvsEta_[1] = ibook.bookProfile("EtlHitQvsEtaZpos","ETL DIGI charge vs #eta (+Z);#eta_{DIGI};Q_{DIGI} [ADC counts]",
				      100, 1.55,  3.05, 0., 1024.);
  meHitTvsPhi_[0] = ibook.bookProfile("EtlHitTvsPhiZneg", "ETL DIGI ToA vs #phi (-Z);#phi_{DIGI} [rad];ToA_{DIGI} [TDC counts]",
				      100, -3.15, 3.15, 0., 1024.);
  meHitTvsPhi_[1] = ibook.bookProfile("EtlHitTvsPhiZpos", "ETL DIGI ToA vs #phi (+Z);#phi_{DIGI} [rad];ToA_{DIGI} [TDC counts]",
				      100, -3.15, 3.15, 0., 1024.);
  meHitTvsEta_[0] = ibook.bookProfile("EtlHitTvsEtaZneg","ETL DIGI ToA vs #eta (-Z);#eta_{DIGI};ToA_{DIGI} [TDC counts]",
				      100, -3.05, -1.55, 0., 1024.);
  meHitTvsEta_[1] = ibook.bookProfile("EtlHitTvsEtaZpos","ETL DIGI ToA vs #eta (+Z);#eta_{DIGI};ToA_{DIGI} [TDC counts]",
				      100, 1.55,  3.05, 0., 1024.);

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EtlDigiHitsValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/ETL/DigiHits");
  desc.add<std::string>("moduleLabel","mix");
  desc.add<std::string>("etlDigiHitsCollection","FTLEndcap");

  descriptions.add("etlDigiHits", desc);

}

DEFINE_FWK_MODULE(EtlDigiHitsValidation);
