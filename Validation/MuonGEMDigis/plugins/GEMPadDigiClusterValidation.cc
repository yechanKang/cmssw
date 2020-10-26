#include "Validation/MuonGEMDigis/plugins/GEMPadDigiClusterValidation.h"
#include <TMath.h>

GEMPadDigiClusterValidation::GEMPadDigiClusterValidation(const edm::ParameterSet& pset)
    : GEMBaseValidation(pset, "GEMPadDigiClusterValidation") {
  const auto& pad_cluster_pset = pset.getParameterSet("gemPadCluster");
  const auto& pad_cluster_tag = pad_cluster_pset.getParameter<edm::InputTag>("inputTag");
  pad_cluster_token_ = consumes<GEMPadDigiClusterCollection>(pad_cluster_tag);

  const auto& simhit_pset = pset.getParameterSet("gemSimHit");
  const auto& simhit_tag = simhit_pset.getParameter<edm::InputTag>("inputTag");
  simhit_token_ = consumes<edm::PSimHitContainer>(simhit_tag);
  
  geomToken_ = esConsumes<GEMGeometry, MuonGeometryRecord>();
  geomTokenBeginRun_ = esConsumes<GEMGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

void GEMPadDigiClusterValidation::bookHistograms(DQMStore::IBooker& booker,
                                                 edm::Run const& Run,
                                                 edm::EventSetup const& setup) {
  const GEMGeometry* gem = &setup.getData(geomTokenBeginRun_);

  // NOTE Occupancy
  booker.setCurrentFolder("MuonGEMDigisV/GEMDigisTask/PadCluster/Occupancy");
  for (const auto& region : gem->regions()) {
    Int_t region_id = region->region();

    me_occ_zr_.emplace(region_id, bookZROccupancy(booker, region_id, "pad", "Pad Cluster"));

    for (const auto& station : region->stations()) {
      Int_t station_id = station->station();
      ME2IdsKey key2{region_id, station_id};

      me_occ_det_[key2] = bookDetectorOccupancy(booker, key2, station, "pad_cluster", "Pad Cluster");
      me_occ_pad_det_[key2] = bookDetectorOccupancy(booker, key2, station, "pad", "Pad from Pad Cluster");
      me_occ_simhit_det_[key2] = bookDetectorOccupancy(booker, key2, station, "simhit", "simhit");
      me_occ_matched_det_[key2] = bookDetectorOccupancy(booker, key2, station, "matched", "mathced with simhit");

      const auto& superChamberVec = station->superChambers();
      if (superChamberVec.empty()) {
        edm::LogError(kLogCategory_) << "Super chambers missing for region = " << region_id
                                     << " and station = " << station_id;
        continue;
      }
      const GEMSuperChamber* super_chamber = superChamberVec.front();
      if (super_chamber == nullptr) {
        edm::LogError(kLogCategory_) << "Failed to find super chamber for region = " << region_id
                                     << " and station = " << station_id;
        continue;
      }
      for (const auto& superChamber : station->superChambers()) {
        for (const auto& chamber : superChamber->chambers()) {
          Int_t layer_id = chamber->id().layer();
          Int_t chamber_id = chamber->id().chamber();
          ME4IdsKey key4{region_id, station_id, chamber_id, layer_id};
          if (detail_plot_) {
            me_detail_pad_cluster_count_[key4] =
                bookHist1D(booker, key4, "pad_cluster_count", "Pad Cluster count", 50, 0, 50, "Pad Cluster Count");
          }
        }
      }
      for (const auto& chamber : super_chamber->chambers()) {
        Int_t layer_id = chamber->id().layer();
        ME3IdsKey key3{region_id, station_id, layer_id};

        const auto& etaPartitionVec = chamber->etaPartitions();
        if (etaPartitionVec.empty() || etaPartitionVec.front() == nullptr) {
          edm::LogError(kLogCategory_)
              << "Eta partition missing or null for region, station, super chamber, chamber = (" << region_id << ", "
              << station_id << ", " << super_chamber->id() << ", " << chamber->id() << ")";
          continue;
        }
        Int_t num_pads = etaPartitionVec.front()->npads();

        if (detail_plot_) {

          me_detail_occ_xy_[key3] = bookXYOccupancy(booker, key3, "pad", "Pad Cluster");

          me_detail_occ_phi_pad_[key3] = bookHist2D(booker,
                                                    key3,
                                                    "occ_phi_pad",
                                                    "Pad Cluster Occupancy",
                                                    280,
                                                    -M_PI,
                                                    M_PI,
                                                    num_pads / 2,
                                                    0,
                                                    num_pads,
                                                    "#phi [rad]",
                                                    "Pad number");

          me_detail_occ_pad_[key3] =
              bookHist1D(booker, key3, "occ_pad", "Pad Cluster Occupancy", num_pads, 0.5, num_pads + 0.5, "Pad number");
          me_detail_pad_cluster_size_[key3] =
              bookHist1D(booker, key3, "pad_cluster_size", "Pad Cluster Size", 100, 0.5, 100 + 0.5, "Pad Cluster Size");
        }
      }  // end loop over layer ids
    }    // end loop over station ids
  }      // end loop over region ids

  // NOTE Bunch Crossing
  if (detail_plot_) {
    booker.setCurrentFolder("MuonGEMDigisV/GEMDigisTask/PadCluster/BunchCrossing");

    for (const auto& region : gem->regions()) {
      Int_t region_id = region->region();
      for (const auto& station : region->stations()) {
        Int_t station_id = station->station();

        const auto& superChamberVec = station->superChambers();
        if (superChamberVec.empty()) {
          edm::LogError(kLogCategory_) << "Super chambers missing for region = " << region_id
                                       << " and station = " << station_id;
          continue;
        }
        const GEMSuperChamber* super_chamber = superChamberVec.front();
        if (super_chamber == nullptr) {
          edm::LogError(kLogCategory_) << "Failed to find super chamber for region = " << region_id
                                       << " and station = " << station_id;
          continue;
        }
        for (const auto& chamber : super_chamber->chambers()) {
          Int_t layer_id = chamber->id().layer();
          ME3IdsKey key3(region_id, station_id, layer_id);
          me_detail_bx_[key3] =
              bookHist1D(booker, key3, "bx", "Pad Cluster Bunch Crossing", 5, -2.5, 2.5, "Bunch crossing");
        }  // chamber loop
      }    // station loop
    }      // region loop
  }        // detail plot
}

GEMPadDigiClusterValidation::~GEMPadDigiClusterValidation() {}

void GEMPadDigiClusterValidation::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  const GEMGeometry* gem = &setup.getData(geomToken_);

  edm::Handle<GEMPadDigiClusterCollection> collection;
  event.getByToken(pad_cluster_token_, collection);
  if (not collection.isValid()) {
    edm::LogError(kLogCategory_) << "Cannot get pads by label GEMPadToken.";
    return;
  }

  edm::Handle<edm::PSimHitContainer> simhit_container;
  event.getByToken(simhit_token_, simhit_container);
  if (not simhit_container.isValid()) {
    edm::LogError(kLogCategory_) << "Failed to get PSimHitContainer." << std::endl;
    return;
  }

  for (auto range_iter = collection->begin(); range_iter != collection->end(); range_iter++) {
    GEMDetId gemid = (*range_iter).first;
    const auto& range = (*range_iter).second;

    if (gem->idToDet(gemid) == nullptr) {
      edm::LogError(kLogCategory_) << "Getting DetId failed. Discard this gem pad hit. "
                                   << "Maybe it comes from unmatched geometry." << std::endl;
      continue;
    }

    const GEMEtaPartition* roll = gem->etaPartition(gemid);
    const BoundPlane& surface = roll->surface();

    Int_t region_id = gemid.region();
    Int_t station_id = gemid.station();
    Int_t layer_id = gemid.layer();
    Int_t chamber_id = gemid.chamber();
    Int_t roll_id = gemid.roll();

    ME2IdsKey key2(region_id, station_id);
    ME3IdsKey key3(region_id, station_id, layer_id);
    ME4IdsKey key4(region_id, station_id, chamber_id, layer_id);

    int count = 0;
    for (auto digi = range.first; digi != range.second; ++digi) {
      // ignore 16-partition GE2/1 pads
      if (gemid.isGE21() and digi->nPartitions() == GEMPadDigiCluster::GE21SplitStrip)
        continue;

      const auto& padsVec = digi->pads();
      if (padsVec.empty()) {
        edm::LogError(kLogCategory_) << "Pads missing for digi from GEM ID = " << gemid;
        continue;
      }
      Int_t pad = padsVec[0];

      // bunch crossing
      Int_t bx = digi->bx();

      const LocalPoint& local_pos = roll->centreOfPad(pad);
      const GlobalPoint& global_pos = surface.toGlobal(local_pos);

      Float_t g_r = global_pos.perp();
      Float_t g_phi = global_pos.phi();
      Float_t g_x = global_pos.x();
      Float_t g_y = global_pos.y();
      Float_t g_abs_z = std::fabs(global_pos.z());

      me_occ_zr_[region_id]->Fill(g_abs_z, g_r);

      Int_t bin_x = getDetOccBinX(chamber_id, layer_id);
      me_occ_det_[key2]->Fill(bin_x, roll_id);

      for (const auto& pad : digi->pads()) {
        me_occ_pad_det_[key2]->Fill(bin_x, roll_id);
      }

      if (detail_plot_) {
        me_detail_occ_xy_[key3]->Fill(g_x, g_y);
        me_detail_occ_phi_pad_[key3]->Fill(g_phi, pad);
        me_detail_occ_pad_[key3]->Fill(pad);
        me_detail_pad_cluster_size_[key3]->Fill((digi->pads()).size());
        me_detail_bx_[key3]->Fill(bx);
        count++;
      }  // detail_plot_
    }
    if (detail_plot_) me_detail_pad_cluster_count_[key4]->Fill(count);
  }  // end loop over range iters

  // NOTE
  for (const auto& simhit : *simhit_container.product()) {
    if (not isMuonSimHit(simhit))
      continue;

    if (gem->idToDet(simhit.detUnitId()) == nullptr) {
      edm::LogError(kLogCategory_) << "SimHit did not match with GEMGeometry." << std::endl;
      continue;
    }

    GEMDetId simhit_gemid(simhit.detUnitId());

    Int_t region_id = simhit_gemid.region();
    Int_t station_id = simhit_gemid.station();
    Int_t layer_id = simhit_gemid.layer();
    Int_t chamber_id = simhit_gemid.chamber();
    Int_t roll_id = simhit_gemid.roll();

    ME2IdsKey key2{region_id, station_id};
    //ME3IdsKey key3{region_id, station_id, layer_id};

    //const GEMEtaPartition* roll = gem->etaPartition(simhit_gemid);

    //const auto& simhit_local_pos = simhit.localPosition();
    //const auto& simhit_global_pos = roll->surface().toGlobal(simhit_local_pos);

    Int_t bin_x = getDetOccBinX(chamber_id, layer_id);
    me_occ_simhit_det_[key2]->Fill(bin_x, roll_id);

    auto range = collection->get(simhit_gemid);
    if (range.first != range.second) me_occ_matched_det_[key2]->Fill(bin_x, roll_id);
  }  // simhit_container loop
}
