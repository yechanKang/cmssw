#include "Validation/MuonGEMHits/plugins/GEMSimHitValidation.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "TDatabasePDG.h"

GEMSimHitValidation::GEMSimHitValidation(const edm::ParameterSet& pset)
    : GEMBaseValidation(pset, "GEMSimHitValidation") {
  const auto& simhit_pset = pset.getParameterSet("gemSimHit");
  const auto& simhit_tag = simhit_pset.getParameter<edm::InputTag>("inputTag");
  simhit_token_ = consumes<edm::PSimHitContainer>(simhit_tag);

  tof_range_ = pset.getUntrackedParameter<std::vector<Double_t> >("TOFRange");
  pids_ = pset.getUntrackedParameter<std::vector<Int_t> >("pids");

  geomToken_ = esConsumes<GEMGeometry, MuonGeometryRecord>();
  geomTokenBeginRun_ = esConsumes<GEMGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

GEMSimHitValidation::~GEMSimHitValidation() {}

void GEMSimHitValidation::bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& setup) {
  const GEMGeometry* gem = &setup.getData(geomTokenBeginRun_);

  // NOTE Time of flight
  booker.setCurrentFolder("MuonGEMHitsV/GEMHitsTask/TimeOfFlight");

  TString tof_xtitle = "Time of flight [ns]";
  TString tof_ytitle = "Entries";

  const auto& regionsVec = gem->regions();
  if (regionsVec.empty() || regionsVec[0] == nullptr) {
    edm::LogError(kLogCategory_) << "Regions missing or null.";
    return;
  } else {
    for (const auto& station : regionsVec[0]->stations()) {
      Int_t station_id = station->station();
      const auto [tof_min, tof_max] = getTOFRange(station_id);
      auto tof_name = TString::Format("tof_muon_GE%d1", station_id);
      auto tof_title = TString::Format("SimHit Time Of Flight (Muon only) : Station %d", station_id);

      me_tof_mu_st_[station_id] = booker.book1D(tof_name, tof_title, 20, tof_min, tof_max);
      if (detail_plot_) {
        for (const auto& pid : pids_) {
          auto tof_name = TString::Format("tof_pid_%d_GE%d1", pid, station_id);
          auto tof_title = TString::Format("SimHit Time Of Flight (PDG %d) : Station %d", pid, station_id);
          me_detail_tof_pid_[pid][station_id] = booker.book1D(tof_name, tof_title, 20, tof_min, tof_max);
        }
        auto tof_name = TString::Format("tof_pid_others_GE%d1", station_id);
        auto tof_title = TString::Format("SimHit Time Of Flight (PDG Others) : Station %d", station_id);
        me_detail_tof_pid_[0][station_id] = booker.book1D(tof_name, tof_title, 20, tof_min, tof_max);
      }
    }  // end for
  }    // end else

  for (const auto& region : gem->regions()) {
    Int_t region_id = region->region();

    for (const auto& station : region->stations()) {
      Int_t station_id = station->station();

      const auto [tof_min, tof_max] = getTOFRange(station_id);
      const auto& superChamberVec = station->superChambers();
      if (superChamberVec.empty() || superChamberVec.front() == nullptr) {
        edm::LogError(kLogCategory_) << "Super chambers missing or null for region = " << region_id
                                     << " and station = " << station_id;
      } else {
        const GEMSuperChamber* super_chamber = superChamberVec.front();

        for (const auto& chamber : super_chamber->chambers()) {
          Int_t layer_id = chamber->id().layer();
          ME3IdsKey key3{region_id, station_id, layer_id};

          me_tof_[key3] = bookHist1D(
              booker, key3, "tof", "Time of Flight of Muon SimHits", 20, tof_min, tof_max, tof_xtitle, tof_ytitle);

          me_tof_mu_[key3] = bookHist1D(
              booker, key3, "tof_muon", "SimHit TOF (Muon only)", 20, tof_min, tof_max, tof_xtitle, tof_ytitle);
        }  // chamber loop
      }    // end else
    }      // station loop
  }        // region loop

  // NOTE Energy Loss
  booker.setCurrentFolder("MuonGEMHitsV/GEMHitsTask/EnergyLoss");

  TString eloss_xtitle = "Energy loss [eV]";
  TString eloss_ytitle = "Entries / 0.5 keV";

  for (const auto& station : gem->regions()[0]->stations()) {
    Int_t station_id = station->station();

    auto eloss_name = TString::Format("eloss_muon_GE%d1", station_id);
    auto eloss_title = TString::Format("SimHit Energy Loss (Muon only) : Station %d", station_id);

    me_eloss_mu_[station_id] =
        booker.book1D(eloss_name, eloss_title + ";" + eloss_xtitle + ";" + eloss_ytitle, 20, 0.0, 10.0);
  }  // station loop

  if (detail_plot_) {
    for (const auto& pid : pids_) {
      auto eloss_name = TString::Format("eloss_pid_%d", pid);
      auto eloss_title = TString::Format("SimHit Energy Loss (PID %d)", pid);
      me_detail_eloss_pid_[pid] =
          booker.book1D(eloss_name, eloss_title + ";" + eloss_xtitle + ";" + eloss_ytitle, 20, 0.0, 10.0);
    }
    auto eloss_name = TString::Format("eloss_pid_others");
    auto eloss_title = TString::Format("SimHit Energy Loss (PID Others)");
    me_detail_eloss_pid_[0] =
        booker.book1D(eloss_name, eloss_title + ";" + eloss_xtitle + ";" + eloss_ytitle, 20, 0.0, 10.0);
  }

  if (detail_plot_) {
    for (const auto& region : gem->regions()) {
      Int_t region_id = region->region();
      for (const auto& station : region->stations()) {
        Int_t station_id = station->station();
        const auto& superChamberVec = station->superChambers();
        if (superChamberVec.empty() || superChamberVec.front() == nullptr) {
          edm::LogError(kLogCategory_) << "Super chambers missing or null for region = " << region_id
                                       << " and station = " << station_id;
        } else {
          for (const auto& chamber : superChamberVec.front()->chambers()) {
            Int_t layer_id = chamber->id().layer();
            ME3IdsKey key3{region_id, station_id, layer_id};

            me_detail_eloss_[key3] =
                bookHist1D(booker, key3, "eloss", "SimHit Energy Loss", 20, 0.0, 10.0, eloss_xtitle, eloss_ytitle);

            me_detail_eloss_mu_[key3] = bookHist1D(
                booker, key3, "eloss_muon", "SimHit Energy Loss (Muon Only)", 20, 0.0, 10.0, eloss_xtitle, eloss_ytitle);

          }  // chamber loop
        }    // end else
      }      // station loop
    }        // region loop
  }          // detail plot

  // NOTE Occupancy
  booker.setCurrentFolder("MuonGEMHitsV/GEMHitsTask/Occupancy");

  me_total_hits_ = booker.book1D("total_hits", "Total number of hits on GEM", 100, -0.5, 99.5);

  for (const auto& region : gem->regions()) {
    Int_t region_id = region->region();

    if (detail_plot_)
      me_detail_occ_zr_[region_id] = bookZROccupancy(booker, region_id, "simhit", "SimHit");

    for (const auto& station : region->stations()) {
      Int_t station_id = station->station();
      ME2IdsKey key2{region_id, station_id};

      if (detail_plot_) {
        me_detail_occ_det_[key2] = bookDetectorOccupancy(booker, key2, station, "simhit", "SimHit");
        me_detail_mu_occ_det_[key2] = bookDetectorOccupancy(booker, key2, station, "muon_simhit", "Muon SimHit");
      }

      const auto& superChamberVec = station->superChambers();
      if (superChamberVec.empty() || superChamberVec.front() == nullptr) {
        edm::LogError(kLogCategory_) << "Super chambers missing or null for region = " << region_id
                                     << " and station = " << station_id;
      } else {
        const GEMSuperChamber* super_chamber = superChamberVec.front();
        for (const auto& chamber : super_chamber->chambers()) {
          Int_t layer_id = chamber->id().layer();
          Int_t num_eta_partitions = chamber->nEtaPartitions();
          ME3IdsKey key3{region_id, station_id, layer_id};

          me_mu_occ_eta_[key3] = bookHist1D(booker,
                                            key3,
                                            "muon_simhit_occ_eta",
                                            "Muon SimHit Eta Occupancy",
                                            16,
                                            eta_range_[station_id * 2 + 0],
                                            eta_range_[station_id * 2 + 1],
                                            "#eta");

          me_mu_occ_phi_[key3] = bookHist1D(
              booker, key3, "muon_simhit_occ_phi", "Muon SimHit Phi Occupancy", 36, -5, 355, "#phi [degrees]");
          if (detail_plot_) {
            me_detail_occ_xy_[key3] = bookXYOccupancy(booker, key3, "simhit", "SimHit");

            me_detail_occ_ieta_[key3] = bookHist1D(booker,
                                                   key3,
                                                   "simhit_occ_ieta",
                                                   "SimHit Occupancy per eta partition",
                                                   num_eta_partitions,
                                                   0.5,
                                                   num_eta_partitions + 0.5,
                                                   "i_{#eta}");

            me_detail_occ_phi_[key3] =
                bookHist1D(booker, key3, "simhit_occ_phi", "SimHit Phi Occupancy", 108, -5, 355, "#phi [degrees]");
          }  // detail plot
        }    // layer loop
      }      // end else
    }        // station loop
  }          // region loop
}

std::tuple<Double_t, Double_t> GEMSimHitValidation::getTOFRange(Int_t station_id) {
  UInt_t start_index = station_id * 2;
  Double_t tof_min = tof_range_[start_index];
  Double_t tof_max = tof_range_[start_index + 1];
  return std::make_tuple(tof_min, tof_max);
}

void GEMSimHitValidation::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  const GEMGeometry* gem = &setup.getData(geomToken_);

  edm::Handle<edm::PSimHitContainer> simhit_container;
  event.getByToken(simhit_token_, simhit_container);
  if (not simhit_container.isValid()) {
    edm::LogError(kLogCategory_) << "Cannot get GEMHits by Token simhitLabel" << std::endl;
    return;
  }

  Int_t total_simhit = 0;
  for (const auto& simhit : *simhit_container.product()) {
    const GEMDetId gemid(simhit.detUnitId());
    total_simhit++;

    if (gem->idToDet(gemid) == nullptr) {
      edm::LogError(kLogCategory_) << "SimHit did not matched with GEM Geometry." << std::endl;
      continue;
    }

    Int_t region_id = gemid.region();
    Int_t station_id = gemid.station();
    Int_t layer_id = gemid.layer();
    Int_t chamber_id = gemid.chamber();
    Int_t roll_id = gemid.roll();
    Int_t num_chambers = gemid.nlayers();

    ME2IdsKey key2{region_id, station_id};
    ME3IdsKey key3{region_id, station_id, layer_id};

    GlobalPoint&& simhit_global_pos = gem->idToDet(gemid)->surface().toGlobal(simhit.localPosition());

    Float_t simhit_g_x = simhit_global_pos.x();
    Float_t simhit_g_y = simhit_global_pos.y();
    Float_t simhit_g_r = simhit_global_pos.perp();
    Float_t simhit_g_abs_z = std::fabs(simhit_global_pos.z());
    Float_t simhit_g_phi = toDegree(simhit_global_pos.phi());
    Float_t simhit_g_eta = std::fabs(simhit_global_pos.eta());

    Float_t energy_loss = kEnergyCF_ * simhit.energyLoss();
    energy_loss = energy_loss > 10 ? 9.9 : energy_loss;
    Float_t tof = simhit.timeOfFlight();
    Int_t pid = std::abs(simhit.particleType());

    // NOTE Fill MonitorElement
    Int_t bin_x = getDetOccBinX(num_chambers, chamber_id, layer_id);

    me_tof_[key3]->Fill(tof);
    me_tof_mu_[key3]->Fill(tof);

    bool is_muon_simhit = isMuonSimHit(simhit);
    if (is_muon_simhit) {
      me_tof_mu_st_[station_id]->Fill(tof);
      me_eloss_mu_[station_id]->Fill(energy_loss);
      me_mu_occ_eta_[key3]->Fill(simhit_g_eta);
      me_mu_occ_phi_[key3]->Fill(simhit_g_phi);
    }

    if (detail_plot_) {
      me_detail_occ_ieta_[key3]->Fill(roll_id);
      me_detail_occ_phi_[key3]->Fill(simhit_g_phi);
      me_detail_occ_xy_[key3]->Fill(simhit_g_x, simhit_g_y);
      me_detail_occ_zr_[region_id]->Fill(simhit_g_abs_z, simhit_g_r);
      me_detail_occ_det_[key2]->Fill(bin_x, roll_id);
      if (!is_muon_simhit) {
        if (find(pids_.begin(), pids_.end(), pid) != pids_.end()) {
          me_detail_eloss_pid_[pid]->Fill(energy_loss);
          me_detail_tof_pid_[pid][station_id]->Fill(tof);
        } else {
          me_detail_eloss_pid_[0]->Fill(energy_loss);
          me_detail_tof_pid_[0][station_id]->Fill(tof);
        }
      }

      if (is_muon_simhit) {
        me_detail_eloss_mu_[key3]->Fill(energy_loss);
        me_detail_mu_occ_det_[key2]->Fill(bin_x, roll_id);
      }

    }  // detail_plot
  }    // simhit loop
  me_total_hits_->Fill(total_simhit);
}
