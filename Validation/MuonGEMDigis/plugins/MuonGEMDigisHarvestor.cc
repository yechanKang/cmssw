#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Validation/MuonGEMDigis/plugins/MuonGEMDigisHarvestor.h"

MuonGEMDigisHarvestor::MuonGEMDigisHarvestor(const edm::ParameterSet& pset)
    : MuonGEMBaseHarvestor(pset, "MuonGEMDigisHarvestor") {
  // to make it compatible to both full geometry and slice test
  region_ids_ = pset.getUntrackedParameter<std::vector<Int_t> >("regionIds");
  station_ids_ = pset.getUntrackedParameter<std::vector<Int_t> >("stationIds");
  layer_ids_ = pset.getUntrackedParameter<std::vector<Int_t> >("layerIds");
}

MuonGEMDigisHarvestor::~MuonGEMDigisHarvestor() {}

void MuonGEMDigisHarvestor::dqmEndJob(DQMStore::IBooker& booker, DQMStore::IGetter& getter) {
  TString occ_folder = "MuonGEMDigisV/GEMDigisTask/Strip/Occupancy/";
  TString eff_folder = "MuonGEMDigisV/GEMDigisTask/Strip/Efficiency/";

  for (Int_t region_id : region_ids_) {
    TString name_suffix_re = GEMUtils::getSuffixName(region_id);
    TString title_suffix_re = GEMUtils::getSuffixTitle(region_id);

    // NOTE eta efficiency
    TString strip_eta_name = "matched_strip_occ_eta" + name_suffix_re;
    TString simhit_eta_name = "muon_simhit_occ_eta" + name_suffix_re;
    TString strip_eta_path = occ_folder + strip_eta_name;
    TString simhit_eta_path = occ_folder + simhit_eta_name;
    TString eff_eta_name = "eff_eta" + name_suffix_re;
    TString eff_eta_title = "Eta Efficiency (Muon Only) :" + title_suffix_re;

    bookEff1D(booker, getter, strip_eta_path, simhit_eta_path, eff_folder, eff_eta_name, eff_eta_title);

    for (Int_t station_id : station_ids_) {
      TString name_suffix_re_st = GEMUtils::getSuffixName(region_id, station_id);
      TString title_suffix_re_st = GEMUtils::getSuffixTitle(region_id, station_id);

      // NOTE phi efficiency
      TString strip_phi_name = "matched_strip_occ_phi" + name_suffix_re_st;
      TString simhit_phi_name = "muon_simhit_occ_phi" + name_suffix_re_st;
      TString strip_phi_path = occ_folder + strip_phi_name;
      TString simhit_phi_path = occ_folder + simhit_phi_name;
      TString eff_phi_name = "eff_phi" + name_suffix_re_st;
      TString eff_phi_title = "Phi Efficiency (Muon Only) :" + title_suffix_re;

      bookEff1D(booker, getter, strip_phi_path, simhit_phi_path, eff_folder, eff_phi_name, eff_phi_title);

      // NOTE Detector Component efficiency
      TString strip_det_name = "matched_strip_occ_det" + name_suffix_re_st;
      TString simhit_det_name = "muon_simhit_occ_det" + name_suffix_re_st;
      TString strip_det_path = occ_folder + strip_det_name;
      TString simhit_det_path = occ_folder + simhit_det_name;
      TString eff_det_name = "eff_det" + name_suffix_re_st;
      TString eff_det_title = "Detector Component Efficiency (Muon Only) :" + title_suffix_re_st;

      bookEff2D(booker, getter, strip_det_path, simhit_det_path, eff_folder, eff_det_name, eff_det_title);

    }  // statino loop
  }    // region loop

  TString pad_occ_folder = "MuonGEMDigisV/GEMDigisTask/Pad/Occupancy/";
  TString pad_eff_folder = "MuonGEMDigisV/GEMDigisTask/Pad/Efficiency/";
  TString cluster_occ_folder = "MuonGEMDigisV/GEMDigisTask/PadCluster/Occupancy/";
  TString cluster_eff_folder = "MuonGEMDigisV/GEMDigisTask/PadCluster/Efficiency/";

  for (Int_t region_id : region_ids_) {
    TString name_suffix_re = GEMUtils::getSuffixName(region_id);
    TString title_suffix_re = GEMUtils::getSuffixTitle(region_id);

    for (Int_t station_id : station_ids_) {
      TString name_suffix_re_st = GEMUtils::getSuffixName(region_id, station_id);
      TString title_suffix_re_st = GEMUtils::getSuffixTitle(region_id, station_id);

      // NOTE Detector Component efficiency
      TString pad_det_name = "pad_occ_det" + name_suffix_re_st;
      TString strip_det_name = "pad_strip_occ_det" + name_suffix_re_st;
      TString pad_det_path = pad_occ_folder + pad_det_name;
      TString strip_det_path = pad_occ_folder + strip_det_name;
      TString eff_det_name = "eff_det" + name_suffix_re_st;
      TString eff_det_title = "Detector Component Efficiency :" + title_suffix_re_st;

      bookEff2D(booker, getter, pad_det_path, strip_det_path, pad_eff_folder, eff_det_name, eff_det_title);

      TString resolved_det_name = "pad_resolved_occ_det" + name_suffix_re_st;
      TString resolved_det_path = cluster_occ_folder + resolved_det_name;
      TString eff_resolved_det_name = "eff_det" + name_suffix_re_st;

      bookEff2D(booker, getter, resolved_det_path, pad_det_path, cluster_eff_folder, eff_resolved_det_name, eff_det_title);
    }  // statino loop
  }    // region loop
}
