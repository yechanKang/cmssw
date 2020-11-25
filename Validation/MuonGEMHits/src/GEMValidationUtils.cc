#include "Validation/MuonGEMHits/interface/GEMValidationUtils.h"

#include "TString.h"

TString GEMUtils::getSuffixName(Int_t region_id) { return TString::Format("_RE%+d", region_id); }

TString GEMUtils::getSuffixName(Int_t region_id, Int_t station_id) {
  char sign = region_id > 0 ? '+' : '-';
  station_id = abs(station_id * 10 + (station_id != 0));
  return TString::Format("_GE%c%d", sign, station_id);
}

TString GEMUtils::getSuffixName(Int_t region_id, Int_t station_id, Int_t layer_id) {
  char sign = region_id > 0 ? '+' : '-';
  station_id = abs(station_id * 10 + (station_id != 0));
  return TString::Format("_GE%c%d_L%d", sign, station_id, layer_id);
}

TString GEMUtils::getSuffixName(Int_t region_id, Int_t station_id, Int_t layer_id, Int_t roll_id) {
  char sign = region_id > 0 ? '+' : '-';
  station_id = abs(station_id * 10 + (station_id != 0));
  return TString::Format("_GE%c%d_L%d_R%d", sign, station_id, layer_id, roll_id);
}

TString GEMUtils::getSuffixName(const ME2IdsKey& key) {
  auto [region_id, station_id] = key;
  return getSuffixName(region_id, station_id);
}

TString GEMUtils::getSuffixName(const ME3IdsKey& key) {
  auto [region_id, station_id, layer_id] = key;
  return getSuffixName(region_id, station_id, layer_id);
}

TString GEMUtils::getSuffixName(const ME4IdsKey& key) {
  auto [region_id, station_id, layer_id, roll_id] = key;
  return getSuffixName(region_id, station_id, layer_id, roll_id);
}

TString GEMUtils::getSuffixTitle(Int_t region_id) { return TString::Format(" Region %+d", region_id); }

TString GEMUtils::getSuffixTitle(Int_t region_id, Int_t station_id) {
  char sign = region_id > 0 ? '+' : '-';
  station_id = abs(station_id * 10 + (station_id != 0));
  return TString::Format(" GE%c%d", sign, station_id);
}

TString GEMUtils::getSuffixTitle(Int_t region_id, Int_t station_id, Int_t layer_id) {
  char sign = region_id > 0 ? '+' : '-';
  station_id = abs(station_id * 10 + (station_id != 0));
  return TString::Format(" GE%c%d_L%d", sign, station_id, layer_id);
}

TString GEMUtils::getSuffixTitle(Int_t region_id, Int_t station_id, Int_t layer_id, Int_t roll_id) {
  char sign = region_id > 0 ? '+' : '-';
  station_id = abs(station_id * 10 + (station_id != 0));
  return TString::Format(" GE%c%d_L%dR_%d", sign, station_id, layer_id, roll_id);
}

TString GEMUtils::getSuffixTitle(const ME2IdsKey& key) {
  auto [region_id, station_id] = key;
  return getSuffixTitle(region_id, station_id);
}

TString GEMUtils::getSuffixTitle(const ME3IdsKey& key) {
  auto [region_id, station_id, layer_id] = key;
  return getSuffixTitle(region_id, station_id, layer_id);
}

TString GEMUtils::getSuffixTitle(const ME4IdsKey& key) {
  auto [region_id, station_id, layer_id, roll_id] = key;
  return getSuffixTitle(region_id, station_id, layer_id, roll_id);
}
