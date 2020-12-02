#ifndef Validation_MuonGEMRecHits_MuonGEMDigisHarvestor_h
#define Validation_MuonGEMRecHits_MuonGEMDigisHarvestor_h

#include "Validation/MuonGEMHits/interface/MuonGEMBaseHarvestor.h"

class MuonGEMRecHitsHarvestor : public MuonGEMBaseHarvestor {
public:
  explicit MuonGEMRecHitsHarvestor(const edm::ParameterSet&);
  ~MuonGEMRecHitsHarvestor() override;
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

  void bookMean1D(DQMStore::IBooker& booker,
                  DQMStore::IGetter& getter,
                  const TString residual_path,
                  const TString folder,
                  const TString name,
                  const TString title);

private:
  // NOTE to make it compatible to both full geometry and slice test
  std::vector<Int_t> region_ids_, station_ids_, layer_ids_;
};

#endif  // Validation_MuonGEMRecHits_MuonGEMDigisHarvestor_h
