#include "SimMuon/GEMDigitizer/interface/GEMNoiseModel.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <cmath>
#include <utility>
#include <map>

GEMNoiseModel::GEMNoiseModel(const edm::ParameterSet& config, GEMDigiModule* digiModule) :
GEMDigiModel(config, digiModule)
, averageNoiseRate_(config.getParameter<double> ("averageNoiseRate"))
, bxwidth_(config.getParameter<int> ("bxwidth"))
, minBunch_(config.getParameter<int> ("minBunch"))
, maxBunch_(config.getParameter<int> ("maxBunch"))
{
}

GEMNoiseModel::~GEMNoiseModel()
{
}

void GEMNoiseModel::simulate(const GEMEtaPartition* roll, const edm::PSimHitContainer&, CLHEP::HepRandomEngine* engine)
{
  const GEMDetId& gemId(roll->id());
  const int nstrips(roll->nstrips());
  double trStripArea(0.0);
  if (gemId.region() == 0)
  {
    throw cms::Exception("Geometry")
        << "GEMSynchronizer::simulate() - this GEM id is from barrel, which cannot happen.";
  }
  const TrapezoidalStripTopology* top_(dynamic_cast<const TrapezoidalStripTopology*>(&(roll->topology())));
  const float striplength(top_->stripLength());
  trStripArea = (roll->pitch()) * striplength;
  const int nBxing(maxBunch_ - minBunch_ + 1);
  //simulate intrinsic noise
  const double aveIntrinsicNoisePerStrip(averageNoiseRate_ * nBxing * bxwidth_ * trStripArea * 1.0e-9);
  for(int j = 0; j < nstrips; ++j)
  {
    CLHEP::RandPoissonQ randPoissonQ(*engine, aveIntrinsicNoisePerStrip);
    const int n_intrHits(randPoissonQ.fire());
  
    for (int k = 0; k < n_intrHits; k++ )
      {
      const int time_hit(static_cast<int>(CLHEP::RandFlat::shoot(engine, nBxing)) + minBunch_);
      std::pair<int, int> digi(k+1,time_hit);
      digiModule_->emplaceStrip(digi);
    }
  }
  //end simulate intrinsic noise

  return;
}
