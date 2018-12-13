#ifndef SimMuon_GEMDigitizer_GEMSignalModel_h
#define SimMuon_GEMDigitizer_GEMSignalModel_h

/** 
 * \class GEMSignalModel
 *
 * Class for the GEM strip response simulation based on a very simple model
 * Originally comes from GEMSimpleModel
 *
 * \author Sven Dildick
 * \modified by Roumyana Hadjiiska
 * \splitted by Yechan Kang
 */

#include "SimMuon/GEMDigitizer/interface/GEMDigiModel.h"

class GEMGeometry;

namespace CLHEP
{
  class HepRandomEngine;
}

class GEMSignalModel: public GEMDigiModel
{
public:

  GEMSignalModel(const edm::ParameterSet&, GEMDigiModule*);

  ~GEMSignalModel() override;

  void simulate(const GEMEtaPartition*, const edm::PSimHitContainer&, CLHEP::HepRandomEngine*) override;

  int getSimHitBx(const PSimHit*, CLHEP::HepRandomEngine*);

  std::vector<std::pair<int,int> > 
    simulateClustering(const GEMEtaPartition*, const PSimHit*, const int, CLHEP::HepRandomEngine*);

private:

  double averageEfficiency_;
  double averageShapingTime_;
  double timeResolution_;
  double timeJitter_;
  double signalPropagationSpeed_;
  int bxwidth_;
  bool digitizeOnlyMuons_;
  double resolutionX_; 

};
#endif


