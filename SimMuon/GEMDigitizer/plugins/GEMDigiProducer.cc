#ifndef SimMuon_GEMDigitizer_GEMDigiProducer_h
#define SimMuon_GEMDigitizer_GEMDigiProducer_h

#include "SimMuon/GEMDigitizer/interface/GEMDigiModule.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include <iostream>
using std::cout;
using std::endl;
#include <random>

#include <sstream>
#include <string>
#include <map>
#include <vector>

namespace CLHEP {
  class HepRandomEngine;
}

class GEMDigiProducer : public edm::stream::EDProducer<> {
public:
  typedef edm::DetSetVector<StripDigiSimLink> StripDigiSimLinks;

  typedef edm::DetSetVector<GEMDigiSimLink> GEMDigiSimLinks;

  explicit GEMDigiProducer(const edm::ParameterSet& ps);

  ~GEMDigiProducer() override;

  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  //Name of Collection used for create the XF
  edm::EDGetTokenT<CrossingFrame<PSimHit> > cf_token;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> geom_token_;

  const GEMGeometry* geometry_;

  std::unique_ptr<GEMDigiModule> gemDigiModule_;
};

GEMDigiProducer::GEMDigiProducer(const edm::ParameterSet& ps) : gemDigiModule_(std::make_unique<GEMDigiModule>(ps)) {
  produces<GEMDigiCollection>();
  produces<StripDigiSimLinks>("GEM");
  produces<GEMDigiSimLinks>("GEM");

  edm::Service<edm::RandomNumberGenerator> rng;
  if (!rng.isAvailable()) {
    throw cms::Exception("Configuration")
        << "GEMDigiProducer::GEMDigiProducer() - RandomNumberGeneratorService is not present in configuration file.\n"
        << "Add the service in the configuration file or remove the modules that require it.";
  }

  std::string mix_(ps.getParameter<std::string>("mixLabel"));
  std::string collection_(ps.getParameter<std::string>("inputCollection"));

  cf_token = consumes<CrossingFrame<PSimHit> >(edm::InputTag(mix_, collection_));
  geom_token_ = esConsumes<GEMGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

GEMDigiProducer::~GEMDigiProducer() = default;

void GEMDigiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("inputCollection", "g4SimHitsMuonGEMHits");
  desc.add<std::string>("mixLabel", "mix");

  desc.add<double>("signalPropagationSpeed", 0.66);
  desc.add<double>("timeResolution", 5.);
  desc.add<double>("timeJitter", 1.0);
  desc.add<double>("averageShapingTime", 50.0);
  desc.add<double>("averageEfficiency", 0.98);
  desc.add<double>("averageNoiseRate", 0.001);
  // intrinsic noise rate (Hz/cm^2)

  // in terms of 25 ns
  desc.add<int>("minBunch", -5);
  desc.add<int>("maxBunch", 3);

  desc.add<bool>("fixedRollRadius", true);
  // Uses fixed radius in the center of the roll
  desc.add<bool>("digitizeOnlyMuons", false);
  desc.add<bool>("simulateBkgNoise", false);
  // false == No background simulation
  desc.add<bool>("simulateNoiseCLS", true);
  desc.add<bool>("simulateElectronBkg", true);
  // flase == simulate only neutral bkg
  desc.add<bool>("simulateIntrinsicNoise", false);

  desc.add<double>("instLumi", 7.5);
  // in units of 1E34 cm^-2 s^-1. Internally the background is parmetrized from FLUKA+GEANT result at 5E+34 (PU 140). We are adding a 1.5 factor for PU 200
  desc.add<double>("rateFact", 1.0);
  // We are adding also a safety factor of 2 to tak into account the new beam pipe effect (not yet known). Hits can be thrown away later at re-digi step. Parameters are kept in sync with the ones used in the GEM digitizer
  desc.add<double>("bxWidth", 25E-9);
  desc.add<double>("referenceInstLumi", 5.);
  // referecne inst. luminosity 5E+34 cm^-2s^-1
  desc.add<double>("resolutionX", 0.03);

  // The follwing parameters are needed to model the background contribution
  // The parameters have been obtained after the fit of th perdicted by FLUKA
  // By default the backgroundmodeling with these parameters should be disabled with
  // the 9_2_X release setting simulateBkgNoise = false
  desc.add<double>("GE11ModNeuBkgParam0", 5710.23);
  desc.add<double>("GE11ModNeuBkgParam1", -43.3928);
  desc.add<double>("GE11ModNeuBkgParam2", 0.0863681);
  desc.add<double>("GE21ModNeuBkgParam0", 1440.44);
  desc.add<double>("GE21ModNeuBkgParam1", -7.48607);
  desc.add<double>("GE21ModNeuBkgParam2", 0.0103078);
  desc.add<double>("GE11ElecBkgParam0", 406.249);
  desc.add<double>("GE11ElecBkgParam1", -2.90939);
  desc.add<double>("GE11ElecBkgParam2", 0.00548191);
  desc.add<double>("GE21ElecBkgParam0", 97.0505);
  desc.add<double>("GE21ElecBkgParam1", -43.3928);
  desc.add<double>("GE21ElecBkgParam2", 00.000550599);

  descriptions.add("simMuonGEMDigisDef", desc);
}

void GEMDigiProducer::beginRun(const edm::Run&, const edm::EventSetup& eventSetup) {
  edm::ESHandle<GEMGeometry> hGeom = eventSetup.getHandle(geom_token_);
  gemDigiModule_->setGeometry(&*hGeom);
  geometry_ = &*hGeom;
}

void GEMDigiProducer::produce(edm::Event& e, const edm::EventSetup& eventSetup) {
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  edm::Handle<CrossingFrame<PSimHit> > cf;
  e.getByToken(cf_token, cf);

  MixCollection<PSimHit> hits{cf.product()};
  
  // Scripts for the double polarization
  
  std::map<GEMDetId, bool> etaMask;
  std::random_device rd; // obtain a random number from hardware
  std::mt19937 eng(rd()); // seed the generator
  std::uniform_real_distribution<> randomFloat(0, 1);

  float VALUE[8] = { 0.212997435,
                     0.27290354,
                     0.29952727,
                     0.379401715,
                     0.4792445,
                     0.59239915,
                     0.965145475,
                     1.25136305
                     };

  // Crosstalk with deadtime multiflied by some number
  float multiplier = 1.; // 50BX
  //float multiplier = 0.5; // 25BX
  //float multiplier = 0.1; // 50BX with 10 times reduced rate
  //float multiplier = 0.05; // 25BX with 10 times reduced rate

  float percentage = 0.01;
  
  for (const auto& roll : geometry_->etaPartitions()) {
    // The values from gemDPG_phase2_20200303.pdf X ===> Updated : see last slide in https://docs.google.com/presentation/d/1AfmQ1uDBZWHX0DdL1rtaHb3ILIKKVu_KFY_vjYsGj00/edit#slide=id.g8a6032e479_0_1
    const GEMDetId detId(roll->id());
    if (detId.station() != 2) continue;
    auto randomNumber = randomFloat(eng);
    int idx = (detId.roll()-1) / 2;
    if ( randomNumber < percentage * VALUE[idx] * multiplier ) etaMask[detId] = true;
  }

  //////////////////// double pol end

  // Create empty output
  auto digis = std::make_unique<GEMDigiCollection>();
  auto stripDigiSimLinks = std::make_unique<StripDigiSimLinks>();
  auto gemDigiSimLinks = std::make_unique<GEMDigiSimLinks>();

  // arrange the hits by eta partition
  std::map<uint32_t, edm::PSimHitContainer> hitMap;
  for (const auto& hit : hits) {
    hitMap[GEMDetId(hit.detUnitId()).rawId()].emplace_back(hit);
  }

  // simulate signal and noise for each eta partition
  const auto& etaPartitions(geometry_->etaPartitions());

  for (const auto& roll : etaPartitions) {
    const GEMDetId detId(roll->id());
    const uint32_t rawId(detId.rawId());
    const auto& simHits(hitMap[rawId]);

    LogDebug("GEMDigiProducer") << "GEMDigiProducer: found " << simHits.size() << " hit(s) in eta partition" << rawId;

    //if (etaMask[detId] == true) cout << "Skip digis on a dead board : " << detId << endl;
    //else {
    if (etaMask[detId] != true) {
      gemDigiModule_->simulate(roll, simHits, engine);
      gemDigiModule_->fillDigis(rawId, *digis);
      (*stripDigiSimLinks).insert(gemDigiModule_->stripDigiSimLinks());
      (*gemDigiSimLinks).insert(gemDigiModule_->gemDigiSimLinks());
    }
  }

  // store them in the event
  e.put(std::move(digis));
  e.put(std::move(stripDigiSimLinks), "GEM");
  e.put(std::move(gemDigiSimLinks), "GEM");
}

DEFINE_FWK_MODULE(GEMDigiProducer);
#endif
