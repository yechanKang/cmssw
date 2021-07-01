
/** Derived from DTGeometryAnalyzer by Nicola Amapane
 *
 *  \author M. Maggi - INFN Bari
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include <memory>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <set>

class GEMAreaAnalyser : public edm::one::EDAnalyzer<> {
public:
  GEMAreaAnalyser(const edm::ParameterSet& pset);

  ~GEMAreaAnalyser() override;

  void beginJob() override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
  void endJob() override {}

private:
  const std::string& myName() { return myName_; }

  const int dashedLineWidth_;
  const std::string dashedLine_;
  const std::string myName_;
  std::ofstream ofos;
};

using namespace std;
GEMAreaAnalyser::GEMAreaAnalyser(const edm::ParameterSet& /*iConfig*/)
    : dashedLineWidth_(104), dashedLine_(std::string(dashedLineWidth_, '-')), myName_("GEMAreaAnalyser") {
  ofos.open("GEMtestOutput.out");
  ofos << "st,chIdx,ly,ieta,trArea" << endl;
}

GEMAreaAnalyser::~GEMAreaAnalyser() {
  ofos.close();
}

void GEMAreaAnalyser::analyze(const edm::Event& /*iEvent*/, const edm::EventSetup& iSetup) {
  edm::ESHandle<GEMGeometry> pDD;
  iSetup.get<MuonGeometryRecord>().get(pDD);

  for (auto station : pDD->regions()[0]->stations()) {
    for (auto ring : station->rings()) {
      auto sch_even = ring->superChamber(2);
      auto sch_odd = ring->superChamber(1);
      for (auto sch : {sch_even, sch_odd}) {
        for (auto ch : sch->chambers()) {
          for (auto roll : ch->etaPartitions()) {
            GEMDetId rId(roll->id());

            ofos << rId.station() << "," << rId.chamber()%2 << "," << rId.layer() << "," << rId.roll() << ",";

            const GEMStripTopology* top_(dynamic_cast<const GEMStripTopology*>(&(roll->topology())));

            float striplength = top_->stripLength();
            float trStripArea = (roll->pitch()) * striplength;

            int nstrips = roll->nstrips();

            float trArea = trStripArea * nstrips;

            ofos << trArea << endl;

          }
        }
      }
    }
  }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GEMAreaAnalyser);
