#include "DetectorDescription/DDCMS/interface/MuonNumbering.h"
#include "Geometry/MuonNumbering/interface/MuonBaseNumber.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

using namespace cms;
using namespace std;
using namespace edm;

const MuonBaseNumber
MuonNumbering::geoHistoryToBaseNumber(const cms::ExpandedNodes &nodes) const {
  MuonBaseNumber num;
  
  int levelPart = get("level");
  int superPart = get("super");
  int basePart  = get("base");
  int startCopyNo = get("xml_starts_with_copyno");

  // some consistency checks
  if(basePart != 1) {
    LogError("Geometry") << "MuonNumbering finds unusual base constant:"
			 << basePart;
  }
  if(superPart < 100) {
    LogError("Geometry") << "MuonNumbering finds unusual super constant:"
			 << superPart;
  }
  if(levelPart < 10*superPart) {
    LogError("Geometry") << "MuonNumbering finds unusual level constant:"
			 << levelPart;
  }
  if((startCopyNo !=0 ) && (startCopyNo != 1)) {
    LogError("Geometry") << "MuonNumbering finds unusual start value for copy numbers:"
			 << startCopyNo;
  }
  int ctr(0);
  for(auto const& it : nodes.tags) {
    int tag = it/levelPart;
    if(tag > 0) {
      int offset = nodes.offsets[ctr];
      int copyno = nodes.copyNos[ctr] + offset%superPart;
      int super = offset/superPart;
      num.addBase(tag, super, copyno - startCopyNo);
    }
    ++ctr;
  }
  return num;
}

const int
MuonNumbering::get(const char* key) const {
  int result(0);
  auto const& it = values.find(key);
  if(it != end(values)) 
    result = it->second;
  return result;
}
