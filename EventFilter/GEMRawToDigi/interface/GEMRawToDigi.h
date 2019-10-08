#ifndef EventFilter_GEMRawToDigi_GEMRawToDigi_h
#define EventFilter_GEMRawToDigi_GEMRawToDigi_h

/** \class GEMRawToDigi
 *  \author J. Lee, Yechan Kang - UoS
 */
#include <memory>
#include "EventFilter/GEMRawToDigi/interface/AMC13Event.h"

class GEMRawToDigi {
public:
  /// Constructor
  GEMRawToDigi(){};

  std::unique_ptr<gem::AMC13Event> convertWordToAMC13Event(const uint64_t* word);

private:
};
#endif
