#ifndef RecoLocalMuon_GEMRecHit_GEMEtaPartitionMask_h
#define RecoLocalMuon_GEMRecHit_GEMEtaPartitionMask_h

#include <bitset>
#include <vector>

// strip numbering start from 1 in simulation. if strip numbering has fixed, it should be fixed.
const int maskSIZE = 769; 
typedef std::bitset<maskSIZE> EtaPartitionMask;

#endif
