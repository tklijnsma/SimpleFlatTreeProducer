#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

namespace nTupler {
  float getNormedIX(const DetId& id);
  float getNormedIY(const DetId& id);
  float getIEta(const DetId& id);
  float getIPhi(const DetId& id);
  float getNrCrysDiffInEta(const DetId& crysId,const DetId& orginId);
  float getNrCrysDiffInPhi(const DetId& crysId,const DetId& orginId);
}

