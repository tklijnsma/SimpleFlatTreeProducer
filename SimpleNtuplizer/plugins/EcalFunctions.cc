#include "EcalFunctions.h"

float nTupler::getNormedIX(const DetId& id)
{
  if(id.det()==DetId::Ecal && id.subdetId()==EcalEndcap){
    EEDetId eeId(id);
    int iXNorm  = eeId.ix()-50;
    if(iXNorm<=0) iXNorm--;
    return iXNorm;
  }
  return 0;
}

float nTupler::getNormedIY(const DetId& id)
{
  if(id.det()==DetId::Ecal && id.subdetId()==EcalEndcap){
    EEDetId eeId(id);
    int iYNorm  = eeId.iy()-50;
    if(iYNorm<=0) iYNorm--;
    return iYNorm;
  }
  return 0;
}


float nTupler::getIEta(const DetId& id)
{
  if(id.det()==DetId::Ecal){
    if(id.subdetId()==EcalBarrel){
      EBDetId ebId(id);
      return ebId.ieta();
    }else if(id.subdetId()==EcalEndcap){
      float iXNorm = nTupler::getNormedIX(id);
      float iYNorm = nTupler::getNormedIY(id);

      return std::sqrt(iXNorm*iXNorm+iYNorm*iYNorm);
    }
  }
  return 0.;
}

float nTupler::getIPhi(const DetId& id)
{
  if(id.det()==DetId::Ecal){
    if(id.subdetId()==EcalBarrel){
      EBDetId ebId(id);
      return ebId.iphi();
    }
  }
  return 0.;
}

//nr crystals crysId is away from orgin id in eta
float nTupler::getNrCrysDiffInEta(const DetId& crysId,const DetId& orginId)
{
  float crysIEta = nTupler::getIEta(crysId);
  float orginIEta = nTupler::getIEta(orginId);
  bool isBarrel = orginId.subdetId()==EcalBarrel;

  float nrCrysDiff = crysIEta-orginIEta;

  //no iEta=0 in barrel, so if go from positive to negative
  //need to reduce abs(detEta) by 1
  if(isBarrel){
    if(crysIEta*orginIEta<0){ // -1 to 1 transition
      if(crysIEta>0) nrCrysDiff--;
      else nrCrysDiff++;
    }
  }
  return nrCrysDiff;
}

//nr crystals crysId is away from orgin id in phi
float nTupler::getNrCrysDiffInPhi(const DetId& crysId,const DetId& orginId)
{
  float crysIPhi = nTupler::getIPhi(crysId);
  float orginIPhi = nTupler::getIPhi(orginId);
  bool isBarrel = orginId.subdetId()==EcalBarrel;

  float nrCrysDiff = crysIPhi-orginIPhi;

  if(isBarrel){ //if barrel, need to map into 0-180
    if (nrCrysDiff > + 180) { nrCrysDiff = nrCrysDiff - 360; }
    if (nrCrysDiff < - 180) { nrCrysDiff = nrCrysDiff + 360; }
  }
  return nrCrysDiff;
}

