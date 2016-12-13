#include "SimpleNtuplizer.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "DataFormats/L1TGlobal/interface/GlobalObjectMap.h"
#include "DataFormats/L1TGlobal/interface/GlobalObjectMapRecord.h"
#include "DataFormats/L1TGlobal/interface/GlobalObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <boost/foreach.hpp>
#include <iostream>

void SimpleNtuplizer::findTag(
        const reco::RecoCandidate& object, float corrToRaw,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

  tp_mll = -1;
  tp_ptll = -1;
  tp_tagpt = -1;
  tp_tageta = -1;
  tp_tagphi = -1;

  edm::Handle<edm::TriggerResults> hTrgRes;
  edm::Handle<pat::TriggerObjectStandAloneCollection> hTrgEvt;

  //  std::cout << "Trying to get trigger info" << std::endl;
  if (isData_) {

    // Get trigger information
    iEvent.getByToken(HLTTag_token_,hTrgRes);    
    iEvent.getByToken(HLTObjTag_token_,hTrgEvt);
    
    // Discard events that did not pass the list of triggers
    const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);  
    for (std::string pattern : elecTrig_) {
      //      std::cout << pattern.c_str() << std::endl;
      if(edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
	std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
	if(matches.empty()) {
	  return;
	}
      }
    }
  }

  //  std::cout << "Trying to get electrons" << std::endl;
  // Get electron collection
  edm::Handle<edm::View<pat::Electron > > electrons;  // For AODSIM
  iEvent.getByToken(electronToken_, electrons);

  // Get electron ID decision
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);

  // Loop over electrons
  size_t elsIndex = 0;
  for( edm::View<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); el++, elsIndex++) {

    // Discard electrons that did not match the trigger object or that did not pass tightID
    bool triggerMatch = false;    
    bool passID = false;
    bool passKinematics = false;

    float rawEnergy = el->superCluster()->rawEnergy();
    float energy = el->energy();
    float uncorr = rawEnergy/energy;
    
    auto electron = (*el);
    passKinematics = electron.pt() > 35. && fabs(electron.eta()) < 1.4;
    //    std::cout << "pass kinematics " << passKinematics << std::endl;
    if (!passKinematics) continue;
    
    const edm::Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );
    passID = (*tight_id_decisions)[elPtr];
    //    std::cout << "pass ID " << passID << std::endl;
    if (!passID) continue;

    if (isData_) {
      pat::TriggerObjectStandAlone TO;
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);
      for ( uint i = 0; i < hTrgEvt->size(); i++ ) {
	TO = hTrgEvt->at(i);
	TO.unpackPathNames(triggerNames);

	bool isTrigger = false;
	bool isFilter = false;
	for (std::string pattern : elecTrig_) {	    
	  if (TO.hasPathName(pattern, false)) { isTrigger = true; break; }
	}
	if (!isTrigger) continue;
	for (std::string pattern : elecFilt_) {
	  if (TO.hasFilterLabel(pattern)) { isFilter = true; break; }
	}
	if (!isFilter) continue;
	if(reco::deltaR(electron.eta(),electron.phi(),TO.eta(),TO.phi()) < 0.2) {
	      triggerMatch = true;
	}
      }
      //      std::cout << "pass trigger " << triggerMatch << std::endl;
      if (!triggerMatch) continue;
    }
    
    float mass = (electron.p4()*uncorr+object.p4()*corrToRaw).mass();
    float pt = (electron.p4()*uncorr+object.p4()*corrToRaw).pt();
    //    std::cout << "Found a tag with mass " << mass << std::endl;
    if (fabs(mass-91.2) < fabs(tp_mll-91.2)) {
      tp_mll = mass;
      tp_ptll = pt;
      tp_tagpt = electron.pt()*uncorr;
      tp_tageta = electron.eta();
      tp_tagphi = electron.phi();
    }
  }
  return;
    
}
