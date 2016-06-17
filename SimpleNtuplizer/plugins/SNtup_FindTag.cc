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
        const reco::RecoCandidate& object,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

  tp_mll = -1;
  tp_tagpt = -1;
  tp_tageta = -1;
  
  edm::Handle<edm::TriggerResults> hTrgRes;
  iEvent.getByToken(HLTTag_token_,hTrgRes);

  edm::Handle<trigger::TriggerEvent> hTrgEvt;
  iEvent.getByToken(HLTObjTag_token_,hTrgEvt);

  // Discard events that did not pass the list of triggers
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);  
  for (std::string pattern : elecTrig_) {
    if(edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
      std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
      if(matches.empty()) {
	return;
      }
    }
  }
    
  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> electrons;  // For AODSIM
  iEvent.getByToken(electronToken_, electrons);

  // Get electron ID decision
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);

  // Loop over electrons
  size_t elsIndex = 0;
  for( reco::GsfElectronCollection::const_iterator el = electrons->begin(); el != electrons->end(); el++, elsIndex++) {

    // Discard electrons that did not match the trigger object or that did not pass tightID
    bool triggerMatch = false;    
    bool passID = false;
    bool passKinematics = false;

    passKinematics = el->pt() > 35. && fabs(el->eta()) < 1.4;
    if (!passKinematics) continue;
    
    const edm::Ptr<reco::GsfElectron> elPtr(electrons, el - electrons->begin() );
    passID = (*tight_id_decisions)[elPtr];
    if (!passID) continue;

    for (std::string pattern : elecFilt_) {
      edm::InputTag filterTag(pattern,"","HLT");
      if(hTrgEvt->filterIndex(filterTag) < hTrgEvt->sizeFilters()) {
	const trigger::TriggerObjectCollection& toc(hTrgEvt->getObjects());
	const trigger::Keys& keys(hTrgEvt->filterKeys(hTrgEvt->filterIndex(filterTag)));
	
	for(size_t hlto=0; hlto<keys.size(); hlto++) {
	  trigger::size_type hltf = keys[hlto];
	  const trigger::TriggerObject& tobj(toc[hltf]);
	  std::cout << el->eta() << " " << el->phi() << " " << tobj.eta() << " " << tobj.phi() << std::endl;
	  if(reco::deltaR(el->eta(),el->phi(),tobj.eta(),tobj.phi()) < 0.2) {
	    triggerMatch = true;
	  }
	}
      }
    }
    if (!triggerMatch) continue;

    float mass = (el->p4()+object.p4()).mass();
    std::cout << "Found a tag with mass " << mass << std::endl;
    if (fabs(mass-91.2) < fabs(tp_mll-91.2)) {
      tp_mll = mass;
      tp_tagpt = el->pt();
      tp_tageta = el->eta();
    }
  }
  return;
    
}
