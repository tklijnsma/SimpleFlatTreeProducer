// -*- C++ -*-
//
// Package:    SimpleNtuplizer
// Class:      SimpleNtuplizer
// 
/**\class SimpleNtuplizer SimpleNtuplizer.cc SimpleNtuplizer/plugins/SimpleNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


#include "SimpleNtuplizer.h"


//######################################
//# Constructors and Destructor
//######################################

// Constructor
SimpleNtuplizer::SimpleNtuplizer(const edm::ParameterSet& iConfig):
    // All tokens given in the python config!
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    electronToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    photonToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
    {

    std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

    edm::Service<TFileService> fs;


    // =====================================
    // Making event variable tree and setting branches
    // Event variables include the quantities of which there are exactly one per event

    eventTree_ = fs->make<TTree> ("EventTree", "Per event data");

    // Event variables
    eventTree_->Branch( "nPV",              &nPV_,             "nPV/I"   );
    eventTree_->Branch( "nElectrons",       &nElectrons_,      "nEle/I"  );
    eventTree_->Branch( "nPhotons",         &nPhotons_,        "nPho/I"  );
    eventTree_->Branch( "nPhotonsMatched_", &nPhotonsMatched_, "nPhoMatched/I"  );


    // =====================================
    // Making electron tree and setting branches

    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    // Electron variables
    electronTree_->Branch( "pt_",                 &pt_    );
    electronTree_->Branch( "etaSC_",              &etaSC_ );
    electronTree_->Branch( "phiSC_",              &phiSC_ );
    electronTree_->Branch( "rawEnergySC_",        &rawEnergySC_ );
    electronTree_->Branch( "etaWidthSC_",         &etaWidthSC_ );
    electronTree_->Branch( "phiWidthSC_",         &phiWidthSC_ );
    electronTree_->Branch( "r9SS_",               &r9SS_ );
    electronTree_->Branch( "seedEnergySS_",       &seedEnergySS_ );
    electronTree_->Branch( "eMaxSS_",             &eMaxSS_ );
    electronTree_->Branch( "e2ndSS_",             &e2ndSS_ );
    electronTree_->Branch( "eHorizontalSS_",      &eHorizontalSS_ );
    electronTree_->Branch( "eVerticalSS_",        &eVerticalSS_ );
    electronTree_->Branch( "sigmaIetaIetaSS_",    &sigmaIetaIetaSS_ );
    electronTree_->Branch( "sigmaIetaIphiSS_",    &sigmaIetaIphiSS_ );
    electronTree_->Branch( "sigmaIphiIphiSS_",    &sigmaIphiIphiSS_ );
    electronTree_->Branch( "preshowerEnergy_",    &preshowerEnergy_ );
    electronTree_->Branch( "numberOfClustersSC_", &numberOfClustersSC_ );
    electronTree_->Branch( "isEB_",               &isEB_ );

    // These for EB electrons
    electronTree_->Branch( "iEtaCoordinate_",     &iEtaCoordinate_ );         
    electronTree_->Branch( "iPhiCoordinate_",     &iPhiCoordinate_ );         
    electronTree_->Branch( "cryEtaCoordinate_",   &cryEtaCoordinate_ );
    electronTree_->Branch( "cryPhiCoordinate_",   &cryPhiCoordinate_ );

    // These for EE electrons
    electronTree_->Branch( "iXCoordinate_",       &iXCoordinate_ );         
    electronTree_->Branch( "iYCoordinate_",       &iYCoordinate_ );         
    electronTree_->Branch( "cryXCoordinate_",     &cryXCoordinate_ );
    electronTree_->Branch( "cryYCoordinate_",     &cryYCoordinate_ );

    // Cluster variables
    electronTree_->Branch( "MaxDRclusterDR_",        &MaxDRclusterDR_ );
    electronTree_->Branch( "MaxDRclusterDPhi_",      &MaxDRclusterDPhi_ );
    electronTree_->Branch( "MaxDRclusterDEta_",      &MaxDRclusterDEta_ );
    electronTree_->Branch( "MaxDRclusterRawEnergy_", &MaxDRclusterRawEnergy_ );

    electronTree_->Branch( "clusterRawEnergy_",      &clusterRawEnergy_ );        
    electronTree_->Branch( "clusterDPhiToSeed_",     &clusterDPhiToSeed_ );         
    electronTree_->Branch( "clusterDEtaToSeed_",     &clusterDEtaToSeed_ );

    electronTree_->Branch( "match_dR",     &match_dR     );
    electronTree_->Branch( "match_dE",     &match_dE     );
    electronTree_->Branch( "match_dRdE",   &match_dRdE   );

    electronTree_->Branch( "gen_pt",     &gen_pt     );
    electronTree_->Branch( "gen_phi",    &gen_phi    );
    electronTree_->Branch( "gen_eta",    &gen_eta    );
    electronTree_->Branch( "gen_M",      &gen_M      );
    electronTree_->Branch( "gen_E",      &gen_E      );
    electronTree_->Branch( "gen_pdgId",  &gen_pdgId  );
    electronTree_->Branch( "gen_status", &gen_status );


    // =====================================
    // Making E-p tree and setting branches

    EpTree_ = fs->make<TTree> ("EpTree", "E-p data");

    EpTree_->Branch( "totEnergyEp_",         &totEnergyEp_ );
    EpTree_->Branch( "epEp_",                &epEp_ );
    EpTree_->Branch( "epRelErrorEp_",        &epRelErrorEp_ );
    EpTree_->Branch( "ecalDrivenEp_",        &ecalDrivenEp_ );
    EpTree_->Branch( "trackerDrivenSeedEp_", &trackerDrivenSeedEp_ );
    EpTree_->Branch( "classificationEp_",    &classificationEp_ );
    EpTree_->Branch( "isEBEp_",              &isEBEp_ );


    // =====================================
    // Making photon tree and setting branches

    photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");

    photonTree_->Branch( "rawEnergy_",        &ph_rawEnergy_        );
    photonTree_->Branch( "r9_",               &ph_r9_               );
    photonTree_->Branch( "etaWidth_",         &ph_etaWidth_         );
    photonTree_->Branch( "phiWidth_",         &ph_phiWidth_         );
    photonTree_->Branch( "numberOfClusters_", &ph_numberOfClusters_ );                     
    photonTree_->Branch( "hadronicOverEm_",   &ph_hadronicOverEm_   );
    photonTree_->Branch( "rhoValue_",         &ph_rhoValue_         );
    photonTree_->Branch( "delEtaSeed_",       &ph_delEtaSeed_       );
    photonTree_->Branch( "delPhiSeed_",       &ph_delPhiSeed_       );
    photonTree_->Branch( "seedEnergy_",       &ph_seedEnergy_       );
    photonTree_->Branch( "3x3_5x5_",          &ph_3x3_5x5_          );
    photonTree_->Branch( "sigmaIetaIeta_",    &ph_sigmaIetaIeta_    );
    photonTree_->Branch( "sigmaIphiIphi_",    &ph_sigmaIphiIphi_    );
    photonTree_->Branch( "sigmaIetaIphi_",    &ph_sigmaIetaIphi_    );
    photonTree_->Branch( "Emax_5x5_",         &ph_Emax_5x5_         );
    photonTree_->Branch( "e2nd_5x5_",         &ph_e2nd_5x5_         );
    photonTree_->Branch( "eTop_5x5_",         &ph_eTop_5x5_         );
    photonTree_->Branch( "eBottom_5x5_",      &ph_eBottom_5x5_      );
    photonTree_->Branch( "eLeft_5x5_",        &ph_eLeft_5x5_        );
    photonTree_->Branch( "eRight_5x5_",       &ph_eRight_5x5_       );
    photonTree_->Branch( "e2x5Max_5x5_",      &ph_e2x5Max_5x5_      );
    photonTree_->Branch( "e2x5Left_5x5_",     &ph_e2x5Left_5x5_     );
    photonTree_->Branch( "e2x5Right_5x5_",    &ph_e2x5Right_5x5_    );
    photonTree_->Branch( "e2x5Top_5x5_",      &ph_e2x5Top_5x5_      );
    photonTree_->Branch( "e2x5Bottom_5x5_",   &ph_e2x5Bottom_5x5_   );

    // Coordinate variables
    photonTree_->Branch( "isEB_",                   &ph_isEB_                   );
    photonTree_->Branch( "5x5_seedEnergy_",         &ph_5x5_seedEnergy_         );
    photonTree_->Branch( "iEtaCoordinate_",         &ph_iEtaCoordinate_         );
    photonTree_->Branch( "iPhiCoordinate_",         &ph_iPhiCoordinate_         );
    photonTree_->Branch( "iEtaMod5_",               &ph_iEtaMod5_               );
    photonTree_->Branch( "iPhiMod2_",               &ph_iPhiMod2_               );
    photonTree_->Branch( "iEtaMod20_",              &ph_iEtaMod20_              );
    photonTree_->Branch( "iPhiMod20_",              &ph_iPhiMod20_              );
    photonTree_->Branch( "preShowerE_rawEnergy_",   &ph_preShowerE_rawEnergy_   );
    photonTree_->Branch( "preShowerEp1_rawEnergy_", &ph_preShowerEp1_rawEnergy_ );
    photonTree_->Branch( "preShowerEp2_rawEnergy_", &ph_preShowerEp2_rawEnergy_ );
    photonTree_->Branch( "iXCoordinate_",           &ph_iXCoordinate_           );
    photonTree_->Branch( "iYCoordinate_",           &ph_iYCoordinate_           );

    // Matching variables

    photonTree_->Branch( "match_dR",     &match_dR     );
    photonTree_->Branch( "match_dE",     &match_dE     );
    photonTree_->Branch( "match_dRdE",   &match_dRdE   );

    photonTree_->Branch( "gen_pt",     &gen_pt     );
    photonTree_->Branch( "gen_phi",    &gen_phi    );
    photonTree_->Branch( "gen_eta",    &gen_eta    );
    photonTree_->Branch( "gen_M",      &gen_M      );
    photonTree_->Branch( "gen_E",      &gen_E      );
    photonTree_->Branch( "gen_pdgId",  &gen_pdgId  );
    photonTree_->Branch( "gen_status", &gen_status );

    }

// Deconstructor
SimpleNtuplizer::~SimpleNtuplizer() {
    //std::cout << ">>>> Inside SimpleNtuplizer::destructor" << std::endl;
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    }





//######################################
//# Member functions
//######################################

// =====================================
// analyze - The method that is executed on every event

void SimpleNtuplizer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ){
    
    using namespace std;
    using namespace edm;
    using namespace reco;

    //######################################
    //# Get all the collections
    //######################################

    // Get vertex collection
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);

    // Get electron collection
    edm::Handle<reco::GsfElectronCollection> electrons;  // For AODSIM
    iEvent.getByToken(electronToken_, electrons);

    // Get photon collection
    edm::Handle<reco::PhotonCollection> photons;
    iEvent.getByToken(photonToken_, photons);

    // Get GenParticle collection
    // Definition moved --> class variable
    // edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken( genParticleToken_, genParticles );


    //######################################
    //# Event specific quantities (not used in regression)
    //######################################

    // Determine number of primary vertices
    if (vertices->empty()) nPV_ = 0;
    else nPV_ = vertices->size();


    //######################################
    //# Analyze electrons
    //######################################

    // Loop over electrons
    nElectrons_ = 0;
    nElectronsMatched_ = 0;
    //for (const pat::Electron &el : *electrons) {
    for (const reco::GsfElectron &el : *electrons) {

        // Fill in the class variables for this particle; also sets E-p variables
        setElectronVariables( el, iEvent, iSetup );

        // Write E-p variables to the E-p tree
        // EpTree_->Fill();

        }

    // for (const reco::GenParticle &genParticle : *genParticles) {
    //     cout << "looping, genParticle pt = " << genParticle.pt() << endl;
    //     }


    // Attempt to match photons to genParticle
    // matchPhotonToGenParticle( photons, genParticles );

    // Loop over photons
    nPhotons_ = 0;
    nPhotonsMatched_ = 0;
    for (const reco::Photon &photon : *photons) {
        // Fill in the class variables for this particle; also sets E-p variables
        setPhotonVariables( photon, iEvent, iSetup );
        // Write E-p variables to the E-p tree
        //EpTree_->Fill();
        }


    // Fill in the event specific variables
    eventTree_->Fill();

    }






//######################################
//# Necessary functions and settings for framework
//######################################

// ------------ method called once each job just before starting event loop  ------------
void SimpleNtuplizer::beginJob() {
    //std::cout << ">>>> Inside SimpleNtuplizer::beginJob" << std::endl;
    }

// ------------ method called once each job just after ending the event loop  ------------
void SimpleNtuplizer::endJob() {
    //std::cout << ">>>> Inside SimpleNtuplizer::endJob" << std::endl;
    }

// // ------------ method called when starting to processes a run  ------------
// void SimpleNtuplizer::beginRun(edm::Run const&, edm::EventSetup const&){
//     }

// // ------------ method called when ending the processing of a run  ------------
// void SimpleNtuplizer::endRun(edm::Run const&, edm::EventSetup const&){
//     }

// // ------------ method called when starting to processes a luminosity block  ------------
// void SimpleNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
//     }

// // ------------ method called when ending the processing of a luminosity block  ------------
// void SimpleNtuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
//     }



//define this as a plug-in
DEFINE_FWK_MODULE(SimpleNtuplizer);
