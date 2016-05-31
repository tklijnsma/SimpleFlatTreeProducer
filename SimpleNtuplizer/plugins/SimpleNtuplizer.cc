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

    // Central event counter (specific to this output tree)
    eventNumber = 0;

    std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

    edm::Service<TFileService> fs;


    // =====================================
    // Making event variable tree and setting branches
    // Event variables include the quantities of which there are exactly one per event

    eventTree_ = fs->make<TTree> ("EventTree", "Per event data");

    // Event variables
    eventTree_->Branch( "nPV",               &nPV_,                "nPV/I"   );
    eventTree_->Branch( "nElectrons",        &nElectrons_,         "nEle/I"  );
    eventTree_->Branch( "nElectronsMatched", &nElectronsMatched_ , "nEleMatched/I"  );
    eventTree_->Branch( "nPhotons",          &nPhotons_,           "nPho/I"  );
    eventTree_->Branch( "nPhotonsMatched",   &nPhotonsMatched_,    "nPhoMatched/I"  );


    // =====================================
    // Making electron tree and setting branches

    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    // Electron variables
    //  - All branch names now without a dash
    //  - Switched SC / SS tags to the front

    electronTree_->Branch( "nPV",                           &nPV_ );
    electronTree_->Branch( "eventNumber",                   &eventNumber );

    electronTree_->Branch( "pt",                            &pt_    );
    electronTree_->Branch( "SC_eta",                        &etaSC_ );
    electronTree_->Branch( "SC_phi",                        &phiSC_ );
    electronTree_->Branch( "SC_rawEnergy",                  &rawEnergySC_ );
    electronTree_->Branch( "SC_etaWidth",                   &etaWidthSC_ );
    electronTree_->Branch( "SC_phiWidth",                   &phiWidthSC_ );
    electronTree_->Branch( "SS_r9",                         &r9SS_ );
    electronTree_->Branch( "SS_seedEnergy_overRaw",         &seedEnergySS_ );
    electronTree_->Branch( "SS_eMax_overRaw",               &eMaxSS_ );
    electronTree_->Branch( "SS_e2nd_overRaw",               &e2ndSS_ );
    electronTree_->Branch( "SS_eHorizontal",                &eHorizontalSS_ );
    electronTree_->Branch( "SS_eVertical",                  &eVerticalSS_ );
    electronTree_->Branch( "SS_sigmaIetaIeta",              &sigmaIetaIetaSS_ );
    electronTree_->Branch( "SS_sigmaIetaIphi",              &sigmaIetaIphiSS_ );
    electronTree_->Branch( "SS_sigmaIphiIphi",              &sigmaIphiIphiSS_ );
    electronTree_->Branch( "preshowerEnergy_overRaw",       &preshowerEnergy_ );
    electronTree_->Branch( "SC_numberOfClustersSC",         &numberOfClustersSC_ );
    electronTree_->Branch( "isEB",                          &isEB_ );

    //electronTree_->Branch( "trkMomentum",                 &trkMomentum );
    //electronTree_->Branch( "trkMomentumError",            &trkMomentumError );

    // These for EB electrons
    electronTree_->Branch( "iEtaCoordinate",                &iEtaCoordinate_ );         
    electronTree_->Branch( "iPhiCoordinate",                &iPhiCoordinate_ );         
    electronTree_->Branch( "cryEtaCoordinate",              &cryEtaCoordinate_ );
    electronTree_->Branch( "cryPhiCoordinate",              &cryPhiCoordinate_ );

    // These for EE electrons
    electronTree_->Branch( "iXCoordinate",                  &iXCoordinate_ );         
    electronTree_->Branch( "iYCoordinate",                  &iYCoordinate_ );         
    electronTree_->Branch( "cryXCoordinate",                &cryXCoordinate_ );
    electronTree_->Branch( "cryYCoordinate",                &cryYCoordinate_ );

    // Max dR Cluster variables (now always 1 entry)
    electronTree_->Branch( "MaxDRclusterDR",                &MaxDRclusterDR );
    electronTree_->Branch( "MaxDRclusterDPhi",              &MaxDRclusterDPhi );
    electronTree_->Branch( "MaxDRclusterDEta",              &MaxDRclusterDEta );
    electronTree_->Branch( "MaxDRclusterRawEnergy_overRaw", &MaxDRclusterRawEnergy );

    // Separate cluster variables; currently contain exactly 3 elements
    electronTree_->Branch( "clusterRawEnergy_overRaw",      &clusterRawEnergy_ );        
    electronTree_->Branch( "clusterDPhiToSeed",             &clusterDPhiToSeed_ );         
    electronTree_->Branch( "clusterDEtaToSeed",             &clusterDEtaToSeed_ );

    electronTree_->Branch( "match_dR",                      &match_dR     );
    electronTree_->Branch( "match_dE",                      &match_dE     );
    electronTree_->Branch( "match_dRdE",                    &match_dRdE   );

    electronTree_->Branch( "gen_pt",                        &gen_pt     );
    electronTree_->Branch( "gen_phi",                       &gen_phi    );
    electronTree_->Branch( "gen_eta",                       &gen_eta    );
    electronTree_->Branch( "gen_M",                         &gen_M      );
    electronTree_->Branch( "gen_E",                         &gen_E      );
    electronTree_->Branch( "gen_pdgId",                     &gen_pdgId  );
    electronTree_->Branch( "gen_status",                    &gen_status );


    // =====================================
    // Making E-p tree and setting branches

    // EpTree_ = fs->make<TTree> ("EpTree", "E-p data");
    // --> Just put this in the Electron tree as well

    electronTree_->Branch( "EP_totEnergy",            &totEnergyEp_ );
    electronTree_->Branch( "EP_trkMomentum",          &epEp_ );
    electronTree_->Branch( "EP_trkMomentumError",     &epErrorEp_ );
    electronTree_->Branch( "EP_trkMomentumRelError",  &epRelErrorEp_ );
    electronTree_->Branch( "EP_ecalDriven",           &ecalDrivenEp_ );
    electronTree_->Branch( "EP_trackerDrivenSeed",    &trackerDrivenSeedEp_ );
    electronTree_->Branch( "EP_classification",       &classificationEp_ );
    electronTree_->Branch( "EP_isEB",                 &isEBEp_ );


    // =====================================
    // Making photon tree and setting branches

    photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");

    photonTree_->Branch( "rawEnergy",        &ph_rawEnergy_        );
    photonTree_->Branch( "r9",               &ph_r9_               );
    photonTree_->Branch( "etaWidth",         &ph_etaWidth_         );
    photonTree_->Branch( "phiWidth",         &ph_phiWidth_         );
    photonTree_->Branch( "numberOfClusters", &ph_numberOfClusters_ );                     
    photonTree_->Branch( "hadronicOverEm",   &ph_hadronicOverEm_   );
    photonTree_->Branch( "rhoValue",         &ph_rhoValue_         );
    photonTree_->Branch( "delEtaSeed",       &ph_delEtaSeed_       );
    photonTree_->Branch( "delPhiSeed",       &ph_delPhiSeed_       );
    photonTree_->Branch( "seedEnergy",       &ph_seedEnergy_       );
    photonTree_->Branch( "3x3_5x5",          &ph_3x3_5x5_          );
    photonTree_->Branch( "sigmaIetaIeta",    &ph_sigmaIetaIeta_    );
    photonTree_->Branch( "sigmaIphiIphi",    &ph_sigmaIphiIphi_    );
    photonTree_->Branch( "sigmaIetaIphi",    &ph_sigmaIetaIphi_    );
    photonTree_->Branch( "Emax_5x5",         &ph_Emax_5x5_         );
    photonTree_->Branch( "e2nd_5x5",         &ph_e2nd_5x5_         );
    photonTree_->Branch( "eTop_5x5",         &ph_eTop_5x5_         );
    photonTree_->Branch( "eBottom_5x5",      &ph_eBottom_5x5_      );
    photonTree_->Branch( "eLeft_5x5",        &ph_eLeft_5x5_        );
    photonTree_->Branch( "eRight_5x5",       &ph_eRight_5x5_       );
    photonTree_->Branch( "e2x5Max_5x5",      &ph_e2x5Max_5x5_      );
    photonTree_->Branch( "e2x5Left_5x5",     &ph_e2x5Left_5x5_     );
    photonTree_->Branch( "e2x5Right_5x5",    &ph_e2x5Right_5x5_    );
    photonTree_->Branch( "e2x5Top_5x5",      &ph_e2x5Top_5x5_      );
    photonTree_->Branch( "e2x5Bottom_5x5",   &ph_e2x5Bottom_5x5_   );

    // Coordinate variables
    photonTree_->Branch( "isEB",                   &ph_isEB_                   );
    photonTree_->Branch( "5x5_seedEnergy",         &ph_5x5_seedEnergy_         );
    photonTree_->Branch( "iEtaCoordinate",         &ph_iEtaCoordinate_         );
    photonTree_->Branch( "iPhiCoordinate",         &ph_iPhiCoordinate_         );
    photonTree_->Branch( "iEtaMod5",               &ph_iEtaMod5_               );
    photonTree_->Branch( "iPhiMod2",               &ph_iPhiMod2_               );
    photonTree_->Branch( "iEtaMod20",              &ph_iEtaMod20_              );
    photonTree_->Branch( "iPhiMod20",              &ph_iPhiMod20_              );
    photonTree_->Branch( "preShowerE_rawEnergy",   &ph_preShowerE_rawEnergy_   );
    photonTree_->Branch( "preShowerEp1_rawEnergy", &ph_preShowerEp1_rawEnergy_ );
    photonTree_->Branch( "preShowerEp2_rawEnergy", &ph_preShowerEp2_rawEnergy_ );
    photonTree_->Branch( "iXCoordinate",           &ph_iXCoordinate_           );
    photonTree_->Branch( "iYCoordinate",           &ph_iYCoordinate_           );

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

    // Central event counter (specific to this output tree)
    eventNumber++;

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
