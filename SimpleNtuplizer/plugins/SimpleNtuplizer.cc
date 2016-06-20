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
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
    ecalRecHitEBToken_(consumes<edm::SortedCollection<EcalRecHit>>(iConfig.getParameter<edm::InputTag>("ecalrechitsEB"))),
    ecalRecHitEEToken_(consumes<edm::SortedCollection<EcalRecHit>>(iConfig.getParameter<edm::InputTag>("ecalrechitsEE")))
    {

    std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

    edm::Service<TFileService> fs;


    //######################################
    //# eventTree
    //######################################

    // Event variables include the quantities of which there are exactly one per event

    eventTree_ = fs->make<TTree> ("EventTree", "Per event data");

    // Event variables
    eventTree_->Branch( "nPV",               &nPV_,                "nPV/I"   );
    eventTree_->Branch( "nElectrons",        &nElectrons_,         "nEle/I"  );
    eventTree_->Branch( "nElectronsMatched", &nElectronsMatched_ , "nEleMatched/I"  );
    eventTree_->Branch( "nPhotons",          &nPhotons_,           "nPho/I"  );
    eventTree_->Branch( "nPhotonsMatched",   &nPhotonsMatched_,    "nPhoMatched/I"  );

    eventTree_->Branch( "NtupID",                        &NtupID_               ); // For convenient cuts
    eventTree_->Branch( "eventNumber",                   &eventNumber_          );
    eventTree_->Branch( "luminosityBlock",               &luminosityBlock_      );
    eventTree_->Branch( "run",                           &run_                  );


    //######################################
    //# electronTree
    //######################################

    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    // Electron variables
    //  - All branch names now without a dash
    //  - Switched SC / SS tags to the front

    electronTree_->Branch( "NtupID",                        &NtupID_               ); // For convenient cuts
    electronTree_->Branch( "eventNumber",                   &eventNumber_          );
    electronTree_->Branch( "luminosityBlock",               &luminosityBlock_      );
    electronTree_->Branch( "run",                           &run_                  );

    electronTree_->Branch( "pt",                            &pt_e       );
    electronTree_->Branch( "scIsEB",                        &isEB_e     );
    electronTree_->Branch( "isMatched",                     &isMatched_ ); // Only matched is saved, so should be 1


    // =====================================
    // All actually used for training

    electronTree_->Branch( "nVtx",                          &nPV_                );
    electronTree_->Branch( "scRawEnergy",                   &rawEnergy_e         );
    electronTree_->Branch( "scEta",                         &eta_e               );
    electronTree_->Branch( "scPhi",                         &phi_e               );
    electronTree_->Branch( "scEtaWidth",                    &etaWidth_e          );
    electronTree_->Branch( "scPhiWidth",                    &phiWidth_e          );
    electronTree_->Branch( "scSeedR9",                      &r9_e                );
    electronTree_->Branch( "scSeedRawEnergy",               &seedEnergy_e        );
    electronTree_->Branch( "scSeedEmax",                    &eMax_e              );
    electronTree_->Branch( "scSeedE2nd",                    &e2nd_e              );
    electronTree_->Branch( "scSeedLeftRightAsym",           &eHorizontal_e       );
    electronTree_->Branch( "scSeedTopBottomAsym",           &eVertical_e         );
    electronTree_->Branch( "scSeedSigmaIetaIeta",           &sigmaIetaIeta_e     );
    electronTree_->Branch( "scSeedSigmaIetaIphi",           &sigmaIetaIphi_e     );
    electronTree_->Branch( "scSeedSigmaIphiIphi",           &sigmaIphiIphi_e     );
    electronTree_->Branch( "N_ECALClusters",                &numberOfClusters_e  );

    electronTree_->Branch( "N_SATURATEDXTALS",              &N_SATURATEDXTALS_e  );
    electronTree_->Branch( "seedIsSaturated",               &seedIsSaturated_e   );
    electronTree_->Branch( "seedCrystalEnergy",             &seedCrystalEnergy_e );


    // Max dR Cluster variables (now always 1 entry)
    electronTree_->Branch( "clusterMaxDR",                  &MaxDRclusterDR_e );
    electronTree_->Branch( "clusterMaxDRDPhi",              &MaxDRclusterDPhi_e );
    electronTree_->Branch( "clusterMaxDRDEta",              &MaxDRclusterDEta_e );
    electronTree_->Branch( "clusterMaxDRRawEnergy",         &MaxDRclusterRawEnergy_e );

    // Separate cluster variables; currently contain exactly 3 elements
    electronTree_->Branch( "clusterRawEnergy",              &clusterRawEnergy_e );        
    electronTree_->Branch( "clusterDPhiToSeed",             &clusterDPhiToSeed_e );         
    electronTree_->Branch( "clusterDEtaToSeed",             &clusterDEtaToSeed_e );

    // These for EB electrons
    electronTree_->Branch( "scSeedCryEta",                  &cryEtaCoordinate_e );         
    electronTree_->Branch( "scSeedCryPhi",                  &cryPhiCoordinate_e );         
    electronTree_->Branch( "scSeedCryIetaV2",               &iEtaCoordinate_e );
    electronTree_->Branch( "scSeedCryIphiV2",               &iPhiCoordinate_e );

    // These for EE electrons
    electronTree_->Branch( "scSeedCryIxV2",                 &iXCoordinate_e );         
    electronTree_->Branch( "scSeedCryIyV2",                 &iYCoordinate_e );         
    electronTree_->Branch( "scPreshowerEnergy",             &preshowerEnergy_e ); // <-- Also saved for EB


    // =====================================
    // Not used for training (but some still needed for e.g. cuts)

    electronTree_->Branch( "scSeedCryX",                    &cryXCoordinate_e );
    electronTree_->Branch( "scSeedCryY",                    &cryYCoordinate_e );
    
    electronTree_->Branch( "IsEcalEnergyCorrected",         &IsEcalEnergyCorrected_e );
    electronTree_->Branch( "CorrectedEcalEnergy",           &CorrectedEcalEnergy_e );
    electronTree_->Branch( "CorrectedEcalEnergyError",      &CorrectedEcalEnergyError_e );

    electronTree_->Branch( "genMatchdR",                    &genMatchdR_e     );
    electronTree_->Branch( "genMatchdE",                    &genMatchdE_e     );
    electronTree_->Branch( "genMatchdRdE",                  &genMatchdRdE_e   );
    electronTree_->Branch( "genPt",                         &genPt_e          );
    electronTree_->Branch( "genPhi",                        &genPhi_e         );
    electronTree_->Branch( "genEta",                        &genEta_e         );
    electronTree_->Branch( "genMass",                       &genMass_e        );
    electronTree_->Branch( "genEnergy",                     &genEnergy_e      );
    electronTree_->Branch( "genPdgId",                      &genPdgId_e       );
    electronTree_->Branch( "genStatus",                     &genStatus_e      );

    electronTree_->Branch( "trkMomentum",                   &trkMomentum_e );
    electronTree_->Branch( "trkMomentumError",              &trkMomentumError_e );
    electronTree_->Branch( "trkMomentumRelError",           &trkMomentumRelError_e );
    electronTree_->Branch( "eleEcalDriven",                 &ecalDriven_e );
    electronTree_->Branch( "eleTrackerDriven",              &trackerDriven_e );
    electronTree_->Branch( "eleClass",                      &classification_e );


    //######################################
    //# photonTree
    //######################################

    photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");

    photonTree_->Branch( "NtupID",                        &NtupID_               ); // For convenient cuts
    photonTree_->Branch( "eventNumber",                   &eventNumber_          );
    photonTree_->Branch( "luminosityBlock",               &luminosityBlock_      );
    photonTree_->Branch( "run",                           &run_                  );

    photonTree_->Branch( "pt",                            &pt_p       );
    photonTree_->Branch( "scIsEB",                        &isEB_p     );
    photonTree_->Branch( "isMatched",                     &isMatched_ ); // Only matched is saved, so should be 1


    // =====================================
    // All actually used for training

    photonTree_->Branch( "nVtx",                          &nPV_                );
    photonTree_->Branch( "scRawEnergy",                   &rawEnergy_p         );
    photonTree_->Branch( "scEta",                         &eta_p               );
    photonTree_->Branch( "scPhi",                         &phi_p               );
    photonTree_->Branch( "scEtaWidth",                    &etaWidth_p          );
    photonTree_->Branch( "scPhiWidth",                    &phiWidth_p          );
    photonTree_->Branch( "scSeedR9",                      &r9_p                );
    photonTree_->Branch( "scSeedRawEnergy",               &seedEnergy_p        );
    // photonTree_->Branch( "scSeedEmax",                    &eMax_p              );
    // photonTree_->Branch( "scSeedE2nd",                    &e2nd_p              );
    photonTree_->Branch( "scSeedLeftRightAsym",           &eHorizontal_p       );
    photonTree_->Branch( "scSeedTopBottomAsym",           &eVertical_p         );
    photonTree_->Branch( "scSeedSigmaIetaIeta",           &sigmaIetaIeta_p     );
    photonTree_->Branch( "scSeedSigmaIetaIphi",           &sigmaIetaIphi_p     );
    photonTree_->Branch( "scSeedSigmaIphiIphi",           &sigmaIphiIphi_p     );
    photonTree_->Branch( "N_ECALClusters",                &numberOfClusters_p  );

    // Extra variables not in electronTree
    photonTree_->Branch( "hadronicOverEm",                &hadronicOverEm_p    );
    photonTree_->Branch( "rhoValue",                      &rhoValue_p          );
    photonTree_->Branch( "delEtaSeed",                    &delEtaSeed_p        );
    photonTree_->Branch( "delPhiSeed",                    &delPhiSeed_p        );
    photonTree_->Branch( "e5x5",                          &e5x5_p              );
    photonTree_->Branch( "e3x3",                          &e3x3_p              );
    photonTree_->Branch( "eMax",                          &eMax_p              );
    photonTree_->Branch( "e2nd",                          &e2nd_p              );
    photonTree_->Branch( "eTop",                          &eTop_p              );
    photonTree_->Branch( "eBottom",                       &eBottom_p           );
    photonTree_->Branch( "eLeft",                         &eLeft_p             );
    photonTree_->Branch( "eRight",                        &eRight_p            );
    photonTree_->Branch( "e2x5Max",                       &e2x5Max_p           );
    photonTree_->Branch( "e2x5Left",                      &e2x5Left_p          );
    photonTree_->Branch( "e2x5Right",                     &e2x5Right_p         );
    photonTree_->Branch( "e2x5Top",                       &e2x5Top_p           );
    photonTree_->Branch( "e2x5Bottom",                    &e2x5Bottom_p        );
    photonTree_->Branch( "scPreshowerEnergy",             &preshowerEnergy_p   );

    photonTree_->Branch( "N_SATURATEDXTALS",              &N_SATURATEDXTALS_p  );
    photonTree_->Branch( "seedIsSaturated",               &seedIsSaturated_p   );
    photonTree_->Branch( "seedCrystalEnergy",             &seedCrystalEnergy_p );

    // Max dR Cluster variables (now always 1 entry)
    photonTree_->Branch( "clusterMaxDR",                  &MaxDRclusterDR_p );
    photonTree_->Branch( "clusterMaxDRDPhi",              &MaxDRclusterDPhi_p );
    photonTree_->Branch( "clusterMaxDRDEta",              &MaxDRclusterDEta_p );
    photonTree_->Branch( "clusterMaxDRRawEnergy",         &MaxDRclusterRawEnergy_p );

    // Separate cluster variables; currently contain exactly 3 elements
    photonTree_->Branch( "clusterRawEnergy",              &clusterRawEnergy_p );        
    photonTree_->Branch( "clusterDPhiToSeed",             &clusterDPhiToSeed_p );         
    photonTree_->Branch( "clusterDEtaToSeed",             &clusterDEtaToSeed_p );

    // EB coordinate variables
    photonTree_->Branch( "iEtaCoordinate",                &iEtaCoordinate_p         );
    photonTree_->Branch( "iPhiCoordinate",                &iPhiCoordinate_p         );
    photonTree_->Branch( "iEtaMod5",                      &iEtaMod5_p               );
    photonTree_->Branch( "iPhiMod2",                      &iPhiMod2_p               );
    photonTree_->Branch( "iEtaMod20",                     &iEtaMod20_p              );
    photonTree_->Branch( "iPhiMod20",                     &iPhiMod20_p              );
    // EE coordinate variables
    photonTree_->Branch( "preshowerEnergyPlane1",         &preshowerEnergyPlane1_p  );
    photonTree_->Branch( "preshowerEnergyPlane2",         &preshowerEnergyPlane2_p  );
    photonTree_->Branch( "iXCoordinate",                  &iXCoordinate_p           );
    photonTree_->Branch( "iYCoordinate",                  &iYCoordinate_p           );


    // 
    // Last minute additions: corrections
    photonTree_->Branch( "scEcalEnergy",              &scEcalEnergy_p            );
    photonTree_->Branch( "scEcalEnergyError",         &scEcalEnergyError_p       );
    photonTree_->Branch( "phoEcalEnergy",             &phoEcalEnergy_p           );
    photonTree_->Branch( "phoEcalEnergyError",        &phoEcalEnergyError_p      );
    photonTree_->Branch( "regression1Energy",         &regression1Energy_p       );
    photonTree_->Branch( "regression1EnergyError",    &regression1EnergyError_p  );
    photonTree_->Branch( "regression2Energy",         &regression2Energy_p       );
    photonTree_->Branch( "regression2EnergyError",    &regression2EnergyError_p  );

    // =====================================
    // Not used for training (but some still needed for e.g. cuts)

    photonTree_->Branch( "genMatchdR",                    &genMatchdR_p     );
    photonTree_->Branch( "genMatchdE",                    &genMatchdE_p     );
    photonTree_->Branch( "genMatchdRdE",                  &genMatchdRdE_p   );
    photonTree_->Branch( "genPt",                         &genPt_p          );
    photonTree_->Branch( "genPhi",                        &genPhi_p         );
    photonTree_->Branch( "genEta",                        &genEta_p         );
    photonTree_->Branch( "genMass",                       &genMass_p        );
    photonTree_->Branch( "genEnergy",                     &genEnergy_p      );
    photonTree_->Branch( "genPdgId",                      &genPdgId_p       );
    photonTree_->Branch( "genStatus",                     &genStatus_p      );

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
    //   Definition moved --> class variable
    //   edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken( genParticleToken_, genParticles_ );
    // iEvent.getByToken( caloClusterToken_, caloClusters_ );


    iEvent.getByToken( ecalRecHitEBToken_, ecalRecHitsEB_ );
    iEvent.getByToken( ecalRecHitEEToken_, ecalRecHitsEE_ );


    //######################################
    //# Event specific quantities (not used in regression)
    //######################################

    // Central event counter (specific to this output tree)
    NtupID_++;

    // Event specific variables
    eventNumber_     = iEvent.id().event();
    luminosityBlock_ = iEvent.id().luminosityBlock();
    run_             = iEvent.id().run();

    // Determine number of primary vertices
    if (vertices->empty()) nPV_ = 0;
    else nPV_ = vertices->size();


    //######################################
    //# Analyze electrons and photons
    //######################################    

    // Loop over electrons
    nElectrons_ = 0;
    nElectronsMatched_ = 0;
    for (const reco::GsfElectron &el : *electrons) {
        setElectronVariables( el, iEvent, iSetup );
        }

    // Loop over photons
    nPhotons_ = 0;
    nPhotonsMatched_ = 0;
    for (const reco::Photon &photon : *photons) {
        setPhotonVariables( photon, iEvent, iSetup );
        }

    // Fill in the event specific variables
    eventTree_->Fill();

    }


//######################################
//# Count saturated crystal; also gets the energy in the seed crystal
//######################################    

void SimpleNtuplizer::SetSaturationVariables( edm::Ptr<reco::CaloCluster> seedCluster, bool isEB, bool isElectron ){

    DetId seedId                                            = seedCluster->seed();
    std::vector< std::pair<DetId, float> > hitsAndFractions = seedCluster->hitsAndFractions();

    DetId hitId;        // ID of the hit in the seedCluster
    DetId ecalRecHitId; // ID of the hit in ecal -- to be compared with with hitId

    bool isSaturated;
    
    // Written-to-tree variables
    Int_t    N_SATURATEDXTALS  = 0;
    bool     seedIsSaturated   = false;
    Double_t seedCrystalEnergy = 0.0;

    // Get the right ecalRecHits collection (different for barrel and encap)
    edm::Handle<edm::SortedCollection<EcalRecHit>> ecalRecHits;
    if (isEB) ecalRecHits = ecalRecHitsEB_ ;
    else      ecalRecHits = ecalRecHitsEE_ ;

    // Loop over all hits in the seedCluster
    for (const std::pair<DetId, float> hitFractionPair : hitsAndFractions) {

        // Get hitId
        hitId    = std::get<0>(hitFractionPair);

        // Loop over all hits in ecal, find the hit that corresponds to the hit in the seedCluster
        for (const EcalRecHit &ecalRecHit : *ecalRecHits) {
            
            ecalRecHitId = ecalRecHit.detid();
            if (!( ecalRecHitId == hitId )) continue;

            isSaturated  = ecalRecHit.checkFlag( EcalRecHit::Flags::kSaturated );

            // Increase the count of saturated crystals
            if (isSaturated) N_SATURATEDXTALS++;

            // Check if this hit concerns the seed of the seedCluster
            if ( ecalRecHitId == seedId ){
                seedIsSaturated   = isSaturated;
                seedCrystalEnergy = ecalRecHit.energy();
                }
            }

        }


    if (seedCrystalEnergy == 0.0){
        std::cout << "  WARNING: Seed crystal energy is 0.0" << std::endl;
        }


    // Arrange output -- depends on type of particle
    if (isElectron){
        N_SATURATEDXTALS_e  = N_SATURATEDXTALS ;
        seedIsSaturated_e   = seedIsSaturated ;
        seedCrystalEnergy_e = seedCrystalEnergy ;
        }
    else {
        N_SATURATEDXTALS_p  = N_SATURATEDXTALS ;
        seedIsSaturated_p   = seedIsSaturated ;
        seedCrystalEnergy_p = seedCrystalEnergy ;
        }

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
