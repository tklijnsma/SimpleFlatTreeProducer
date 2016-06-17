#include "SimpleNtuplizer.h"


// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setElectronVariables(
        const reco::GsfElectron& electron,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

    using namespace std;
    using namespace edm;
    using namespace reco;


    // =====================================
    // Gen matching

    // Increase count of photons in event
    nElectrons_++;

    // Try to match to genParticle; quit function if electron is not matched
    if (!isData_) {
      bool successful_match = matchElectronToGenParticle( electron );
      if(!successful_match) return;
    }
    
    // Increase count of matched photons in event
    nElectronsMatched_++;


    // =====================================
    // Fill other variables

    //cout << "Setting class variables for type " << typeid(electron).name() << endl;

    // Convenience definitions
    double rawEnergy;       // Raw energy of the super cluster per electron
    int numberOfClusters;   // Number of (sub)clusters in the super cluster

    // Open variable superCluster
    const reco::SuperClusterRef& superCluster = electron.superCluster();

    // Raw energy of the super cluster per electron
    rawEnergy = superCluster->rawEnergy();

    // Write electron variables to class variables
    pt_e               = electron.pt() ;
    rawEnergy_e        = rawEnergy ;
    eta_e              = superCluster->eta() ;
    phi_e              = superCluster->phi() ;
    etaWidth_e         = superCluster->etaWidth() ;
    phiWidth_e         = superCluster->phiWidth() ;
    r9_e               = electron.showerShape().r9 ;
    seedEnergy_e       = superCluster->seed()->energy() ;
    eMax_e             = electron.showerShape().eMax ;
    e2nd_e             = electron.showerShape().e2nd ;
    eHorizontal_e      = electron.showerShape().eLeft + electron.showerShape().eRight != 0.f  
                                    ? ( electron.showerShape().eLeft - electron.showerShape().eRight ) /
                                      ( electron.showerShape().eLeft + electron.showerShape().eRight ) : 0.f  ;
    eVertical_e        = electron.showerShape().eTop + electron.showerShape().eBottom != 0.f
                                    ? ( electron.showerShape().eTop - electron.showerShape().eBottom ) /
                                      ( electron.showerShape().eTop + electron.showerShape().eBottom ) : 0.f  ;
    sigmaIetaIeta_e    = electron.showerShape().sigmaIetaIeta ;
    sigmaIetaIphi_e    = electron.showerShape().sigmaIetaIphi ;
    sigmaIphiIphi_e    = electron.showerShape().sigmaIphiIphi ;
    preshowerEnergy_e  = superCluster->preshowerEnergy() ;        

    numberOfClusters   = std::max( 0, int (superCluster->clusters().size()) );
    numberOfClusters_e = numberOfClusters ;
    isEB_e             = electron.isEB() ;

    // To compare to current 'official' regression
    IsEcalEnergyCorrected_e    = electron.corrections().isEcalEnergyCorrected;
    CorrectedEcalEnergy_e      = electron.corrections().correctedEcalEnergy;
    CorrectedEcalEnergyError_e = electron.corrections().correctedEcalEnergyError;

    // =====================================
    // Cluster variables (subs of the superCluster)

    // Clear the std::vectors from the previous electron
    clusterRawEnergy_e.clear();
    clusterDPhiToSeed_e.clear();
    clusterDEtaToSeed_e.clear();

    // Fill with zeroes for the needed amount of clusters
    //   make sure at least indices 0, 1 and 2 are filled
    clusterRawEnergy_e.resize(std::max(3, numberOfClusters), 0);
    clusterDPhiToSeed_e.resize(std::max(3, numberOfClusters), 0);
    clusterDEtaToSeed_e.resize(std::max(3, numberOfClusters), 0);

    // Default values
    MaxDRclusterDR_e        = 0.; // Changed from 999.; Make sure it's considered
    MaxDRclusterDPhi_e      = 0.; // Changed from 999.; Make sure it's considered
    MaxDRclusterDEta_e      = 0.; // Changed from 999.; Make sure it's considered
    MaxDRclusterRawEnergy_e = 0.;

    // compared throughout loop
    float maxDR             = 0.;

    // Define cluster
    edm::Ptr<reco::CaloCluster> cluster;

    // Loop over all clusters that aren't the seed  
    size_t i_cluster = 0;
    for( auto pcluster = superCluster->clustersBegin(); pcluster != superCluster->clustersEnd(); ++pcluster ) {
        
        // Dereference to get the actual object
        cluster = *pcluster;

        // Continue if this is the seed
        if( cluster == superCluster->seed() ) continue;

        // Set basic cluster quantities in vectors
        // clusterRawEnergy_  .push_back( cluster->energy() / rawEnergy );
        // clusterDPhiToSeed_ .push_back( reco::deltaPhi( cluster->phi(), superCluster->seed()->phi() ) );
        // clusterDEtaToSeed_ .push_back( cluster->eta() - superCluster->seed()->eta() );
        clusterRawEnergy_e[i_cluster]  = cluster->energy() ;
        clusterDPhiToSeed_e[i_cluster] = reco::deltaPhi( cluster->phi(), superCluster->seed()->phi() );
        clusterDEtaToSeed_e[i_cluster] = cluster->eta() - superCluster->seed()->eta();

        // Find the cluster that has maximum delR to the seed
        const auto deltaR = reco::deltaR( *cluster, *superCluster->seed() );
        if( deltaR > maxDR) {
            maxDR = deltaR;
            MaxDRclusterDR_e        = maxDR;
            MaxDRclusterDPhi_e      = clusterDPhiToSeed_e[i_cluster];
            MaxDRclusterDEta_e      = clusterDEtaToSeed_e[i_cluster];
            MaxDRclusterRawEnergy_e = clusterRawEnergy_e[i_cluster];
            }

        i_cluster++;

        // If cutting off after a certain amount of clusters, set this limit here
        if(i_cluster == 3) break;
        }


    // =====================================
    // Coordinate variables
    // Does different things for when the electron is in the barrel or endcap

    // Open up temporary variables
    int iPhi, iEta, iX, iY; float cryPhi, cryEta, cryX, cryY, dummy;
    EcalClusterLocal ecalLocal;

    // Clear up values from previous electron
    iEtaCoordinate_e.clear();
    iPhiCoordinate_e.clear();
    cryEtaCoordinate_e.clear();
    cryPhiCoordinate_e.clear();

    iXCoordinate_e.clear();
    iYCoordinate_e.clear();
    cryXCoordinate_e.clear();
    cryYCoordinate_e.clear();

    if( electron.isEB() ){
        ecalLocal.localCoordsEB( *superCluster->seed(), iSetup,
                                  cryEta, cryPhi, iEta, iPhi, dummy, dummy );
        iEtaCoordinate_e   .push_back( iEta );
        iPhiCoordinate_e   .push_back( iPhi );
        cryEtaCoordinate_e .push_back( cryEta );
        cryPhiCoordinate_e .push_back( cryPhi );
        }

    else{
        ecalLocal.localCoordsEE( *superCluster->seed(), iSetup,
                                  cryX, cryY, iX, iY, dummy, dummy );
        iXCoordinate_e     .push_back( iX );
        iYCoordinate_e     .push_back( iY );
        cryXCoordinate_e   .push_back( cryX );
        cryYCoordinate_e   .push_back( cryY );
        }


    //######################################
    //# Analyze EP
    //######################################

    trkMomentum_e          = electron.trackMomentumAtVtx().R();
    trkMomentumError_e     = electron.trackMomentumError();
    trkMomentumRelError_e  = electron.trackMomentumError() / electron.trackMomentumAtVtx().R();
    ecalDriven_e           = electron.ecalDriven();
    trackerDriven_e    = electron.trackerDrivenSeed();
    classification_e       = int(electron.classification());

    // Write class variables to the output EpTree_
    electronTree_->Fill();

    }


bool SimpleNtuplizer::matchElectronToGenParticle(
        // const reco::PhotonCollection& photons,
        // const reco::GenParticleCollection& genPartices
        //const edm::Handle<reco::PhotonCollection> photons,
        const reco::GsfElectron& electron
        // const edm::Handle<reco::GenParticleCollection> genParticles
        ){

    // int nTempCounter = 0;
    // std::cout << "In matchElectronToGenParticle" << std::endl;
    // for (const reco::GenParticle &genParticle : *genParticles) {
    //     nTempCounter++;
    //     std::cout << "    genParticle " << nTempCounter << std::endl;
    //     std::cout << "        pt     = " << genParticle.pt() << std::endl;
    //     std::cout << "        eta    = " << genParticle.eta() << std::endl;
    //     std::cout << "        phi    = " << genParticle.phi() << std::endl;
    //     std::cout << "        pdgId  = " << genParticle.pdgId() << std::endl;
    //     std::cout << "        status = " << genParticle.status() << std::endl;
    //     }

    //######################################
    //# Start matching
    //######################################

    // =====================================
    // Setting local matching variables

    // Maximum match radius
    double match_MaxDR = 0.5;

    // Keep track of minimum dX's
    double minDr   = 1e6;
    double minDe   = 1e6;
    double minDeDr = 1e6;

    // dX's of the match between current electron and genParticle
    double this_dr;
    double this_de;
    double this_dedr;

    // Only use the electron if it's matched successfully
    bool successful_match = false;
    const reco::GenParticle* matched_genParticle;


    // =====================================
    // Loop over genParticles

    for (const reco::GenParticle &genParticle : *genParticles_) {

        // Continue if pdgId is not 22 or status is not 1
        if(!( abs(genParticle.pdgId())==11 && genParticle.status()==1 ))
            continue;

        // Calculate distance variables
        this_dr   = reco::deltaR( genParticle, electron );
        this_de   = fabs( genParticle.energy()- electron.energy() ) / genParticle.energy();
        this_dedr = sqrt( this_dr*this_dr + this_de*this_de );

        if( this_dr < match_MaxDR
            // && this_dr<minDr       // matching type 1
            // && this_de<minDe       // matching type 2
            && this_dedr < minDeDr    // matching type 3
            ){

            minDr   = this_dr;
            minDe   = this_de;
            minDeDr = this_dedr;

            successful_match = true;
            matched_genParticle = &genParticle;

            }
        }

    // std::cout << "        minDr   = " << minDr << std::endl;
    // std::cout << "        minDe   = " << minDe << std::endl;
    // std::cout << "        minDeDr = " << minDeDr << std::endl;

    // Return if particle could not be matched
    if(!successful_match) return successful_match;


    // =====================================
    // Fill necessary branches

    genMatchdR_e    = minDr;
    genMatchdE_e    = minDe;
    genMatchdRdE_e  = minDeDr;
    genPt_e         = matched_genParticle->pt();
    genPhi_e        = matched_genParticle->phi();
    genEta_e        = matched_genParticle->eta();
    genMass_e       = matched_genParticle->mass();
    genEnergy_e     = matched_genParticle->energy();
    genPdgId_e      = matched_genParticle->pdgId();
    genStatus_e     = matched_genParticle->status();

    // Return successful match value (should be true)
    return successful_match;

    }
