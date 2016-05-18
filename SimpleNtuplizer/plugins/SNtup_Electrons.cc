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
    bool successful_match = matchElectronToGenParticle( electron );
    if(!successful_match) return;

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
    pt_              = electron.pt() ;
    rawEnergySC_     = rawEnergy ;
    etaSC_           = superCluster->eta() ;
    phiSC_           = superCluster->phi() ;
    etaWidthSC_      = superCluster->etaWidth() ;
    phiWidthSC_      = superCluster->phiWidth() ;
    r9SS_            = electron.showerShape().r9 ;
    seedEnergySS_    = superCluster->seed()->energy() / rawEnergy ;
    eMaxSS_          = electron.showerShape().eMax / rawEnergy ;
    e2ndSS_          = electron.showerShape().e2nd / rawEnergy ;
    eHorizontalSS_   = electron.showerShape().eLeft + electron.showerShape().eRight != 0.f  
                                 ? ( electron.showerShape().eLeft - electron.showerShape().eRight ) /
                                   ( electron.showerShape().eLeft + electron.showerShape().eRight ) : 0.f  ;
    eVerticalSS_     = electron.showerShape().eTop + electron.showerShape().eBottom != 0.f
                                 ? ( electron.showerShape().eTop - electron.showerShape().eBottom ) /
                                   ( electron.showerShape().eTop + electron.showerShape().eBottom ) : 0.f  ;
    sigmaIetaIetaSS_ = electron.showerShape().sigmaIetaIeta ;
    sigmaIetaIphiSS_ = electron.showerShape().sigmaIetaIphi ;
    sigmaIphiIphiSS_ = electron.showerShape().sigmaIphiIphi ;
    preshowerEnergy_ = superCluster->preshowerEnergy() / superCluster->rawEnergy() ;        

    numberOfClusters = std::max( 0, int (superCluster->clusters().size()) );
    numberOfClustersSC_ = numberOfClusters ;
    isEB_               = electron.isEB() ;


    // =====================================
    // Cluster variables (subs of the superCluster)

    // Clear the std::vectors from the previous electron
    clusterRawEnergy_.clear();
    clusterDPhiToSeed_.clear();
    clusterDEtaToSeed_.clear();

    // These have either 0 or 1 entry
    MaxDRclusterDR_.clear();
    MaxDRclusterDPhi_.clear();
    MaxDRclusterDEta_.clear();
    MaxDRclusterRawEnergy_.clear();

    // Default values
    float MaxDRclusterDR        = 999.;
    float MaxDRclusterDPhi      = 999.;
    float MaxDRclusterDEta      = 999.;
    float MaxDRclusterRawEnergy = 0.;
    float maxDR                 = 0.;

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
        clusterRawEnergy_  .push_back( cluster->energy() / rawEnergy );
        clusterDPhiToSeed_ .push_back( reco::deltaPhi( cluster->phi(), superCluster->seed()->phi() ) );
        clusterDEtaToSeed_ .push_back( cluster->eta() - superCluster->seed()->eta() );

        // Find the cluster that has maximum delR to the seed
        const auto deltaR = reco::deltaR( *cluster, *superCluster->seed() );
        if( deltaR > maxDR) {
            maxDR = deltaR;
            MaxDRclusterDR        = maxDR;
            MaxDRclusterDPhi      = clusterDPhiToSeed_[i_cluster];
            MaxDRclusterDEta      = clusterDEtaToSeed_[i_cluster];
            MaxDRclusterRawEnergy = clusterRawEnergy_[i_cluster];
            }

        i_cluster++;

        // If cutting off after a certain amount of clusters, set this limit here
        if(i_cluster == 3) continue;
        }

    // Write cluster variables to class member vectors, only if needed
    if( i_cluster > 0 ) {
        // Still use vectors, since these vars may be empty
        MaxDRclusterDR_        .push_back( MaxDRclusterDR );    
        MaxDRclusterDPhi_      .push_back( MaxDRclusterDPhi );    
        MaxDRclusterDEta_      .push_back( MaxDRclusterDEta );    
        MaxDRclusterRawEnergy_ .push_back( MaxDRclusterRawEnergy );
        }


    // =====================================
    // Coordinate variables
    // Does different things for when the electron is in the barrel or endcap

    // Open up temporary variables
    int iPhi, iEta, iX, iY; float cryPhi, cryEta, cryX, cryY, dummy;
    EcalClusterLocal ecalLocal;

    // Clear up values from previous electron
    iEtaCoordinate_.clear();
    iPhiCoordinate_.clear();
    cryEtaCoordinate_.clear();
    cryPhiCoordinate_.clear();

    iXCoordinate_.clear();
    iYCoordinate_.clear();
    cryXCoordinate_.clear();
    cryYCoordinate_.clear();

    if( electron.isEB() ){
        ecalLocal.localCoordsEB( *superCluster->seed(), iSetup,
                                  cryEta, cryPhi, iEta, iPhi, dummy, dummy );
        iEtaCoordinate_   .push_back( iEta );
        iPhiCoordinate_   .push_back( iPhi );
        cryEtaCoordinate_ .push_back( cryEta );
        cryPhiCoordinate_ .push_back( cryPhi );
        }

    else{
        ecalLocal.localCoordsEE( *superCluster->seed(), iSetup,
                                  cryX, cryY, iX, iY, dummy, dummy );
        iXCoordinate_     .push_back( iX );
        iYCoordinate_     .push_back( iY );
        cryXCoordinate_   .push_back( cryX );
        cryYCoordinate_   .push_back( cryY );
        }


    //######################################
    //# Analyze EP
    //######################################

    totEnergyEp_         = superCluster->rawEnergy() + superCluster->preshowerEnergy();
    epEp_                = electron.trackMomentumAtVtx().R();
    epRelErrorEp_        = electron.trackMomentumError() / electron.trackMomentumAtVtx().R();
    ecalDrivenEp_        = electron.ecalDriven();
    trackerDrivenSeedEp_ = electron.trackerDrivenSeed();
    classificationEp_    = int(electron.classification());
    isEBEp_              = electron.isEB();



    // Write class variables to the output EpTree_
    electronTree_->Fill();

    // Write E-p variables to the E-p tree
    EpTree_->Fill();

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

    for (const reco::GenParticle &genParticle : *genParticles) {

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

    match_dR     = minDr;
    match_dE     = minDe;
    match_dRdE   = minDeDr;
    gen_pt     = matched_genParticle->pt();
    gen_phi    = matched_genParticle->phi();
    gen_eta    = matched_genParticle->eta();
    gen_M      = matched_genParticle->mass();
    gen_E      = matched_genParticle->energy();
    gen_pdgId  = matched_genParticle->pdgId();
    gen_status = matched_genParticle->status();

    // Return successful match value (should be true)
    return successful_match;

    }

