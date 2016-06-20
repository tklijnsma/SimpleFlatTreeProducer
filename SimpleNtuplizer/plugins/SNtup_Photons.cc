#include "SimpleNtuplizer.h"

// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setPhotonVariables(
        const reco::Photon& photon,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

    // =====================================
    // Gen matching

    // Increase count of photons in event
    nPhotons_++;

    // Try to match to genParticle; quit function if photon is not matched
    bool successful_match = matchPhotonToGenParticle( photon );
    if(!successful_match) return;

    // Increase count of matched photons in event
    nPhotonsMatched_++;


    // =====================================
    // Fill main variables

    // Set convenience variables
    const reco::SuperClusterRef& superCluster = photon.superCluster();
    const edm::Ptr<reco::CaloCluster>& seedCluster = superCluster->seed();
    int numberOfClusters;

    // T: Not sure if this should be done!!
    // const bool missing_clusters = !superCluster->clusters()[numberOfClusters-1].isAvailable();
    // if( missing_clusters ) return ; // do not apply corrections in case of missing info (slimmed MiniAOD electrons)

    const double rawEnergy = superCluster->rawEnergy(); 
    const auto& showerShape = photon.showerShapeVariables();

    // Get rho
    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken_,rhoH);
    Float_t rho = *rhoH;


    // Write photon variables to class variables
    pt_p              = photon.pt() ;
    rawEnergy_p     = rawEnergy ;
    eta_p           = superCluster->eta() ;
    phi_p           = superCluster->phi() ;
    etaWidth_p      = superCluster->etaWidth() ;
    phiWidth_p      = superCluster->phiWidth() ;
    r9_p            = photon.r9() ;
    seedEnergy_p    = superCluster->seed()->energy() ;

    // Vars NOT in electronTree:
    hadronicOverEm_p   = photon.hadronicOverEm();
    rhoValue_p         = rho;
    delEtaSeed_p       = seedCluster->eta()-superCluster->position().Eta();
    delPhiSeed_p       = reco::deltaPhi(seedCluster->phi(),superCluster->position().Phi());
    e5x5_p             = showerShape.e5x5       ;
    e3x3_p             = showerShape.e3x3       ;
    eMax_p             = showerShape.maxEnergyXtal ; // <-- Also in electronTree
    e2nd_p             = showerShape.e2nd       ; // <-- Also in electronTree
    eTop_p             = showerShape.eTop       ;
    eBottom_p          = showerShape.eBottom    ;
    eLeft_p            = showerShape.eLeft      ;
    eRight_p           = showerShape.eRight     ;  
    e2x5Max_p          = showerShape.e2x5Max    ;
    e2x5Left_p         = showerShape.e2x5Left   ;
    e2x5Right_p        = showerShape.e2x5Right  ;
    e2x5Top_p          = showerShape.e2x5Top    ;
    e2x5Bottom_p       = showerShape.e2x5Bottom ;

    eHorizontal_p   = showerShape.eLeft + showerShape.eRight != 0.f  
                                 ? ( showerShape.eLeft - showerShape.eRight ) /
                                   ( showerShape.eLeft + showerShape.eRight ) : 0.f  ;
    eVertical_p     = showerShape.eTop + showerShape.eBottom != 0.f
                                 ? ( showerShape.eTop - showerShape.eBottom ) /
                                   ( showerShape.eTop + showerShape.eBottom ) : 0.f  ;
    sigmaIetaIeta_p = showerShape.sigmaIetaIeta ;
    sigmaIetaIphi_p = showerShape.sigmaIetaIphi ;
    sigmaIphiIphi_p = showerShape.sigmaIphiIphi ;
    preshowerEnergy_p = superCluster->preshowerEnergy() ;        

    numberOfClusters   = std::max( 0, int (superCluster->clusters().size()) );
    numberOfClusters_p = numberOfClusters ;
    isEB_p              = photon.isEB() ;


    // Which one is the current 74 regression energy?
    scEcalEnergy_p            = photon.energyCorrections().scEcalEnergy ;
    scEcalEnergyError_p       = photon.energyCorrections().scEcalEnergyError ;
    phoEcalEnergy_p           = photon.energyCorrections().phoEcalEnergy ;
    phoEcalEnergyError_p      = photon.energyCorrections().phoEcalEnergyError ;
    regression1Energy_p       = photon.energyCorrections().regression1Energy ;
    regression1EnergyError_p  = photon.energyCorrections().regression1EnergyError ;
    regression2Energy_p       = photon.energyCorrections().regression2Energy ;
    regression2EnergyError_p  = photon.energyCorrections().regression2EnergyError ;

    // std::cout << "Some vars for the photon:"                      << std::endl;
    // std::cout << "rawEnergy_p         = " <<  rawEnergy_p         << std::endl ;    
    // std::cout << "scEcalEnergy_p      = " <<  scEcalEnergy_p      << std::endl ;
    // std::cout << "phoEcalEnergy_p     = " <<  phoEcalEnergy_p     << std::endl ;
    // std::cout << "regression1Energy_p = " <<  regression1Energy_p << std::endl ;
    // std::cout << "regression2Energy_p = " <<  regression2Energy_p << std::endl ;

    // std::cout << "scEcalEnergyError_p      = " <<  scEcalEnergyError_p      << std::endl ;
    // std::cout << "phoEcalEnergyError_p     = " <<  phoEcalEnergyError_p     << std::endl ;
    // std::cout << "regression1EnergyError_p = " <<  regression1EnergyError_p << std::endl ;
    // std::cout << "regression2EnergyError_p = " <<  regression2EnergyError_p << std::endl ;

    // Double_t the_energy = photon.energy();
    // std::cout << "photon.energy()     = " <<  the_energy << std::endl << std::endl;
    

    // Saturation variables
    SetSaturationVariables( superCluster->seed(), isEB_p, false );


    // =====================================
    // Cluster variables

    // Clear the std::vectors from the previous electron
    clusterRawEnergy_p.clear();
    clusterDPhiToSeed_p.clear();
    clusterDEtaToSeed_p.clear();

    // Fill with zeroes for the needed amount of clusters
    //   make sure at least indices 0, 1 and 2 are filled
    clusterRawEnergy_p.resize(std::max(3, numberOfClusters), 0);
    clusterDPhiToSeed_p.resize(std::max(3, numberOfClusters), 0);
    clusterDEtaToSeed_p.resize(std::max(3, numberOfClusters), 0);

    // These have either 0 or 1 entry
    // MaxDRclusterDR_.clear();
    // MaxDRclusterDPhi_.clear();
    // MaxDRclusterDEta_.clear();
    // MaxDRclusterRawEnergy_.clear();

    // Default values
    MaxDRclusterDR_p        = 0.; // Changed from 999.; Make sure it's considered
    MaxDRclusterDPhi_p      = 0.; // Changed from 999.; Make sure it's considered
    MaxDRclusterDEta_p      = 0.; // Changed from 999.; Make sure it's considered
    MaxDRclusterRawEnergy_p = 0.;

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
        clusterRawEnergy_p[i_cluster]  = cluster->energy() ;
        clusterDPhiToSeed_p[i_cluster] = reco::deltaPhi( cluster->phi(), superCluster->seed()->phi() );
        clusterDEtaToSeed_p[i_cluster] = cluster->eta() - superCluster->seed()->eta();

        // Find the cluster that has maximum delR to the seed
        const auto deltaR = reco::deltaR( *cluster, *superCluster->seed() );
        if( deltaR > maxDR) {
            maxDR = deltaR;
            MaxDRclusterDR_p        = maxDR;
            MaxDRclusterDPhi_p      = clusterDPhiToSeed_p[i_cluster];
            MaxDRclusterDEta_p      = clusterDEtaToSeed_p[i_cluster];
            MaxDRclusterRawEnergy_p = clusterRawEnergy_p[i_cluster];
            }

        i_cluster++;

        // If cutting off after a certain amount of clusters, set this limit here
        if(i_cluster == 3) break;
        }


    // =====================================
    // Coordinate variables

    // Clear the std::vectors from last photon
    // EB
    iEtaCoordinate_p         .clear() ;
    iPhiCoordinate_p         .clear() ;
    iEtaMod5_p               .clear() ;
    iPhiMod2_p               .clear() ;
    iEtaMod20_p              .clear() ;
    iPhiMod20_p              .clear() ;
    // EE
    preshowerEnergyPlane1_p  .clear() ;
    preshowerEnergyPlane2_p  .clear() ;
    iXCoordinate_p           .clear() ;
    iYCoordinate_p           .clear() ;

    if ( photon.isEB() ) {
        
        EBDetId ebseedid( seedCluster->seed());
        
        int ieta = ebseedid.ieta();
        int iphi = ebseedid.iphi();
        int signieta = ieta > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training

        iEtaCoordinate_p .push_back( ieta );
        iPhiCoordinate_p .push_back( iphi );
        iEtaMod5_p       .push_back( (ieta-signieta)%5 );
        iPhiMod2_p       .push_back( (iphi-1)%2 );
        iEtaMod20_p      .push_back( abs(ieta)<=25 ? ieta-signieta : (ieta - 26*signieta) % 20  );
        iPhiMod20_p      .push_back( (iphi-1)%20 );
        }
    
    else {
        EEDetId eeseedid( seedCluster->seed());
        preshowerEnergyPlane1_p  .push_back( superCluster->preshowerEnergyPlane1() / rawEnergy );
        preshowerEnergyPlane2_p  .push_back( superCluster->preshowerEnergyPlane2() / rawEnergy );
        iXCoordinate_p           .push_back( eeseedid.ix() );
        iYCoordinate_p           .push_back( eeseedid.iy() );
        }

    // Write class variables to the output tree
    photonTree_->Fill();

    }





bool SimpleNtuplizer::matchPhotonToGenParticle(
        // const reco::PhotonCollection& photons,
        // const reco::GenParticleCollection& genPartices
        //const edm::Handle<reco::PhotonCollection> photons,

        const reco::Photon& photon
        // const edm::Handle<reco::GenParticleCollection> genParticles
        ){

    // int nTempCounter = 0;
    // std::cout << "In matchPhotonToGenParticle" << std::endl;
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

    // dX's of the match between current photon and genParticle
    double this_dr;
    double this_de;
    double this_dedr;

    // Only use the photon if it's matched successfully
    bool successful_match = false;
    const reco::GenParticle* matched_genParticle;


    // =====================================
    // Loop over genParticles

    for (const reco::GenParticle &genParticle : *genParticles_) {

        // Continue if pdgId is not 22 or status is not 1
        if(!( abs(genParticle.pdgId())==22 && genParticle.status()==1 ))
            continue;

        // Calculate distance variables
        this_dr   = reco::deltaR( genParticle, photon );
        this_de   = fabs( genParticle.energy()- photon.energy() ) / genParticle.energy();
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

    genMatchdR_p   = minDr;
    genMatchdE_p   = minDe;
    genMatchdRdE_p = minDeDr;
    genPt_p        = matched_genParticle->pt();
    genPhi_p       = matched_genParticle->phi();
    genEta_p       = matched_genParticle->eta();
    genMass_p      = matched_genParticle->mass();
    genEnergy_p    = matched_genParticle->energy();
    genPdgId_p     = matched_genParticle->pdgId();
    genStatus_p    = matched_genParticle->status();

    // Return successful match value (should be true)
    return successful_match;

    }

