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
    if (!isData_) {
      bool successful_match = matchPhotonToGenParticle( photon );
      if(!successful_match) return;
    }
    
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

    // Get rho
    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken_,rhoH);
    Float_t rho = *rhoH;


    //######################################
    //# Start filling branch variables
    //######################################

    // Check if it's barrel or not
    isEB_p             = photon.isEB() ;

    // Write photon variables to class variables
    pt_p               = photon.pt() ;
    rawEnergy_p        = rawEnergy ;
    eta_p              = superCluster->eta() ;
    phi_p              = superCluster->phi() ;
    etaWidth_p         = superCluster->etaWidth() ;
    phiWidth_p         = superCluster->phiWidth() ;

    seedEnergy_p       = superCluster->seed()->energy() ;

    // Vars NOT in electronTree:
    hadronicOverEm_p   = photon.hadronicOverEm();
    rhoValue_p         = rho;
    delEtaSeed_p       = seedCluster->eta()-superCluster->position().Eta();
    delPhiSeed_p       = reco::deltaPhi(seedCluster->phi(),superCluster->position().Phi());

    preshowerEnergy_p = superCluster->preshowerEnergy() ;        

    numberOfClusters   = std::max( 0, int (superCluster->clusters().size()) );
    numberOfClusters_p = numberOfClusters ;


    // =====================================
    // Showershape variables

    const auto& showerShape = photon.showerShapeVariables();
    
    r9_p               = photon.r9() ;
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

    eHorizontal_p      = showerShape.eLeft + showerShape.eRight != 0.f  
                                 ? ( showerShape.eLeft - showerShape.eRight ) /
                                   ( showerShape.eLeft + showerShape.eRight ) : 0.f  ;
    eVertical_p        = showerShape.eTop + showerShape.eBottom != 0.f
                                 ? ( showerShape.eTop - showerShape.eBottom ) /
                                   ( showerShape.eTop + showerShape.eBottom ) : 0.f  ;
    sigmaIetaIeta_p    = showerShape.sigmaIetaIeta ;
    sigmaIetaIphi_p    = showerShape.sigmaIetaIphi ;
    sigmaIphiIphi_p    = showerShape.sigmaIphiIphi ;

    // -------------------------------
    // Repeat for the full 5x5

    const auto& full5x5_showerShape = photon.full5x5_showerShapeVariables();
    
    full5x5_r9_p               = photon.full5x5_r9() ;
    full5x5_e5x5_p             = full5x5_showerShape.e5x5       ;
    full5x5_e3x3_p             = full5x5_showerShape.e3x3       ;
    full5x5_eMax_p             = full5x5_showerShape.maxEnergyXtal ; // <-- Also in electronTree
    full5x5_e2nd_p             = full5x5_showerShape.e2nd       ; // <-- Also in electronTree
    full5x5_eTop_p             = full5x5_showerShape.eTop       ;
    full5x5_eBottom_p          = full5x5_showerShape.eBottom    ;
    full5x5_eLeft_p            = full5x5_showerShape.eLeft      ;
    full5x5_eRight_p           = full5x5_showerShape.eRight     ;  
    full5x5_e2x5Max_p          = full5x5_showerShape.e2x5Max    ;
    full5x5_e2x5Left_p         = full5x5_showerShape.e2x5Left   ;
    full5x5_e2x5Right_p        = full5x5_showerShape.e2x5Right  ;
    full5x5_e2x5Top_p          = full5x5_showerShape.e2x5Top    ;
    full5x5_e2x5Bottom_p       = full5x5_showerShape.e2x5Bottom ;

    full5x5_eHorizontal_p      = full5x5_showerShape.eLeft + full5x5_showerShape.eRight != 0.f  
                                 ? ( full5x5_showerShape.eLeft - full5x5_showerShape.eRight ) /
                                   ( full5x5_showerShape.eLeft + full5x5_showerShape.eRight ) : 0.f  ;
    full5x5_eVertical_p        = full5x5_showerShape.eTop + full5x5_showerShape.eBottom != 0.f
                                 ? ( full5x5_showerShape.eTop - full5x5_showerShape.eBottom ) /
                                   ( full5x5_showerShape.eTop + full5x5_showerShape.eBottom ) : 0.f  ;
    full5x5_sigmaIetaIeta_p    = full5x5_showerShape.sigmaIetaIeta ;
    full5x5_sigmaIetaIphi_p    = full5x5_showerShape.sigmaIetaIphi ;
    full5x5_sigmaIphiIphi_p    = full5x5_showerShape.sigmaIphiIphi ;


    // Which one is the current 74 regression energy?
    // scEcalEnergy_p            = photon.energyCorrections().scEcalEnergy ;
    // scEcalEnergyError_p       = photon.energyCorrections().scEcalEnergyError ;
    // phoEcalEnergy_p           = photon.energyCorrections().phoEcalEnergy ;
    // phoEcalEnergyError_p      = photon.energyCorrections().phoEcalEnergyError ;
    // regression1Energy_p       = photon.energyCorrections().regression1Energy ;
    // regression1EnergyError_p  = photon.energyCorrections().regression1EnergyError ;
    // regression2Energy_p       = photon.energyCorrections().regression2Energy ;
    // regression2EnergyError_p  = photon.energyCorrections().regression2EnergyError ;
    corrEnergy74X_p       = photon.energy();
    corrEnergy74XError_p  = photon.energyCorrections().regression1EnergyError ;


    // Saturation variables
    edm::Handle<edm::SortedCollection<EcalRecHit>> ecalRecHits;
    if (isEB_p) ecalRecHits = ecalRecHitsEB_ ;
    else        ecalRecHits = ecalRecHitsEE_ ;
    SetSaturationVariables( superCluster->seed(), ecalRecHits, false );


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

    // // Clear the std::vectors from last photon
    // // EB
    // iEtaCoordinate_p         .clear() ;
    // iPhiCoordinate_p         .clear() ;
    // iEtaMod5_p               .clear() ;
    // iPhiMod2_p               .clear() ;
    // iEtaMod20_p              .clear() ;
    // iPhiMod20_p              .clear() ;
    // // EE
    // preshowerEnergyPlane1_p  .clear() ;
    // preshowerEnergyPlane2_p  .clear() ;
    // iXCoordinate_p           .clear() ;
    // iYCoordinate_p           .clear() ;

    // if ( photon.isEB() ) {
        
    //     EBDetId ebseedid( seedCluster->seed());
        
    //     int ieta = ebseedid.ieta();
    //     int iphi = ebseedid.iphi();
    //     int signieta = ieta > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training

    //     iEtaCoordinate_p .push_back( ieta );
    //     iPhiCoordinate_p .push_back( iphi );
    //     iEtaMod5_p       .push_back( (ieta-signieta)%5 );
    //     iPhiMod2_p       .push_back( (iphi-1)%2 );
    //     iEtaMod20_p      .push_back( abs(ieta)<=25 ? ieta-signieta : (ieta - 26*signieta) % 20  );
    //     iPhiMod20_p      .push_back( (iphi-1)%20 );
    //     }
    
    // else {
        
    //     EEDetId eeseedid( seedCluster->seed());

    //     preshowerEnergyPlane1_p  .push_back( superCluster->preshowerEnergyPlane1() / rawEnergy );
    //     preshowerEnergyPlane2_p  .push_back( superCluster->preshowerEnergyPlane2() / rawEnergy );
    //     iXCoordinate_p           .push_back( eeseedid.ix() );
    //     iYCoordinate_p           .push_back( eeseedid.iy() );
    //     }


    // Open up temporary variables
    int iPhi, iEta, iX, iY; float cryPhi, cryEta, cryX, cryY, dummy;
    EcalClusterLocal ecalLocal;

    // Clear the std::vectors from last electron
    // EB
    cryEtaCoordinate_p       .clear();
    cryPhiCoordinate_p       .clear();
    iEtaCoordinate_p         .clear() ;
    iPhiCoordinate_p         .clear() ;
    iEtaMod5_p               .clear() ;
    iPhiMod2_p               .clear() ;
    iEtaMod20_p              .clear() ;
    iPhiMod20_p              .clear() ;


    // EE
    cryXCoordinate_p         .clear();
    cryYCoordinate_p         .clear();
    iXCoordinate_p           .clear() ;
    iYCoordinate_p           .clear() ;
    preshowerEnergyPlane1_p  .clear() ;
    preshowerEnergyPlane2_p  .clear() ;



    if( photon.isEB() ){
        
        ecalLocal.localCoordsEB( *superCluster->seed(), iSetup,
                                  cryEta, cryPhi, iEta, iPhi, dummy, dummy );

        iEtaCoordinate_p   .push_back( iEta );
        iPhiCoordinate_p   .push_back( iPhi );
        cryEtaCoordinate_p .push_back( cryEta );
        cryPhiCoordinate_p .push_back( cryPhi );


        int signiEta = iEta > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training

        iEtaMod5_p       .push_back( (iEta-signiEta)%5 );
        iPhiMod2_p       .push_back( (iPhi-1)%2 );
        iEtaMod20_p      .push_back( abs(iEta)<=25 ? iEta-signiEta : (iEta - 26*signiEta) % 20  );
        iPhiMod20_p      .push_back( (iPhi-1)%20 );

        }


    else{
        
        ecalLocal.localCoordsEE( *superCluster->seed(), iSetup,
                                  cryX, cryY, iX, iY, dummy, dummy );
        
        iXCoordinate_p     .push_back( iX );
        iYCoordinate_p     .push_back( iY );
        cryXCoordinate_p   .push_back( cryX );
        cryYCoordinate_p   .push_back( cryY );

        preshowerEnergyPlane1_p  .push_back( superCluster->preshowerEnergyPlane1() / rawEnergy );
        preshowerEnergyPlane2_p  .push_back( superCluster->preshowerEnergyPlane2() / rawEnergy );

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

