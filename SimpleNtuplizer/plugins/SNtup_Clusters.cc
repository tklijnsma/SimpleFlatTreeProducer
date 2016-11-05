#include "SimpleNtuplizer.h"


// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setClusterVariables(
        const reco::CaloCluster& cluster,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

    using namespace std;
    using namespace edm;
    using namespace reco;


    // =====================================
    // Gen matching

    // Increase count of photons in event
    nClusters_++;

    // Try to match to genParticle; quit function if cluster is not matched
    if (!isData_) {
      bool successful_match = matchCaloClusterToGenParticle( cluster );
      if(!successful_match) return;
    }
    
    // Increase count of matched photons in event
    nClustersMatched_++;


    // =====================================
    // Fill other variables

    //cout << "Setting class variables for type " << typeid(cluster).name() << endl;

    // Convenience definitions
    double rawEnergy;       // Raw energy of the super cluster per cluster

    // Raw energy of the super cluster per electron
    rawEnergy = cluster.energy();

    // Get rho
    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken_,rhoH);
    Float_t rho = *rhoH;


    //######################################
    //# Start filling branch variables
    //######################################

    // Check if it's barrel or not
    isEB_c             = (cluster.seed().subdetId()==EcalBarrel) ;

    // Get the right ecalRecHits collection (different for barrel and encap)
    edm::Handle<edm::SortedCollection<EcalRecHit>> ecalRecHits;
    if (isEB_c) ecalRecHits = ecalRecHitsEB_ ;
    else        ecalRecHits = ecalRecHitsEE_ ;


    // Write cluster variables to class variables
    rawEnergy_c        = rawEnergy ;
    eta_c              = cluster.eta() ;
    phi_c              = cluster.phi() ;
    pt_c               = rawEnergy/cosh(eta_c) ;

    rhoValue_c         = rho;

    // =====================================
    // Showershape variables

    //    const auto& showerShape = cluster.showerShape();

    std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances( cluster, &*ecalRecHits, topology_ );
    sigmaIetaIeta_c = sqrt(localCovariances[0]);
    sigmaIetaIphi_c = localCovariances[1];
    sigmaIphiIphi_c = sqrt(localCovariances[2]);
    
    eMax_c             = EcalClusterToolsT<false>::eMax( cluster, &*ecalRecHits );
    e2nd_c             = EcalClusterToolsT<false>::e2nd( cluster, &*ecalRecHits );
    eTop_c             = EcalClusterToolsT<false>::eTop( cluster, &*ecalRecHits, topology_ );
    eBottom_c          = EcalClusterToolsT<false>::eBottom( cluster, &*ecalRecHits, topology_ );
    eLeft_c            = EcalClusterToolsT<false>::eLeft( cluster, &*ecalRecHits, topology_ );
    eRight_c           = EcalClusterToolsT<false>::eRight( cluster, &*ecalRecHits, topology_ );
    eHorizontal_c      = (eLeft_c + eRight_c) != 0.f  
                                    ? ( eLeft_c - eRight_c ) /
                                      ( eLeft_c + eRight_c ) : 0.f  ;
    eVertical_c        = eTop_c + eBottom_c != 0.f
                                    ? ( eTop_c - eBottom_c ) /
                                      ( eTop_c + eBottom_c ) : 0.f  ;

    e5x5_c             = EcalClusterToolsT<false>::e5x5( cluster, &*ecalRecHits, topology_ );
    e2x5Max_c          = EcalClusterToolsT<false>::e2x5Max( cluster, &*ecalRecHits, topology_ );
    
    // EcalClusterToolsT<noZS>; noZS = full5x5, ZS = weighted
    e3x3_c             = EcalClusterToolsT<false>::e3x3( cluster, &*ecalRecHits, topology_ );
    r9_c               = e3x3_c/rawEnergy_c;
    e2x5Left_c         = EcalClusterToolsT<false>::e2x5Left(   cluster, &*ecalRecHits, topology_ );
    e2x5Right_c        = EcalClusterToolsT<false>::e2x5Right(  cluster, &*ecalRecHits, topology_ );
    e2x5Top_c          = EcalClusterToolsT<false>::e2x5Top(    cluster, &*ecalRecHits, topology_ );
    e2x5Bottom_c       = EcalClusterToolsT<false>::e2x5Bottom( cluster, &*ecalRecHits, topology_ );
    

    // -------------------------------
    // Repeat for the full 5x5

    std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances( cluster, &*ecalRecHits, topology_ );
    full5x5_sigmaIetaIeta_c = sqrt(localCovariances[0]);
    full5x5_sigmaIetaIphi_c = localCovariances[1];
    full5x5_sigmaIphiIphi_c = sqrt(localCovariances[2]);
    
    full5x5_eMax_c             = EcalClusterToolsT<true>::eMax( cluster, &*ecalRecHits );
    full5x5_e2nd_c             = EcalClusterToolsT<true>::e2nd( cluster, &*ecalRecHits );
    full5x5_eTop_c             = EcalClusterToolsT<true>::eTop( cluster, &*ecalRecHits, topology_ );
    full5x5_eBottom_c          = EcalClusterToolsT<true>::eBottom( cluster, &*ecalRecHits, topology_ );
    full5x5_eLeft_c            = EcalClusterToolsT<true>::eLeft( cluster, &*ecalRecHits, topology_ );
    full5x5_eRight_c           = EcalClusterToolsT<true>::eRight( cluster, &*ecalRecHits, topology_ );
    full5x5_eHorizontal_c      = (eLeft_c + eRight_c) != 0.f  
      ? ( eLeft_c - eRight_c ) /
      ( eLeft_c + eRight_c ) : 0.f  ;
    full5x5_eVertical_c        = eTop_c + eBottom_c != 0.f
                                    ? ( eTop_c - eBottom_c ) /
      ( eTop_c + eBottom_c ) : 0.f  ;
    
    full5x5_e5x5_c             = EcalClusterToolsT<true>::e5x5( cluster, &*ecalRecHits, topology_ );
    full5x5_e2x5Max_c          = EcalClusterToolsT<true>::e2x5Max( cluster, &*ecalRecHits, topology_ );
    
    // EcalClusterToolsT<noZS>; noZS = full5x5, ZS = weighted
    full5x5_e3x3_c             = EcalClusterToolsT<true>::e3x3( cluster, &*ecalRecHits, topology_ );
    full5x5_r9_c               = e3x3_c/rawEnergy_c;
    full5x5_e2x5Left_c         = EcalClusterToolsT<true>::e2x5Left(   cluster, &*ecalRecHits, topology_ );
    full5x5_e2x5Right_c        = EcalClusterToolsT<true>::e2x5Right(  cluster, &*ecalRecHits, topology_ );
    full5x5_e2x5Top_c          = EcalClusterToolsT<true>::e2x5Top(    cluster, &*ecalRecHits, topology_ );
    full5x5_e2x5Bottom_c       = EcalClusterToolsT<true>::e2x5Bottom( cluster, &*ecalRecHits, topology_ );

    // =====================================
    // Saturation variables

    SetClusterSaturationVariables( cluster, ecalRecHits );


    // =====================================
    // Coordinate variables
    // Does different things for when the cluster is in the barrel or endcap

    // Open up temporary variables
    int iPhi, iEta, iX, iY; float cryPhi, cryEta, cryX, cryY, dummy;
    EcalClusterLocal ecalLocal;


    // Clear the std::vectors from last cluster
    // EB
    cryEtaCoordinate_c       .clear();
    cryPhiCoordinate_c       .clear();
    iEtaCoordinate_c         .clear() ;
    iPhiCoordinate_c         .clear() ;
    iEtaMod5_c               .clear() ;
    iPhiMod2_c               .clear() ;
    iEtaMod20_c              .clear() ;
    iPhiMod20_c              .clear() ;


    // EE
    cryXCoordinate_c         .clear();
    cryYCoordinate_c         .clear();
    iXCoordinate_c           .clear() ;
    iYCoordinate_c           .clear() ;


    if( isEB_c ){
        
        ecalLocal.localCoordsEB(  cluster, iSetup,
                                  cryEta, cryPhi, iEta, iPhi, dummy, dummy );

        iEtaCoordinate_c   .push_back( iEta );
        iPhiCoordinate_c   .push_back( iPhi );
        cryEtaCoordinate_c .push_back( cryEta );
        cryPhiCoordinate_c .push_back( cryPhi );


        int signiEta = iEta > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training

        iEtaMod5_c       .push_back( (iEta-signiEta)%5 );
        iPhiMod2_c       .push_back( (iPhi-1)%2 );
        iEtaMod20_c      .push_back( abs(iEta)<=25 ? iEta-signiEta : (iEta - 26*signiEta) % 20  );
        iPhiMod20_c      .push_back( (iPhi-1)%20 );

        }


    else{
        
        ecalLocal.localCoordsEE(  cluster, iSetup,
                                  cryX, cryY, iX, iY, dummy, dummy );
        
        iXCoordinate_c     .push_back( iX );
        iYCoordinate_c     .push_back( iY );
        cryXCoordinate_c   .push_back( cryX );
        cryYCoordinate_c   .push_back( cryY );


        }



    // Write class variables to the output EpTree_
    clusterTree_->Fill();

    }


bool SimpleNtuplizer::matchCaloClusterToGenParticle(
        const reco::CaloCluster& cluster
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

    // dX's of the match between current cluster and genParticle
    double this_dr;
    double this_de;
    double this_dedr;

    // Only use the cluster if it's matched successfully
    bool successful_match = false;
    const reco::GenParticle* matched_genParticle;

    // =====================================
    // Loop over genParticles

    for (const reco::GenParticle &genParticle : *genParticles_) {

        // Continue if pdgId is not 11 or status is not 1
        if(!( abs(genParticle.pdgId())==11 && genParticle.status()==1 ))
            continue;

        // Calculate distance variables
        this_dr   = reco::deltaR( genParticle, cluster.position() );
        this_de   = fabs( genParticle.energy()- cluster.energy() ) / genParticle.energy();
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

    genMatchdR_c    = minDr;
    genMatchdE_c    = minDe;
    genMatchdRdE_c  = minDeDr;
    genPt_c         = matched_genParticle->pt();
    genPhi_c        = matched_genParticle->phi();
    genEta_c        = matched_genParticle->eta();
    genMass_c       = matched_genParticle->mass();
    genEnergy_c     = matched_genParticle->energy();
    genPdgId_c      = matched_genParticle->pdgId();
    genStatus_c     = matched_genParticle->status();

    // Return successful match value (should be true)
    return successful_match;

    }
