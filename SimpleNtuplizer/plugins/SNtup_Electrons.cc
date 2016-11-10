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

    // Get rho
    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken_,rhoH);
    Float_t rho = *rhoH;
    
    //######################################
    //# Start filling branch variables
    //######################################

    // Check if it's barrel or not
    isEB_e             = electron.isEB() ;

    // Get the right ecalRecHits collection (different for barrel and encap)
    edm::Handle<edm::SortedCollection<EcalRecHit>> ecalRecHits;
    if (isEB_e) ecalRecHits = ecalRecHitsEB_ ;
    else        ecalRecHits = ecalRecHitsEE_ ;


    // Write electron variables to class variables
    pt_e               = electron.pt() ;
    rawEnergy_e        = rawEnergy ;
    eta_e              = superCluster->eta() ;
    phi_e              = superCluster->phi() ;
    etaWidth_e         = superCluster->etaWidth() ;
    phiWidth_e         = superCluster->phiWidth() ;
    seedEnergy_e       = superCluster->seed()->energy() ;

    preshowerEnergy_e  = superCluster->preshowerEnergy() ;

    numberOfClusters   = std::max( 0, int (superCluster->clusters().size()) );
    numberOfClusters_e = numberOfClusters ;

    hadronicOverEm_e   = electron.hcalDepth1OverEcalBc() + electron.hcalDepth2OverEcalBc();
    hadronic1OverEm_e  = electron.hcalDepth1OverEcalBc();
    hadronic2OverEm_e  = electron.hcalDepth2OverEcalBc();
    rhoValue_e         = rho;
    delEtaSeed_e       = superCluster->seed()->eta() - superCluster->position().Eta();
    delPhiSeed_e       = reco::deltaPhi( superCluster->seed()->phi(),superCluster->position().Phi());


    // To compare to current 'official' regression
    // IsEcalEnergyCorrected_e    = electron.corrections().isEcalEnergyCorrected;
    // CorrectedEcalEnergy_e      = electron.corrections().correctedEcalEnergy;
    // CorrectedEcalEnergyError_e = electron.corrections().correctedEcalEnergyError;
    corrEnergy74X_e       = electron.corrections().correctedEcalEnergy;
    corrEnergy74XError_e  = electron.corrections().correctedEcalEnergyError;


    // =====================================
    // Showershape variables

    const auto& showerShape = electron.showerShape();

    r9_e               = showerShape.r9 ;
    eMax_e             = showerShape.eMax ;
    e2nd_e             = showerShape.e2nd ;
    eHorizontal_e      = showerShape.eLeft + showerShape.eRight != 0.f  
                                    ? ( showerShape.eLeft - showerShape.eRight ) /
                                      ( showerShape.eLeft + showerShape.eRight ) : 0.f  ;
    eVertical_e        = showerShape.eTop + showerShape.eBottom != 0.f
                                    ? ( showerShape.eTop - showerShape.eBottom ) /
                                      ( showerShape.eTop + showerShape.eBottom ) : 0.f  ;
    sigmaIetaIeta_e    = showerShape.sigmaIetaIeta ;
    sigmaIetaIphi_e    = showerShape.sigmaIetaIphi ;
    sigmaIphiIphi_e    = showerShape.sigmaIphiIphi ;

    e5x5_e             = showerShape.e5x5       ;
    eTop_e             = showerShape.eTop       ;
    eBottom_e          = showerShape.eBottom    ;
    eLeft_e            = showerShape.eLeft      ;
    eRight_e           = showerShape.eRight     ;  
    e2x5Max_e          = showerShape.e2x5Max    ;
    
    // EcalClusterToolsT<noZS>; noZS = full5x5, ZS = weighted
    e3x3_e             = EcalClusterToolsT<false>::e3x3(       *superCluster->seed(), &*ecalRecHits, topology_ );
    e2x5Left_e         = EcalClusterToolsT<false>::e2x5Left(   *superCluster->seed(), &*ecalRecHits, topology_ );
    e2x5Right_e        = EcalClusterToolsT<false>::e2x5Right(  *superCluster->seed(), &*ecalRecHits, topology_ );
    e2x5Top_e          = EcalClusterToolsT<false>::e2x5Top(    *superCluster->seed(), &*ecalRecHits, topology_ );
    e2x5Bottom_e       = EcalClusterToolsT<false>::e2x5Bottom( *superCluster->seed(), &*ecalRecHits, topology_ );


    // -------------------------------
    // Repeat for the full 5x5

    const auto& full5x5_showerShape = electron.full5x5_showerShape();

    full5x5_r9_e               = full5x5_showerShape.r9 ;
    full5x5_eMax_e             = full5x5_showerShape.eMax ;
    full5x5_e2nd_e             = full5x5_showerShape.e2nd ;
    full5x5_eHorizontal_e      = full5x5_showerShape.eLeft + full5x5_showerShape.eRight != 0.f  
                                    ? ( full5x5_showerShape.eLeft - full5x5_showerShape.eRight ) /
                                      ( full5x5_showerShape.eLeft + full5x5_showerShape.eRight ) : 0.f  ;
    full5x5_eVertical_e        = full5x5_showerShape.eTop + full5x5_showerShape.eBottom != 0.f
                                    ? ( full5x5_showerShape.eTop - full5x5_showerShape.eBottom ) /
                                      ( full5x5_showerShape.eTop + full5x5_showerShape.eBottom ) : 0.f  ;
    full5x5_sigmaIetaIeta_e    = full5x5_showerShape.sigmaIetaIeta ;
    full5x5_sigmaIetaIphi_e    = full5x5_showerShape.sigmaIetaIphi ;
    full5x5_sigmaIphiIphi_e    = full5x5_showerShape.sigmaIphiIphi ;

    full5x5_e5x5_e             = full5x5_showerShape.e5x5       ;
    full5x5_eTop_e             = full5x5_showerShape.eTop       ;
    full5x5_eBottom_e          = full5x5_showerShape.eBottom    ;
    full5x5_eLeft_e            = full5x5_showerShape.eLeft      ;
    full5x5_eRight_e           = full5x5_showerShape.eRight     ;  
    full5x5_e2x5Max_e          = full5x5_showerShape.e2x5Max    ;
    
    // EcalClusterToolsT<noZS>; noZS = full5x5, ZS = weighted
    full5x5_e3x3_e             = EcalClusterToolsT<true>::e3x3(       *superCluster->seed(), &*ecalRecHits, topology_ );
    full5x5_e2x5Left_e         = EcalClusterToolsT<true>::e2x5Left(   *superCluster->seed(), &*ecalRecHits, topology_ );
    full5x5_e2x5Right_e        = EcalClusterToolsT<true>::e2x5Right(  *superCluster->seed(), &*ecalRecHits, topology_ );
    full5x5_e2x5Top_e          = EcalClusterToolsT<true>::e2x5Top(    *superCluster->seed(), &*ecalRecHits, topology_ );
    full5x5_e2x5Bottom_e       = EcalClusterToolsT<true>::e2x5Bottom( *superCluster->seed(), &*ecalRecHits, topology_ );


    // =====================================
    // Saturation variables

    SetSaturationVariables( superCluster->seed(), ecalRecHits, true );


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


    // Clear the std::vectors from last electron
    // EB
    cryEtaCoordinate_e       .clear();
    cryPhiCoordinate_e       .clear();
    iEtaCoordinate_e         .clear() ;
    iPhiCoordinate_e         .clear() ;
    iEtaMod5_e               .clear() ;
    iPhiMod2_e               .clear() ;
    iEtaMod20_e              .clear() ;
    iPhiMod20_e              .clear() ;


    // EE
    cryXCoordinate_e         .clear();
    cryYCoordinate_e         .clear();
    iXCoordinate_e           .clear() ;
    iYCoordinate_e           .clear() ;
    preshowerEnergyPlane1_e  .clear() ;
    preshowerEnergyPlane2_e  .clear() ;



    if( electron.isEB() ){
        
        ecalLocal.localCoordsEB( *superCluster->seed(), iSetup,
                                  cryEta, cryPhi, iEta, iPhi, dummy, dummy );

        iEtaCoordinate_e   .push_back( iEta );
        iPhiCoordinate_e   .push_back( iPhi );
        cryEtaCoordinate_e .push_back( cryEta );
        cryPhiCoordinate_e .push_back( cryPhi );


        int signiEta = iEta > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training

        iEtaMod5_e       .push_back( (iEta-signiEta)%5 );
        iPhiMod2_e       .push_back( (iPhi-1)%2 );
        iEtaMod20_e      .push_back( abs(iEta)<=25 ? iEta-signiEta : (iEta - 26*signiEta) % 20  );
        iPhiMod20_e      .push_back( (iPhi-1)%20 );

        }


    else{
        
        ecalLocal.localCoordsEE( *superCluster->seed(), iSetup,
                                  cryX, cryY, iX, iY, dummy, dummy );
        
        iXCoordinate_e     .push_back( iX );
        iYCoordinate_e     .push_back( iY );
        cryXCoordinate_e   .push_back( cryX );
        cryYCoordinate_e   .push_back( cryY );

        preshowerEnergyPlane1_e  .push_back( superCluster->preshowerEnergyPlane1() / rawEnergy );
        preshowerEnergyPlane2_e  .push_back( superCluster->preshowerEnergyPlane2() / rawEnergy );

        }


    //######################################
    //# Analyze EP
    //######################################
    
    auto el_track          = electron.gsfTrack();

    trkMomentum_e          = el_track->pMode();
    trkEta_e               = el_track->etaMode();
    trkPhi_e               = el_track->phiMode();

    float ptMode       = el_track->ptMode();
    float ptModeErrror = el_track->ptModeError();    
    float etaModeError = el_track->etaModeError();
    // p = pT cosh(eta) -> dp = sqrt ( pTerr^2 * cosh^2(eta) + pT^2 * sinh^2(eta) * etaerr^2)
    float pModeError   = sqrt(ptModeErrror*ptModeErrror*cosh(trkEta_e)*cosh(trkEta_e) + ptMode*ptMode*sinh(trkEta_e)*sinh(trkEta_e)*etaModeError*etaModeError);

    trkMomentumError_e     = pModeError;
    trkMomentumRelError_e  = trkMomentumError_e/trkMomentum_e;

    gsfchi2_e              = el_track->chi2();
    gsfndof_e              = el_track->ndof();
    gsfnhits_e             = el_track->numberOfValidHits();

    ecalDriven_e           = electron.ecalDriven();
    trackerDriven_e        = electron.trackerDrivenSeed();
    classification_e       = int(electron.classification());
    fbrem_e                = electron.fbrem();

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

        // Continue if pdgId is not 11 or status is not 1
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
