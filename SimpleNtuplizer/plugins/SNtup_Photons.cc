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
    // Fill other variables

    // Set convenience variables
    const reco::SuperClusterRef& superCluster = photon.superCluster();
    const edm::Ptr<reco::CaloCluster>& seedCluster = superCluster->seed();
    const int numberOfClusters =  superCluster->clusters().size();

    // T: Not sure if this should be done!!
    // const bool missing_clusters = !superCluster->clusters()[numberOfClusters-1].isAvailable();
    // if( missing_clusters ) return ; // do not apply corrections in case of missing info (slimmed MiniAOD electrons)

    const double rawEnergy = superCluster->rawEnergy(); 
    const auto& showerShape = photon.showerShapeVariables();

    // Get rho
    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken_,rhoH);
    Float_t rho = *rhoH;

    // Fill in variables that will be written to the tree
    ph_rawEnergy_         = rawEnergy;
    ph_r9_                = photon.r9();
    ph_etaWidth_          = superCluster->etaWidth();
    ph_phiWidth_          = superCluster->phiWidth(); 
    ph_numberOfClusters_  = std::max(0,numberOfClusters - 1);
    ph_hadronicOverEm_    = photon.hadronicOverEm();
    ph_rhoValue_          = rho;
    ph_delEtaSeed_        = seedCluster->eta()-superCluster->position().Eta();
    ph_delPhiSeed_        = reco::deltaPhi(seedCluster->phi(),superCluster->position().Phi());
    ph_seedEnergy_        = seedCluster->energy()/rawEnergy;
    ph_3x3_5x5_           = showerShape.e3x3 / showerShape.e5x5;
    ph_sigmaIetaIeta_     = showerShape.sigmaIetaIeta;  
    ph_sigmaIphiIphi_     = showerShape.sigmaIphiIphi;
    ph_sigmaIetaIphi_     = showerShape.sigmaIetaIphi / (showerShape.sigmaIphiIphi*showerShape.sigmaIetaIeta);
    ph_Emax_5x5_          = showerShape.maxEnergyXtal / showerShape.e5x5;
    ph_e2nd_5x5_          = showerShape.e2nd       / showerShape.e5x5;
    ph_eTop_5x5_          = showerShape.eTop       / showerShape.e5x5;
    ph_eBottom_5x5_       = showerShape.eBottom    / showerShape.e5x5;
    ph_eLeft_5x5_         = showerShape.eLeft      / showerShape.e5x5;
    ph_eRight_5x5_        = showerShape.eRight     / showerShape.e5x5;  
    ph_e2x5Max_5x5_       = showerShape.e2x5Max    / showerShape.e5x5;
    ph_e2x5Left_5x5_      = showerShape.e2x5Left   / showerShape.e5x5;
    ph_e2x5Right_5x5_     = showerShape.e2x5Right  / showerShape.e5x5;
    ph_e2x5Top_5x5_       = showerShape.e2x5Top    / showerShape.e5x5;
    ph_e2x5Bottom_5x5_    = showerShape.e2x5Bottom / showerShape.e5x5;

    // Clear the std::vectors from last photon
    ph_5x5_seedEnergy_         .clear();
    ph_iEtaCoordinate_         .clear();
    ph_iPhiCoordinate_         .clear();
    ph_iEtaMod5_               .clear();
    ph_iPhiMod2_               .clear();
    ph_iEtaMod20_              .clear();
    ph_iPhiMod20_              .clear();
    ph_preShowerE_rawEnergy_   .clear();
    ph_preShowerEp1_rawEnergy_ .clear();
    ph_preShowerEp2_rawEnergy_ .clear();
    ph_iXCoordinate_           .clear();
    ph_iYCoordinate_           .clear();

    // Coordinate variables
    ph_isEB_ = photon.isEB();
    if ( photon.isEB() ) {
        
        ph_5x5_seedEnergy_ .push_back( photon.e5x5() / seedCluster->energy() ); // Why does photon.e5x5() work here? Up it's ss.e5x5

        EBDetId ebseedid( seedCluster->seed());
        
        int ieta = ebseedid.ieta();
        int iphi = ebseedid.iphi();
        int signieta = ieta > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training

        ph_iEtaCoordinate_ .push_back( ieta );
        ph_iPhiCoordinate_ .push_back( iphi );
        ph_iEtaMod5_       .push_back( (ieta-signieta)%5 );
        ph_iPhiMod2_       .push_back( (iphi-1)%2 );
        ph_iEtaMod20_      .push_back( abs(ieta)<=25 ? ieta-signieta : (ieta - 26*signieta) % 20  );
        ph_iPhiMod20_      .push_back( (iphi-1)%20 );
        }
    
    else {
        EEDetId eeseedid( seedCluster->seed());
        ph_preShowerE_rawEnergy_   .push_back( superCluster->preshowerEnergy() / rawEnergy );
        ph_preShowerEp1_rawEnergy_ .push_back( superCluster->preshowerEnergyPlane1() / rawEnergy );
        ph_preShowerEp2_rawEnergy_ .push_back( superCluster->preshowerEnergyPlane2() / rawEnergy );
        ph_iXCoordinate_           .push_back( eeseedid.ix() );
        ph_iYCoordinate_           .push_back( eeseedid.iy() );
        }


    // Write class variables to the output tree
    photonTree_->Fill();

    // TODO: What to do with the Ep variables?

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

    for (const reco::GenParticle &genParticle : *genParticles) {

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

    match_dR     = minDr;
    match_dE     = minDe;
    match_dRdE   = minDeDr;
    gen_pt       = matched_genParticle->pt();
    gen_phi      = matched_genParticle->phi();
    gen_eta      = matched_genParticle->eta();
    gen_M        = matched_genParticle->mass();
    gen_E        = matched_genParticle->energy();
    gen_pdgId    = matched_genParticle->pdgId();
    gen_status   = matched_genParticle->status();

    // Return successful match value (should be true)
    return successful_match;

    }

