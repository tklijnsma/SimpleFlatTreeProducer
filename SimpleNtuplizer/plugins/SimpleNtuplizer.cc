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

//######################################
//# Includes
//######################################

// system include files
#include <memory>
#include <vector>

// To check types - only use if needed
#include <typeinfo>

// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

// From Rafaels standard imports
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"


//######################################
//# Class declaration
//######################################

class SimpleNtuplizer : public edm::EDAnalyzer {
    public:
        explicit SimpleNtuplizer(const edm::ParameterSet&);
        ~SimpleNtuplizer();

        enum ElectronMatchType{
            UNMATCHED = 0, 
            TRUE_PROMPT_ELECTRON, 
            TRUE_ELECTRON_FROM_TAU,
            TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        void ClearBranchVariables();

        // =====================================
        // Setting tokens
        
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;

        // Defining electronToken from the electron collection
        //edm::EDGetTokenT<pat::ElectronCollection> electronToken_;   // For miniAOD samples
        edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_; // For AODSIM samples

        edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
        edm::EDGetTokenT<double> rhoToken_; 


        // =====================================
        // Class variables

        TTree *electronTree_;

        // Vars for PVs
        //Int_t pvNTracks_;

        // Number of reconsrtucted primary vertices
        Int_t nPV_;

        // electron variables
        Int_t nElectrons_;

        // Probably still needed in case of pt cuts
        std::vector<Float_t> pt_;

        // Basic electron variables
        std::vector<Float_t> rawEnergySC_;
        std::vector<Float_t> etaSC_;
        std::vector<Float_t> phiSC_;
        std::vector<Float_t> etaWidthSC_;
        std::vector<Float_t> phiWidthSC_;
        std::vector<Float_t> r9SS_;
        std::vector<Float_t> seedEnergySS_;
        std::vector<Float_t> eMaxSS_;
        std::vector<Float_t> e2ndSS_;
        std::vector<Float_t> eHorizontalSS_;
        std::vector<Float_t> eVerticalSS_;
        std::vector<Float_t> sigmaIetaIetaSS_;
        std::vector<Float_t> sigmaIetaIphiSS_;
        std::vector<Float_t> sigmaIphiIphiSS_;
        std::vector<Int_t> numberOfClustersSC_;

        // Cluster variables
        std::vector<Float_t> MaxDRclusterDR_;
        std::vector<Float_t> MaxDRclusterDPhi_;
        std::vector<Float_t> MaxDRclusterDEta_;
        std::vector<Float_t> MaxDRclusterRawEnergy_;

        std::vector<Float_t> clusterRawEnergy_;
        std::vector<Float_t> clusterDPhiToSeed_;
        std::vector<Float_t> clusterDEtaToSeed_;

        std::vector<Float_t> iEtaCoordinate_;
        std::vector<Float_t> iPhiCoordinate_;
        std::vector<Float_t> cryEtaCoordinate_;
        std::vector<Float_t> cryPhiCoordinate_;

        std::vector<Float_t> preshowerEnergy_;
    };


//######################################
//# Constructors and Destructor
//######################################

// Constructor
SimpleNtuplizer::SimpleNtuplizer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    pileupToken_(consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup"))),

    // For miniAOD:
    //electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    // For AODSIM:
    electronToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),

    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
    {

    std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

    edm::Service<TFileService> fs;

    // Making tree and setting branches
    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    electronTree_->Branch( "pt",   &pt_    );
    electronTree_->Branch( "nPV",  &nPV_        , "nPV/I");
    electronTree_->Branch( "nEle", &nElectrons_ , "nEle/I");

    electronTree_->Branch( "etaSC",               &etaSC_ );
    electronTree_->Branch( "phiSC",               &phiSC_ );
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
    electronTree_->Branch( "numberOfClustersSC_", &numberOfClustersSC_ );

    // Cluster variables
    electronTree_->Branch( "MaxDRclusterDR_",        &MaxDRclusterDR_ );
    electronTree_->Branch( "MaxDRclusterDPhi_",      &MaxDRclusterDPhi_ );
    electronTree_->Branch( "MaxDRclusterDEta_",      &MaxDRclusterDEta_ );
    electronTree_->Branch( "MaxDRclusterRawEnergy_", &MaxDRclusterRawEnergy_ );

    electronTree_->Branch( "clusterRawEnergy_",      &clusterRawEnergy_ );        
    electronTree_->Branch( "clusterDPhiToSeed_",     &clusterDPhiToSeed_ );         
    electronTree_->Branch( "clusterDEtaToSeed_",     &clusterDEtaToSeed_ );         

    electronTree_->Branch( "iEtaCoordinate_",   &iEtaCoordinate_ );         
    electronTree_->Branch( "iPhiCoordinate_",   &iPhiCoordinate_ );         
    electronTree_->Branch( "cryEtaCoordinate_", &cryEtaCoordinate_ );           
    electronTree_->Branch( "cryPhiCoordinate_", &cryPhiCoordinate_ );

    electronTree_->Branch( "preshowerEnergy_",  &preshowerEnergy_ );

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

void SimpleNtuplizer::analyze(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup
    ){
    
    using namespace std;
    using namespace edm;
    using namespace reco;

    // // Pruned particles are the one containing "important" stuff
    // Handle<edm::View<reco::GenParticle> > prunedGenParticles;
    // iEvent.getByToken(prunedGenToken_,prunedGenParticles);

    // // Packed particles are all the status 1, so usable to remake jets
    // // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    // Handle<edm::View<pat::PackedGenParticle> > packedGenParticles;
    // iEvent.getByToken(packedGenToken_,packedGenParticles);


    // =====================================
    // Determine number of primary vertices

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);

    if (vertices->empty()) nPV_ = 0;
    else nPV_ = vertices->size();


    //######################################
    //# Analyze electrons
    //######################################

    // Get electron collection
    //edm::Handle<pat::ElectronCollection> electrons;    // For miniAOD
    edm::Handle<reco::GsfElectronCollection> electrons;  // For AODSIM
    iEvent.getByToken(electronToken_, electrons);

    // Clear variables from previous event
    ClearBranchVariables();

    // Temporary definitions
    double rawEnergy;       // Raw energy of the super cluster per electron
    int numberOfClusters;   // Number of (sub)clusters in the super cluster

    // Loop over electrons
    //for (const pat::Electron &el : *electrons) {
    for (const reco::GsfElectron &el : *electrons) {

        // Increase count of electrons in event
        nElectrons_++;

        // Open variable superCluster
        const reco::SuperClusterRef& superCluster = el.superCluster();

        // Raw energy of the super cluster per electron
        rawEnergy = superCluster->rawEnergy();
        rawEnergySC_     .push_back( rawEnergy );

        pt_              .push_back( el.pt() );
        etaSC_           .push_back( superCluster->eta() );
        phiSC_           .push_back( superCluster->phi() );
        etaWidthSC_      .push_back( superCluster->etaWidth() );
        phiWidthSC_      .push_back( superCluster->phiWidth() );

        r9SS_            .push_back( el.showerShape().r9 );

        seedEnergySS_    .push_back( superCluster->seed()->energy() / rawEnergy );

        eMaxSS_          .push_back( el.showerShape().eMax / rawEnergy );
        e2ndSS_          .push_back( el.showerShape().e2nd / rawEnergy );

        eHorizontalSS_   .push_back(
            el.showerShape().eLeft + el.showerShape().eRight != 0.f  
            ? ( el.showerShape().eLeft - el.showerShape().eRight ) /
                ( el.showerShape().eLeft + el.showerShape().eRight ) : 0.f  );

        eVerticalSS_     .push_back(
            el.showerShape().eTop + el.showerShape().eBottom != 0.f  
            ? ( el.showerShape().eTop - el.showerShape().eBottom ) /
                ( el.showerShape().eTop + el.showerShape().eBottom ) : 0.f  );

        sigmaIetaIetaSS_ .push_back( el.showerShape().sigmaIetaIeta );
        sigmaIetaIphiSS_ .push_back( el.showerShape().sigmaIetaIphi );
        sigmaIphiIphiSS_ .push_back( el.showerShape().sigmaIphiIphi );

        numberOfClusters = std::max( 0, int (superCluster->clusters().size()) );
        numberOfClustersSC_ .push_back( numberOfClusters );


        // =====================================
        // Cluster variables (subs of the superCluster)

        // Default values
        float MaxDRclusterDR        = 999.;
        float MaxDRclusterDPhi      = 999.;
        float MaxDRclusterDEta      = 999.;
        float MaxDRclusterRawEnergy = 0.;

        // Current maximum deltaR; initially 0.0
        float maxDR = 0;

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

            // Find the cluster that has maximum deltaR to the seed
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

        // Write cluster variables to class member vectors
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

        // Currently, this code writes coordinates if the event is EB, and
        // preshower energy if the event is EE, which is in line with the example code,
        // but seems silly.

        // Open up temporary variables
        int iPhi, iEta; float cryPhi, cryEta, dummy;
        EcalClusterLocal ecalLocal;

        if( el.isEB() ){
            ecalLocal.localCoordsEB( *superCluster->seed(), iSetup,
                                      cryEta, cryPhi, iEta, iPhi, dummy, dummy );
            
            iEtaCoordinate_   .push_back( iEta );
            iPhiCoordinate_   .push_back( iPhi );
            cryEtaCoordinate_ .push_back( cryEta );
            cryPhiCoordinate_ .push_back( cryPhi );
            }
        else{
            ecalLocal.localCoordsEE( *superCluster->seed(), iSetup,
                                      cryEta, cryPhi, iEta, iPhi, dummy, dummy );

            preshowerEnergy_  .push_back( superCluster->preshowerEnergy() / superCluster->rawEnergy() );
            }


        //######################################
        //# Analyze EP
        //######################################

        // eval_ep[0] = tot_energy;
        // eval_ep[2] = ep; 
        // eval_ep[3] = trkMomentumRelError;
        // eval_ep[7] = ele.ecalDriven();
        // eval_ep[8] = ele.trackerDrivenSeed();
        // eval_ep[9] = int(ele.classification());//eleClass;
        // eval_ep[10] = iseb;

        // const float tot_energy = superCluster->rawEnergy() + superCluster->preshowerEnergy();
        // totalEnergyMean_ .push_back( (the_sc->rawEnergy()+the_sc->preshowerEnergy()) * mean )

        //const float ep = ele.trackMomentumAtVtx().R()

        }

    // Save this electron's info
    electronTree_->Fill();

    }


//######################################
//# Class methods
//######################################

void SimpleNtuplizer::ClearBranchVariables(){

    nElectrons_ = 0;
    pt_.clear();
    etaSC_.clear();
    phiSC_.clear();
    
    rawEnergySC_.clear();
    etaWidthSC_.clear();
    phiWidthSC_.clear();
    r9SS_.clear();
    seedEnergySS_.clear();
    eMaxSS_.clear();
    e2ndSS_.clear();
    eHorizontalSS_.clear();
    eVerticalSS_.clear();
    sigmaIetaIetaSS_.clear();
    sigmaIetaIphiSS_.clear();
    sigmaIphiIphiSS_.clear();
    numberOfClustersSC_.clear();

    // Cluster variables
    MaxDRclusterDR_.clear();
    MaxDRclusterDPhi_.clear();
    MaxDRclusterDEta_.clear();
    MaxDRclusterRawEnergy_.clear();

    clusterRawEnergy_.clear();
    clusterDPhiToSeed_.clear();
    clusterDEtaToSeed_.clear();

    iEtaCoordinate_.clear();
    iPhiCoordinate_.clear();
    cryEtaCoordinate_.clear();
    cryPhiCoordinate_.clear();

    preshowerEnergy_.clear();

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

