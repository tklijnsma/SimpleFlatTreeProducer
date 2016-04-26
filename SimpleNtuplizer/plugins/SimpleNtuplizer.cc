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
#include "DataFormats/EgammaCandidates/interface/Photon.h"

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

        void setElectronVariables( const reco::GsfElectron&, const edm::Event&, const edm::EventSetup& );
        void setPhotonVariables(   const reco::Photon&,      const edm::Event&, const edm::EventSetup& );

        enum ElectronMatchType{
            UNMATCHED = 0, 
            TRUE_PROMPT_ELECTRON, 
            TRUE_ELECTRON_FROM_TAU,
            TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // =====================================
        // Setting tokens
        
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;

        // Defining electronToken from the electron collection
        //edm::EDGetTokenT<pat::ElectronCollection> electronToken_;   // For miniAOD samples
        edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_; // For AODSIM samples

        edm::EDGetTokenT<reco::PhotonCollection> photonToken_;

        edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
        edm::EDGetTokenT<double> rhoToken_; 


        // =====================================
        // Event variables

        TTree *eventTree_;

        Int_t nPV_;         // Number of reconsrtucted primary vertices
        Int_t nElectrons_;
        Int_t nPhotons_;


        // =====================================
        // Electron tree variables

        TTree *electronTree_;

        // Basic electron variables
        Float_t pt_;
        Float_t rawEnergySC_;
        Float_t etaSC_;
        Float_t phiSC_;
        Float_t etaWidthSC_;
        Float_t phiWidthSC_;
        Float_t r9SS_;
        Float_t seedEnergySS_;
        Float_t eMaxSS_;
        Float_t e2ndSS_;
        Float_t eHorizontalSS_;
        Float_t eVerticalSS_;
        Float_t sigmaIetaIetaSS_;
        Float_t sigmaIetaIphiSS_;
        Float_t sigmaIphiIphiSS_;
        Int_t numberOfClustersSC_;
        Int_t isEB_;
        Float_t preshowerEnergy_;

        // Currently either 0 or 1 entry, depending on whether event is EB or EE
        std::vector<Float_t> iEtaCoordinate_;
        std::vector<Float_t> iPhiCoordinate_;
        std::vector<Float_t> cryEtaCoordinate_;
        std::vector<Float_t> cryPhiCoordinate_;

        std::vector<Float_t> iXCoordinate_;
        std::vector<Float_t> iYCoordinate_;
        std::vector<Float_t> cryXCoordinate_;
        std::vector<Float_t> cryYCoordinate_;


        // Cluster variables

        // These contain either 0 or 1 entries
        std::vector<Float_t> MaxDRclusterDR_;
        std::vector<Float_t> MaxDRclusterDPhi_;
        std::vector<Float_t> MaxDRclusterDEta_;
        std::vector<Float_t> MaxDRclusterRawEnergy_;

        std::vector<Float_t> clusterRawEnergy_;
        std::vector<Float_t> clusterDPhiToSeed_;
        std::vector<Float_t> clusterDEtaToSeed_;


        // =====================================
        // E-p tree variables

        TTree *EpTree_;

        Float_t totEnergyEp_;
        Float_t epEp_;
        Float_t epRelErrorEp_;
        Float_t ecalDrivenEp_;
        Float_t trackerDrivenSeedEp_;
        Int_t classificationEp_;
        Int_t isEBEp_;


        // =====================================
        // Photon tree variables

        TTree *photonTree_;

        Float_t ph_rawEnergy_           ;
        Float_t ph_r9_                  ;
        Float_t ph_etaWidth_            ;
        Float_t ph_phiWidth_            ;
        Int_t   ph_numberOfClusters_    ;
        Float_t ph_hadronicOverEm_      ;
        Float_t ph_rhoValue_            ;
        Float_t ph_delEtaSeed_          ;
        Float_t ph_delPhiSeed_          ;
        Float_t ph_seedEnergy_          ;
        Float_t ph_3x3_5x5_             ;
        Float_t ph_sigmaIetaIeta_       ;
        Float_t ph_sigmaIphiIphi_       ;
        Float_t ph_sigmaIetaIphi_       ;
        Float_t ph_Emax_5x5_            ;
        Float_t ph_e2nd_5x5_            ;
        Float_t ph_eTop_5x5_            ;
        Float_t ph_eBottom_5x5_         ;
        Float_t ph_eLeft_5x5_           ;
        Float_t ph_eRight_5x5_          ;
        Float_t ph_e2x5Max_5x5_         ;
        Float_t ph_e2x5Left_5x5_        ;
        Float_t ph_e2x5Right_5x5_       ;
        Float_t ph_e2x5Top_5x5_         ;
        Float_t ph_e2x5Bottom_5x5_      ;

        Int_t ph_isEB_;
        std::vector<Float_t> ph_5x5_seedEnergy_         ;
        std::vector<Int_t>   ph_iEtaCoordinate_         ;
        std::vector<Int_t>   ph_iPhiCoordinate_         ;
        std::vector<Int_t>   ph_iEtaMod5_               ;
        std::vector<Int_t>   ph_iPhiMod2_               ;
        std::vector<Int_t>   ph_iEtaMod20_              ;
        std::vector<Int_t>   ph_iPhiMod20_              ;
        std::vector<Float_t> ph_preShowerE_rawEnergy_   ;
        std::vector<Float_t> ph_preShowerEp1_rawEnergy_ ;
        std::vector<Float_t> ph_preShowerEp2_rawEnergy_ ;
        std::vector<Int_t>   ph_iXCoordinate_           ;
        std::vector<Int_t>   ph_iYCoordinate_           ;

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

    photonToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),

    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
    {

    std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

    edm::Service<TFileService> fs;


    // =====================================
    // Making event variable tree and setting branches
    // Event variables include the quantities of which there are exactly one per event

    eventTree_ = fs->make<TTree> ("EventTree", "Per event data");

    // Event variables
    eventTree_->Branch( "nPV",        &nPV_,        "nPV/I"   );
    eventTree_->Branch( "nElectrons", &nElectrons_, "nEle/I"  );
    eventTree_->Branch( "nPhotons",   &nPhotons_,   "nPho/I"  );


    // =====================================
    // Making electron tree and setting branches

    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    // Electron variables
    electronTree_->Branch( "pt_",                 &pt_    );
    electronTree_->Branch( "etaSC_",              &etaSC_ );
    electronTree_->Branch( "phiSC_",              &phiSC_ );
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
    electronTree_->Branch( "preshowerEnergy_",    &preshowerEnergy_ );
    electronTree_->Branch( "numberOfClustersSC_", &numberOfClustersSC_ );
    electronTree_->Branch( "isEB_",               &isEB_ );

    // These for EB electrons
    electronTree_->Branch( "iEtaCoordinate_",     &iEtaCoordinate_ );         
    electronTree_->Branch( "iPhiCoordinate_",     &iPhiCoordinate_ );         
    electronTree_->Branch( "cryEtaCoordinate_",   &cryEtaCoordinate_ );
    electronTree_->Branch( "cryPhiCoordinate_",   &cryPhiCoordinate_ );

    // These for EE electrons
    electronTree_->Branch( "iXCoordinate_",       &iXCoordinate_ );         
    electronTree_->Branch( "iYCoordinate_",       &iYCoordinate_ );         
    electronTree_->Branch( "cryXCoordinate_",     &cryXCoordinate_ );
    electronTree_->Branch( "cryYCoordinate_",     &cryYCoordinate_ );

    // Cluster variables
    electronTree_->Branch( "MaxDRclusterDR_",        &MaxDRclusterDR_ );
    electronTree_->Branch( "MaxDRclusterDPhi_",      &MaxDRclusterDPhi_ );
    electronTree_->Branch( "MaxDRclusterDEta_",      &MaxDRclusterDEta_ );
    electronTree_->Branch( "MaxDRclusterRawEnergy_", &MaxDRclusterRawEnergy_ );

    electronTree_->Branch( "clusterRawEnergy_",      &clusterRawEnergy_ );        
    electronTree_->Branch( "clusterDPhiToSeed_",     &clusterDPhiToSeed_ );         
    electronTree_->Branch( "clusterDEtaToSeed_",     &clusterDEtaToSeed_ );


    // =====================================
    // Making E-p tree and setting branches

    EpTree_ = fs->make<TTree> ("EpTree", "E-p data");

    EpTree_->Branch( "totEnergyEp_",         &totEnergyEp_ );
    EpTree_->Branch( "epEp_",                &epEp_ );
    EpTree_->Branch( "epRelErrorEp_",        &epRelErrorEp_ );
    EpTree_->Branch( "ecalDrivenEp_",        &ecalDrivenEp_ );
    EpTree_->Branch( "trackerDrivenSeedEp_", &trackerDrivenSeedEp_ );
    EpTree_->Branch( "classificationEp_",    &classificationEp_ );
    EpTree_->Branch( "isEBEp_",              &isEBEp_ );


    // =====================================
    // Making photon tree and setting branches

    photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");

    photonTree_->Branch( "rawEnergy_",        &ph_rawEnergy_        );
    photonTree_->Branch( "r9_",               &ph_r9_               );
    photonTree_->Branch( "etaWidth_",         &ph_etaWidth_         );
    photonTree_->Branch( "phiWidth_",         &ph_phiWidth_         );
    photonTree_->Branch( "numberOfClusters_", &ph_numberOfClusters_ );                     
    photonTree_->Branch( "hadronicOverEm_",   &ph_hadronicOverEm_   );
    photonTree_->Branch( "rhoValue_",         &ph_rhoValue_         );
    photonTree_->Branch( "delEtaSeed_",       &ph_delEtaSeed_       );
    photonTree_->Branch( "delPhiSeed_",       &ph_delPhiSeed_       );
    photonTree_->Branch( "seedEnergy_",       &ph_seedEnergy_       );
    photonTree_->Branch( "3x3_5x5_",          &ph_3x3_5x5_          );
    photonTree_->Branch( "sigmaIetaIeta_",    &ph_sigmaIetaIeta_    );
    photonTree_->Branch( "sigmaIphiIphi_",    &ph_sigmaIphiIphi_    );
    photonTree_->Branch( "sigmaIetaIphi_",    &ph_sigmaIetaIphi_    );
    photonTree_->Branch( "Emax_5x5_",         &ph_Emax_5x5_         );
    photonTree_->Branch( "e2nd_5x5_",         &ph_e2nd_5x5_         );
    photonTree_->Branch( "eTop_5x5_",         &ph_eTop_5x5_         );
    photonTree_->Branch( "eBottom_5x5_",      &ph_eBottom_5x5_      );
    photonTree_->Branch( "eLeft_5x5_",        &ph_eLeft_5x5_        );
    photonTree_->Branch( "eRight_5x5_",       &ph_eRight_5x5_       );
    photonTree_->Branch( "e2x5Max_5x5_",      &ph_e2x5Max_5x5_      );
    photonTree_->Branch( "e2x5Left_5x5_",     &ph_e2x5Left_5x5_     );
    photonTree_->Branch( "e2x5Right_5x5_",    &ph_e2x5Right_5x5_    );
    photonTree_->Branch( "e2x5Top_5x5_",      &ph_e2x5Top_5x5_      );
    photonTree_->Branch( "e2x5Bottom_5x5_",   &ph_e2x5Bottom_5x5_   );

    // Coordinate variables
    photonTree_->Branch( "isEB_",                   &ph_isEB_                   );
    photonTree_->Branch( "5x5_seedEnergy_",         &ph_5x5_seedEnergy_         );
    photonTree_->Branch( "iEtaCoordinate_",         &ph_iEtaCoordinate_         );
    photonTree_->Branch( "iPhiCoordinate_",         &ph_iPhiCoordinate_         );
    photonTree_->Branch( "iEtaMod5_",               &ph_iEtaMod5_               );
    photonTree_->Branch( "iPhiMod2_",               &ph_iPhiMod2_               );
    photonTree_->Branch( "iEtaMod20_",              &ph_iEtaMod20_              );
    photonTree_->Branch( "iPhiMod20_",              &ph_iPhiMod20_              );
    photonTree_->Branch( "preShowerE_rawEnergy_",   &ph_preShowerE_rawEnergy_   );
    photonTree_->Branch( "preShowerEp1_rawEnergy_", &ph_preShowerEp1_rawEnergy_ );
    photonTree_->Branch( "preShowerEp2_rawEnergy_", &ph_preShowerEp2_rawEnergy_ );
    photonTree_->Branch( "iXCoordinate_",           &ph_iXCoordinate_           );
    photonTree_->Branch( "iYCoordinate_",           &ph_iYCoordinate_           );


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

    // // Get electron collection
    // //edm::Handle<pat::ElectronCollection> electrons;    // For miniAOD
    // edm::Handle<reco::GsfElectronCollection> electrons;  // For AODSIM
    // iEvent.getByToken(electronToken_, electrons);

    // // Loop over electrons
    // nElectrons_ = 0;
    // //for (const pat::Electron &el : *electrons) {
    // for (const reco::GsfElectron &el : *electrons) {

    //     // Increase count of electrons in event
    //     nElectrons_++;

    //     // Fill in the class variables for this particle; also sets E-p variables
    //     setElectronVariables( el, iEvent, iSetup );

    //     // Write class variables to the output tree
    //     electronTree_->Fill();

    //     // Write E-p variables to the E-p tree
    //     EpTree_->Fill();

    //     }

    // Get photon collection
    edm::Handle<reco::PhotonCollection> photons;
    iEvent.getByToken(photonToken_, photons);

     // Loop over photons
    nPhotons_ = 0;
    for (const reco::Photon &photon : *photons) {

        // Increase count of photons in event
        nPhotons_++;

        // Fill in the class variables for this particle; also sets E-p variables
        setPhotonVariables( photon, iEvent, iSetup );

        // Write class variables to the output tree
        photonTree_->Fill();

        // Write E-p variables to the E-p tree
        //EpTree_->Fill();

        }


    // Fill in the event specific variables
    eventTree_->Fill();

    }


//######################################
//# Class methods
//######################################

// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setElectronVariables(
        const reco::GsfElectron& electron,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

    using namespace std;
    using namespace edm;
    using namespace reco;

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

    }


// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setPhotonVariables(
        const reco::Photon& photon,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

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

