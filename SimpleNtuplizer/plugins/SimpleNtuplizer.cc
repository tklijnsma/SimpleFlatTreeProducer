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

#include <typeinfo>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
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

// Added by Thomas
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"


//######################################
//# Class declaration
//######################################

class SimpleNtuplizer : public edm::EDAnalyzer {
    public:
        explicit SimpleNtuplizer(const edm::ParameterSet&);
        ~SimpleNtuplizer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        enum ElectronMatchType{
            UNMATCHED = 0, 
            TRUE_PROMPT_ELECTRON, 
            TRUE_ELECTRON_FROM_TAU,
            TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        // Define methods
        bool checkAncestor(const reco::Candidate *gen, int ancestorPid);
        int matchToTruth(const pat::Electron &el, const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
        void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
        int matchToTruthAlternative(const pat::Electron &el);

        void printAllZeroMothers(const reco::Candidate *particle);

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
        edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
        edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
        edm::EDGetTokenT<double> rhoToken_; 

        TTree *electronTree_;

        // Vars for PVs
        Int_t pvNTracks_;

        // Vars for pile-up
        Int_t nPUTrue_;    // true pile-up
        Int_t nPU_;        // generated pile-up
        Int_t nPV_;        // number of reconsrtucted primary vertices
        Float_t rho_;      // the rho variable

        // From example code, nVtx needed:
        // evt.getByToken(vtxToken_, vtxH_);
        // nVtx_ = vtxH_->size();


        // electron variables
        Int_t nElectrons_;


        // Old variables; still needed to compile functions

        std::vector<Float_t> pt_;
        
        std::vector<Float_t> dEtaIn_;
        std::vector<Float_t> dPhiIn_;
        std::vector<Float_t> hOverE_;
        std::vector<Float_t> sigmaIetaIeta_;
        std::vector<Float_t> full5x5_sigmaIetaIeta_;
        std::vector<Float_t> isoChargedHadrons_;
        std::vector<Float_t> isoNeutralHadrons_;
        std::vector<Float_t> isoPhotons_;
        std::vector<Float_t> isoChargedFromPU_;
        std::vector<Float_t> relIsoWithDBeta_;
        std::vector<Float_t> ooEmooP_;
        std::vector<Float_t> d0_;
        std::vector<Float_t> dz_;
        std::vector<Int_t>   expectedMissingInnerHits_;
        std::vector<Int_t>   passConversionVeto_;     
        std::vector<Int_t>   isTrueElectron_;
        std::vector<Int_t>   isTrueElectronAlternative_; 


        //######################################
        // Variables from Thomas
        //######################################

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


        // For the cluster variables

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

    };


//######################################
//# Constructors and Destructor
//######################################

// Constructor
SimpleNtuplizer::SimpleNtuplizer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    pileupToken_(consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
    {

    std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

    edm::Service<TFileService> fs;

    // Making tree and setting branches
    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    // electronTree_->Branch( "pvNTracks",  &pvNTracks_ , "pvNTracks/I");
    electronTree_->Branch( "nPV",  &nPV_     , "nPV/I");
    // electronTree_->Branch( "nPU",  &nPU_     , "nPU/I");
    // electronTree_->Branch( "nPUTrue",  &nPUTrue_ , "nPUTrue/I");
    // electronTree_->Branch( "rho",  &rho_ , "rho/F");
    electronTree_->Branch( "nEle",  &nElectrons_ , "nEle/I");
    // electronTree_->Branch( "pt",  &pt_    );
    // electronTree_->Branch( "dEtaIn",  &dEtaIn_);
    // electronTree_->Branch( "dPhiIn",  &dPhiIn_);
    // electronTree_->Branch( "hOverE",  &hOverE_);
    // electronTree_->Branch( "sigmaIetaIeta",         &sigmaIetaIeta_);
    // electronTree_->Branch( "full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
    // electronTree_->Branch( "isoChargedHadrons", &isoChargedHadrons_);
    // electronTree_->Branch( "isoNeutralHadrons", &isoNeutralHadrons_);
    // electronTree_->Branch( "isoPhotons", &isoPhotons_);
    // electronTree_->Branch( "isoChargedFromPU", &isoChargedFromPU_);
    // electronTree_->Branch( "relIsoWithDBeta", &relIsoWithDBeta_);
    // electronTree_->Branch( "ooEmooP", &ooEmooP_);
    // electronTree_->Branch( "d0", &d0_);
    // electronTree_->Branch( "dz", &dz_);
    // electronTree_->Branch( "expectedMissingInnerHits", &expectedMissingInnerHits_);
    // electronTree_->Branch( "passConversionVeto", &passConversionVeto_);
    // electronTree_->Branch( "isTrueElectron", &isTrueElectron_);
    // electronTree_->Branch( "isTrueElectronAlternative", &isTrueElectronAlternative_);

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
// analyze - is executed every event

void SimpleNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace std;
    using namespace edm;
    using namespace reco;

    // Pruned particles are the one containing "important" stuff
    Handle<edm::View<reco::GenParticle> > prunedGenParticles;
    iEvent.getByToken(prunedGenToken_,prunedGenParticles);

    // Packed particles are all the status 1, so usable to remake jets
    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    Handle<edm::View<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken_,packedGenParticles);


    // =====================================
    // Get Pileup

    // THIS IS CURRENTLY BROKEN - Should check if this is still needed
    // Get Pileup info
    // Handle<edm::View<PileupSummaryInfo> > pileupHandle;
    // cout << ">>>> Trying getByToken call (pileup)" << endl;
    // iEvent.getByToken(pileupToken_, pileupHandle);
    // for( auto & puInfoElement : *pileupHandle){
    //   if( puInfoElement.getBunchCrossing() == 0 ){
    //     nPU_    = puInfoElement.getPU_NumInteractions();
    //     nPUTrue_= puInfoElement.getTrueNumInteractions();
    //   }
    // }
    // vector<PileupSummaryInfo>::const_iterator PVI;
    // for (PVI = pileupHandle->begin(); PVI != pileupHandle->end(); ++PVI) {
    //   if (PVI->getBunchCrossing() == 0) {
    //     nPU_     = PVI->getPU_NumInteractions();
    //     nPUTrue_ = PVI->getTrueNumInteractions();
    //   }
    // }


    // =====================================
    // Get Primary Vertex

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    
    if (vertices->empty()) return; // skip the event if no PV found
    //const reco::Vertex &pv = vertices->front();
    nPV_ = vertices->size();


    // =====================================
    // Loop over vertices - probably only nPV needed

    // VertexCollection::const_iterator firstGoodVertex = vertices->end();
    // int firstGoodVertexIdx = 0;
    // for (VertexCollection::const_iterator vtx = vertices->begin(); 
    //     vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    //     // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    //     // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    //     if (
    //         // !vtx->isFake() &&
    //         !(vtx->chi2()==0 && vtx->ndof()==0) 
    //         &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
    //         && fabs(vtx->position().Z())<=24.0  ){
    //             firstGoodVertex = vtx;
    //             break;
    //             }
    //     }

    // if ( firstGoodVertex==vertices->end() )
    // return; // skip event if there are no good PVs

    // // Seems always zero. Not stored in miniAOD...?
    // pvNTracks_ = firstGoodVertex->nTracks();


    // =====================================
    // Get rho - probably no longer needed

    // edm::Handle< double > rhoH;
    // iEvent.getByToken(rhoToken_,rhoH);
    // rho_ = *rhoH;


    // =====================================
    // Analyze actual electrons

    // Get electron collection
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);

    // Why are we clearing this?
    nElectrons_ = 0;
    pt_.clear();
    etaSC_.clear();
    phiSC_.clear();
    dEtaIn_.clear();
    dPhiIn_.clear();
    hOverE_.clear();
    sigmaIetaIeta_.clear();
    full5x5_sigmaIetaIeta_.clear();
    isoChargedHadrons_.clear();
    isoNeutralHadrons_.clear();
    isoPhotons_.clear();
    isoChargedFromPU_.clear();
    relIsoWithDBeta_.clear();
    ooEmooP_.clear();
    d0_.clear();
    dz_.clear();
    expectedMissingInnerHits_.clear();
    passConversionVeto_.clear();     
    isTrueElectron_.clear();
    isTrueElectronAlternative_.clear(); 

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


    double rawEnergy;       // Raw energy of the super cluster per electron
    int numberOfClusters;   // Number of (sub)clusters in the super cluster

    // Loop over electrons
    for (const pat::Electron &el : *electrons) {

        // Simple pt cut
        if( el.pt() < 10 ) continue;

        // printf( "\n-----------------------------------\n" );
        // printf( "Analyzing electron %d\n", nElectrons_ );

        // Increase count of electrons in event
        nElectrons_++;

        // Open variable superCluster
        const reco::SuperClusterRef& superCluster = el.superCluster();

        // Raw energy of the super cluster per electron
        rawEnergy = superCluster->rawEnergy();
        rawEnergySC_     .push_back( rawEnergy );

        //pt_.push_back( el.pt() );
        etaSC_           .push_back( superCluster->eta() );
        phiSC_           .push_back( superCluster->phi() );
        etaWidthSC_      .push_back( superCluster->etaWidth() );
        phiWidthSC_      .push_back( superCluster->phiWidth() );

        r9SS_            .push_back( el.showerShape().r9 );

        //theseed->energy()/raw_energy
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

        // std::cout << "Found etaWidthSC_ = " << etaWidthSC_[ etaWidthSC_.size()-1 ] << std::endl;
        // std::cout << "Found phiWidthSC_ = " << phiWidthSC_[ phiWidthSC_.size()-1 ] << std::endl;

        // std::cout << "Found sigmaIetaIetaSS_ = "
        //           << sigmaIetaIetaSS_[ sigmaIetaIetaSS_.size()-1 ]
        //           << std::endl;

        // std::cout << "Found eHorizontalSS_ = "
        //           << eHorizontalSS_[ eHorizontalSS_.size()-1 ]
        //           << std::endl;


        // =====================================
        // Cluster variables (subs of the superCluster)

        // --> HIER VERDER
        // die max 3 is wrs beter variabel te maken, en dan gewoon push_back gebruiken
        // naar die vectors i.p.v. direct setten (wat nu gebeurt)

        // Open vectors that contain a quantity per cluster (waarom max 3?)
        std::vector<float> clusterRawEnergy;
        //clusterRawEnergy.resize( std::max(3, numberOfClusters), 0);
        std::vector<float> clusterDEtaToSeed;
        //clusterDEtaToSeed.resize(std::max(3, numberOfClusters), 0);
        std::vector<float> clusterDPhiToSeed;
        //clusterDPhiToSeed.resize(std::max(3, numberOfClusters), 0);

        // Default values
        float MaxDRclusterDR        = 999.;
        float MaxDRclusterDPhi      = 999.;
        float MaxDRclusterDEta      = 999.;
        float MaxDRclusterRawEnergy = 0.;


        // ---------------------------
        // Looping over clusters

        size_t i_cluster = 0;
        
        // Current maximum deltaR; initially 0.0
        float maxDR = 0;
        
        // Define cluster
        edm::Ptr<reco::CaloCluster> cluster;
        
        // Loop over all clusters that aren't the seed  
        for( auto pcluster = superCluster->clustersBegin(); pcluster != superCluster->clustersEnd(); ++pcluster ) {
            
            // Dereference to get the actual object
            cluster = *pcluster;

            // Continue if this is the seed
            if( cluster == superCluster->seed() ) continue;

            // cout << "  Currently in cluster " << i_cluster << endl;
            // cout << "    clusterRawEnergy_ = " << cluster->energy() / rawEnergy << endl;

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
            // Misschien een break als i_cluster > cut-1 ? (vorige cut was 3)

            //if i_cluster == 3
            }

        // Write cluster variables to class member vectors
        if( i_cluster > 0 ) {

            // Still use vectors, since these vars may be empty
            MaxDRclusterDR_        .push_back( MaxDRclusterDR );    
            MaxDRclusterDPhi_      .push_back( MaxDRclusterDPhi );    
            MaxDRclusterDEta_      .push_back( MaxDRclusterDEta );    
            MaxDRclusterRawEnergy_ .push_back( MaxDRclusterRawEnergy );

            // These should be correctly set:
            // clusterRawEnergy_
            // clusterDPhiToSeed_
            // clusterDEtaToSeed_

            }

        // Coordinate variables
        // Does different things for when the electron is in the barrel or endcap
        // Not exactly sure what is what


        // THIS IS BROKEN

        // Open up temporary variables
        int iPhi, iEta; float cryPhi, cryEta, dummy;
        EcalClusterLocal _ecalLocal;


        if( el.isEB() ){

            _ecalLocal.localCoordsEB( *superCluster->seed(), iSetup,
                                      cryEta, cryPhi, iEta, iPhi, dummy, dummy );
            
            iEtaCoordinate_   .push_back( iEta );
            iPhiCoordinate_   .push_back( iPhi );
            cryEtaCoordinate_ .push_back( cryEta );
            cryPhiCoordinate_ .push_back( cryPhi );

            }
        else {

            // Does this not include also barrel-barrel events?
            // Seems strange to do 'localCoordsEE' when the event may be BB

            _ecalLocal.localCoordsEE( *superCluster->seed(), iSetup,
                                      cryEta, cryPhi, iEta, iPhi, dummy, dummy );

            // This part does not really make sense. See example code below.

            }


        // =======================================================================================
        // From example code

        // // calculate coordinate variables
        // const bool iseb = ele.isEB();  
        // float dummy;
        // int iPhi;
        // int iEta;
        // float cryPhi;
        // float cryEta;
        // EcalClusterLocal _ecalLocal;
        // if (ele.isEB()) 
        // _ecalLocal.localCoordsEB(*theseed, *iSetup_, cryEta, cryPhi, iEta, iPhi, dummy, dummy);
        // else 
        // _ecalLocal.localCoordsEE(*theseed, *iSetup_, cryEta, cryPhi, iEta, iPhi, dummy, dummy);

        // if (iseb) {
        // eval[29] = cryEta;
        // eval[30] = cryPhi;
        // eval[31] = iEta;
        // eval[32] = iPhi;
        // } else {
        // eval[29] = the_sc->preshowerEnergy()/the_sc->rawEnergy();
        // }

        }

    // Save this electron's info
    electronTree_->Fill();
    }


// =====================================
// method fills 'descriptions' with the allowed parameters for the module

void SimpleNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
    }

bool SimpleNtuplizer::checkAncestor(const reco::Candidate *gen, int ancestorPid){

    // General sanity check
    if( gen == 0 ){
        printf("SimpleNtuplizer::checkAncestor: ERROR null particle is passed in, ignore it.\n");
        return false;
        }

    // If this is true, we found our target ancestor
    if( abs( gen->pdgId() ) == ancestorPid )
        return true;

    // Go deeper and check all mothers
    for(size_t i=0;i< gen->numberOfMothers();i++) {
        if ( checkAncestor( gen->mother(i), ancestorPid) )
            return true;
        }

    // Return false if not found
    return false;
    }


// =====================================
// Matches to truth level
// Explicit loop and geometric matching method (advised by Josh Bendavid)
// T: does not seem to work currently - alternative does work

int SimpleNtuplizer::matchToTruth(
    const pat::Electron &el, 
    const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles
    ){

    // Find the closest status 1 gen electron to the reco electron
    double dR = 999;
    const reco::Candidate *closestElectron = 0;
    
    for(size_t i=0; i<prunedGenParticles->size();i++){

        // Get the particle
        const reco::Candidate *particle = &(*prunedGenParticles)[i];

        // Drop everything that is not electron or not status 1
        if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
            continue;

        // Calculate delta R
        double dRtmp = ROOT::Math::VectorUtil::DeltaR( el.p4(), particle->p4() );

        if( dRtmp < dR ){
            dR = dRtmp;
            closestElectron = particle;
            }
        }

    // See if the closest electron (if it exists) is close enough.
    // If not, no match found.
    if( !(closestElectron != 0 && dR < 0.1) ) {
        return UNMATCHED;
    }


    int ancestorPID = -999; 
    int ancestorStatus = -999;
    findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

    if( ancestorPID == -999 && ancestorStatus == -999 ){
        // No non-electron parent??? This should never happen.
        // Complain.
        printf("SimpleNtuplizer: ERROR! Electron does not apper to have a non-electron parent\n");
        return UNMATCHED;
        }

    if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
        return TRUE_NON_PROMPT_ELECTRON;

    if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
        return TRUE_ELECTRON_FROM_TAU;

    // What remains is true prompt electrons
        return TRUE_PROMPT_ELECTRON;

    }


// =====================================
// Alternative truth matching

int SimpleNtuplizer::matchToTruthAlternative( const pat::Electron &el ){

    int result = UNMATCHED;

    const reco::GenParticle * gen = el.genParticle();
    
    if( gen != 0 ){
        
        int pid = gen->pdgId();
        int status = gen->status();
        bool isFromZ   = checkAncestor(gen, 23);
        bool isFromW   = checkAncestor(gen, 24);
        bool isFromTau = checkAncestor(gen, 15);
        
        // Check if it is a true prompt electron
        if( abs( pid ) == 11 // this is electron
            && (status == 1 || status == 22 || status == 23 ) // NOTE: Pythia8 status here 22/23 (for Pythia6 would be 3)
            && (isFromZ || isFromW ) && !isFromTau // comes from Z or W+-, but not from tau
            ){
                result = TRUE_PROMPT_ELECTRON;
                }

        else if (
            abs( pid ) == 11 
            && (status == 1 || status == 22 || status == 23 ) 
            && (isFromTau ) 
            ){
                // This is a true electron, but it comes from tau
                result = TRUE_ELECTRON_FROM_TAU;
                }

        else if ( abs( pid ) == 11 ) {
            // This is a true electron, but it comes from something else
            const reco::Candidate *mom = el.mother(0);
            int momPid = -999;

            if ( mom != 0 ) momPid = mom->pdgId();
            
            printf("pid= %d  status= %d isFromZ= %d isFromW= %d  isFromTau= %d  momPid= %d\n",
                pid,  status, isFromZ, isFromW, isFromTau, momPid);
            
            result = TRUE_NON_PROMPT_ELECTRON;
            }

        else {
            printf("The reco electron has a truth match with pid= %d\n", pid);
            }

        }

    return result;
    }


// =====================================
// Find first non-electron mother

void SimpleNtuplizer::findFirstNonElectronMother(
    const reco::Candidate *particle,
    int &ancestorPID,
    int &ancestorStatus
    ){

    if( particle == 0 ){
        printf("SimpleNtuplizer: ERROR! null candidate pointer, this should never happen\n");
        return;
        }

    // Is this the first non-electron parent? If yes, return, otherwise
    // go deeper into recursion
    if( abs(particle->pdgId()) == 11 ){
        findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
        }
    else {
        ancestorPID = particle->pdgId();
        ancestorStatus = particle->status();
        }

    return;
    }


// =====================================
// Print all zero mothers

void SimpleNtuplizer::printAllZeroMothers(const reco::Candidate *particle){

    if( particle == 0 ){
        printf("SimpleNtuplizer::printAllZeroMothers: reached the top of the decay tree\n");
        return;
        }

    printf("SimpleNtuplizer::printAllZeroMothers: ancestor ID= %d, status= %d\n",
        particle->pdgId(), particle->status() );

    printAllZeroMothers( particle->mother(0) );

    return;
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

