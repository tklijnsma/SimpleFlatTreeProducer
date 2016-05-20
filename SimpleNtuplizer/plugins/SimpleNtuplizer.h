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

// For the GenPhoton -- Probably only one of these is necessary
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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

        // void matchPhotonToGenParticle( const reco::PhotonCollection&, const reco::GenParticleCollection& );
        // void matchPhotonToGenParticle( const edm::Handle<reco::PhotonCollection>, const edm::Handle<reco::GenParticleCollection> );
        bool matchPhotonToGenParticle( const reco::Photon& );
        bool matchElectronToGenParticle( const reco::GsfElectron& );

        enum ElectronMatchType{
            UNMATCHED = 0, 
            TRUE_PROMPT_ELECTRON, 
            TRUE_ELECTRON_FROM_TAU,
            TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // Central event counter (specific to this output tree)
        int eventNumber;

        // =====================================
        // Setting tokens
        
        edm::EDGetTokenT<reco::VertexCollection>      vtxToken_;
        edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_; // For AODSIM samples
        edm::EDGetTokenT<reco::PhotonCollection>      photonToken_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
        edm::EDGetTokenT<double>                      rhoToken_; 

        edm::Handle<reco::GenParticleCollection> genParticles;


        // =====================================
        // Event variables

        TTree *eventTree_;

        Int_t nPV_;         // Number of reconsrtucted primary vertices
        Int_t nElectrons_;
        Int_t nElectronsMatched_;
        Int_t nPhotons_;
        Int_t nPhotonsMatched_;


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

        Float_t trkMomentum;
        Float_t trkMomentumError;

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
        // std::vector<Float_t> MaxDRclusterDR_;
        // std::vector<Float_t> MaxDRclusterDPhi_;
        // std::vector<Float_t> MaxDRclusterDEta_;
        // std::vector<Float_t> MaxDRclusterRawEnergy_;

        // Now always contains one entry, but this can be 999.!
        Float_t MaxDRclusterDR;
        Float_t MaxDRclusterDPhi;
        Float_t MaxDRclusterDEta;
        Float_t MaxDRclusterRawEnergy;

        std::vector<Float_t> clusterRawEnergy_;
        std::vector<Float_t> clusterDPhiToSeed_;
        std::vector<Float_t> clusterDEtaToSeed_;


        // =====================================
        // E-p tree variables

        TTree *EpTree_;

        Float_t totEnergyEp_;
        Float_t epEp_;
        Float_t epErrorEp_;
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

        // Matching variables
        Float_t match_dR;
        Float_t match_dE;
        Float_t match_dRdE;

        Float_t gen_pt;
        Float_t gen_phi;
        Float_t gen_eta;
        Float_t gen_M;
        Float_t gen_E;
        Float_t gen_pdgId;
        Float_t gen_status;

    };

