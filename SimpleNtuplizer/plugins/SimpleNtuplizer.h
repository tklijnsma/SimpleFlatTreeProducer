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


        // =====================================
        // Setting tokens
        
        edm::EDGetTokenT<reco::VertexCollection>      vtxToken_;
        edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_; // For AODSIM samples
        edm::EDGetTokenT<reco::PhotonCollection>      photonToken_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
        edm::EDGetTokenT<double>                      rhoToken_; 

        edm::Handle<reco::GenParticleCollection> genParticles_;


        // =====================================
        // Event variables

        TTree *eventTree_;
        TTree *electronTree_;
        TTree *photonTree_;

        // Central event counter (specific to this output tree)
        int NtupID_ = 0;

        Int_t eventNumber_;
        Int_t luminosityBlock_;
        Int_t run_;

        Int_t isMatched_ = 1; // Always 1 but it's there so old regression .root files can be ran

        Int_t nPV_;         // Number of reconstructed primary vertices
        Int_t nElectrons_;
        Int_t nElectronsMatched_;
        Int_t nPhotons_;
        Int_t nPhotonsMatched_;


        //######################################
        //# Variable list for both electrons and photons
        //######################################

        /*
        // Definitive list of variables
        // Copy this for the photon and electron tree

        // -----------------------------
        // Variables used in training

        Float_t pt_;
        Float_t rawEnergy_;
        Float_t eta_;
        Float_t phi_;
        Float_t etaWidth_;
        Float_t phiWidth_;
        Float_t r9_;
        Float_t seedEnergy_;
        Float_t eMax_;
        Float_t e2nd_;
        Float_t eHorizontal_;
        Float_t eVertical_;
        Float_t sigmaIetaIeta_;
        Float_t sigmaIetaIphi_;
        Float_t sigmaIphiIphi_;
        Int_t   numberOfClusters_;
        Int_t   isEB_;
        Float_t preshowerEnergy_;

        // For photon (except eMax and e2nd, which are also in electronTree)
        Float_t hadronicOverEm_;
        Float_t rhoValue_;
        Float_t delEtaSeed_;
        Float_t delPhiSeed_;
        Float_t e5x5_;
        Float_t e3x3_;
        Float_t eMax_;
        Float_t e2nd_;
        Float_t eTop_;
        Float_t eBottom_;
        Float_t eLeft_;
        Float_t eRight_;
        Float_t e2x5Max_;
        Float_t e2x5Left_;
        Float_t e2x5Right_;
        Float_t e2x5Top_;
        Float_t e2x5Bottom_;

        // -----------------------------
        // Coordinate variables

        // EB
        std::vector<Int_t>   iEtaCoordinate_;
        std::vector<Int_t>   iPhiCoordinate_;
        std::vector<Int_t>   cryEtaCoordinate_;
        std::vector<Int_t>   cryPhiCoordinate_;
        // EE
        std::vector<Int_t>   iXCoordinate_;
        std::vector<Int_t>   iYCoordinate_;
        std::vector<Int_t>   cryXCoordinate_;
        std::vector<Int_t>   cryYCoordinate_;

        // Additional coordinate variables for photon
        std::vector<Int_t>   iEtaMod5_;
        std::vector<Int_t>   iPhiMod2_;
        std::vector<Int_t>   iEtaMod20_;
        std::vector<Int_t>   iPhiMod20_;
        std::vector<Float_t> preshowerEnergyPlane1_;
        std::vector<Float_t> preshowerEnergyPlane2_;

        // -----------------------------
        // Cluster variables

        Float_t MaxDRclusterDR_;
        Float_t MaxDRclusterDPhi_;
        Float_t MaxDRclusterDEta_;
        Float_t MaxDRclusterRawEnergy_;

        std::vector<Float_t> clusterRawEnergy_;
        std::vector<Float_t> clusterDPhiToSeed_;
        std::vector<Float_t> clusterDEtaToSeed_;

        // Only for electron
        Int_t   IsEcalEnergyCorrected_;
        Float_t CorrectedEcalEnergy_;
        Float_t CorrectedEcalEnergyError_;

        // -----------------------------
        // Ep variables (only for electron)

        Float_t trkMomentum_;
        Float_t trkMomentumError_;
        Float_t trkMomentumRelError_;
        Float_t ecalDriven_;
        Float_t trackerDriven_;
        Float_t classification_;

        Float_t genMatchdR_;
        Float_t genMatchdE_;
        Float_t genMatchdRdE_;
        Float_t genPt_;
        Float_t genPhi_;
        Float_t genEta_;
        Float_t genMass_;
        Float_t genEnergy_;
        Int_t   genPdgId_;
        Int_t   genStatus_;
        */


        // =====================================
        // Electron tree variables

        // -----------------------------
        // Variables used in training

        Float_t pt_e;
        Float_t rawEnergy_e;
        Float_t eta_e;
        Float_t phi_e;
        Float_t etaWidth_e;
        Float_t phiWidth_e;
        Float_t r9_e;
        Float_t seedEnergy_e;
        Float_t eHorizontal_e;
        Float_t eVertical_e;
        Float_t sigmaIetaIeta_e;
        Float_t sigmaIetaIphi_e;
        Float_t sigmaIphiIphi_e;
        Int_t   numberOfClusters_e;
        Int_t   isEB_e;
        Float_t preshowerEnergy_e;

        // For photon (except eMax and e2nd, which are also in electronTree)
        Float_t hadronicOverEm_e;
        Float_t rhoValue_e;
        Float_t delEtaSeed_e;
        Float_t delPhiSeed_e;
        Float_t e5x5_e;
        Float_t e3x3_e;
        Float_t eMax_e;
        Float_t e2nd_e;
        Float_t eTop_e;
        Float_t eBottom_e;
        Float_t eLeft_e;
        Float_t eRight_e;
        Float_t e2x5Max_e;
        Float_t e2x5Left_e;
        Float_t e2x5Right_e;
        Float_t e2x5Top_e;
        Float_t e2x5Bottom_e;

        // -----------------------------
        // Coordinate variables

        // EB
        std::vector<Int_t>   iEtaCoordinate_e;
        std::vector<Int_t>   iPhiCoordinate_e;
        std::vector<Int_t>   cryEtaCoordinate_e;
        std::vector<Int_t>   cryPhiCoordinate_e;
        // EE
        std::vector<Int_t>   iXCoordinate_e;
        std::vector<Int_t>   iYCoordinate_e;
        std::vector<Int_t>   cryXCoordinate_e;
        std::vector<Int_t>   cryYCoordinate_e;

        // Additional coordinate variables for photon
        std::vector<Int_t>   iEtaMod5_e;
        std::vector<Int_t>   iPhiMod2_e;
        std::vector<Int_t>   iEtaMod20_e;
        std::vector<Int_t>   iPhiMod20_e;
        std::vector<Float_t> preshowerEnergyPlane1_e;
        std::vector<Float_t> preshowerEnergyPlane2_e;

        // -----------------------------
        // Cluster variables

        Float_t MaxDRclusterDR_e;
        Float_t MaxDRclusterDPhi_e;
        Float_t MaxDRclusterDEta_e;
        Float_t MaxDRclusterRawEnergy_e;

        std::vector<Float_t> clusterRawEnergy_e;
        std::vector<Float_t> clusterDPhiToSeed_e;
        std::vector<Float_t> clusterDEtaToSeed_e;

        // Only for electron
        Int_t   IsEcalEnergyCorrected_e;
        Float_t CorrectedEcalEnergy_e;
        Float_t CorrectedEcalEnergyError_e;

        // -----------------------------
        // Ep variables (only for electron)

        Float_t trkMomentum_e;
        Float_t trkMomentumError_e;
        Float_t trkMomentumRelError_e;
        Float_t ecalDriven_e;
        Float_t trackerDriven_e;
        Float_t classification_e;

        Float_t genMatchdR_e;
        Float_t genMatchdE_e;
        Float_t genMatchdRdE_e;
        Float_t genPt_e;
        Float_t genPhi_e;
        Float_t genEta_e;
        Float_t genMass_e;
        Float_t genEnergy_e;
        Int_t   genPdgId_e;
        Int_t   genStatus_e;


        // =====================================
        // Photon tree variables

        // Last minute addition (add same type of objects to electrons tree too)
        Float_t scEcalEnergy_p      ;
        Float_t phoEcalEnergy_p     ;
        Float_t regression1Energy_p ;
        Float_t regression2Energy_p ;
        Float_t scEcalEnergyError_p      ;
        Float_t phoEcalEnergyError_p     ;
        Float_t regression1EnergyError_p ;
        Float_t regression2EnergyError_p ;

        // -----------------------------
        // Variables used in training

        Float_t pt_p;
        Float_t rawEnergy_p;
        Float_t eta_p;
        Float_t phi_p;
        Float_t etaWidth_p;
        Float_t phiWidth_p;
        Float_t r9_p;
        Float_t seedEnergy_p;
        Float_t eHorizontal_p;
        Float_t eVertical_p;
        Float_t sigmaIetaIeta_p;
        Float_t sigmaIetaIphi_p;
        Float_t sigmaIphiIphi_p;
        Int_t   numberOfClusters_p;
        Int_t   isEB_p;
        Float_t preshowerEnergy_p;

        // For photon (except eMax and e2nd, which are also in electronTree)
        Float_t hadronicOverEm_p;
        Float_t rhoValue_p;
        Float_t delEtaSeed_p;
        Float_t delPhiSeed_p;
        Float_t e5x5_p;
        Float_t e3x3_p;
        Float_t eMax_p;
        Float_t e2nd_p;
        Float_t eTop_p;
        Float_t eBottom_p;
        Float_t eLeft_p;
        Float_t eRight_p;
        Float_t e2x5Max_p;
        Float_t e2x5Left_p;
        Float_t e2x5Right_p;
        Float_t e2x5Top_p;
        Float_t e2x5Bottom_p;

        // -----------------------------
        // Coordinate variables

        // EB
        std::vector<Int_t>   iEtaCoordinate_p;
        std::vector<Int_t>   iPhiCoordinate_p;
        std::vector<Int_t>   cryEtaCoordinate_p;
        std::vector<Int_t>   cryPhiCoordinate_p;
        // EE
        std::vector<Int_t>   iXCoordinate_p;
        std::vector<Int_t>   iYCoordinate_p;
        std::vector<Int_t>   cryXCoordinate_p;
        std::vector<Int_t>   cryYCoordinate_p;

        // Additional coordinate variables for photon
        std::vector<Int_t>   iEtaMod5_p;
        std::vector<Int_t>   iPhiMod2_p;
        std::vector<Int_t>   iEtaMod20_p;
        std::vector<Int_t>   iPhiMod20_p;
        std::vector<Float_t> preshowerEnergyPlane1_p;
        std::vector<Float_t> preshowerEnergyPlane2_p;

        // -----------------------------
        // Cluster variables

        Float_t MaxDRclusterDR_p;
        Float_t MaxDRclusterDPhi_p;
        Float_t MaxDRclusterDEta_p;
        Float_t MaxDRclusterRawEnergy_p;

        std::vector<Float_t> clusterRawEnergy_p;
        std::vector<Float_t> clusterDPhiToSeed_p;
        std::vector<Float_t> clusterDEtaToSeed_p;

        // Only for electron
        Int_t   IsEcalEnergyCorrected_p;
        Float_t CorrectedEcalEnergy_p;
        Float_t CorrectedEcalEnergyError_p;

        // -----------------------------
        // Ep variables (only for electron)

        Float_t trkMomentum_p;
        Float_t trkMomentumError_p;
        Float_t trkMomentumRelError_p;
        Float_t ecalDriven_p;
        Float_t trackerDriven_p;
        Float_t classification_p;

        Float_t genMatchdR_p;
        Float_t genMatchdE_p;
        Float_t genMatchdRdE_p;
        Float_t genPt_p;
        Float_t genPhi_p;
        Float_t genEta_p;
        Float_t genMass_p;
        Float_t genEnergy_p;
        Int_t   genPdgId_p;
        Int_t   genStatus_p;

    };

