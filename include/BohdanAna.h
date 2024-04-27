#ifndef BohdanAna_h
#define BohdanAna_h 1

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "EventDisplayer.h"
#include "TrackLength.h"



class BohdanAna : public marlin::Processor, EventDisplayer {
    friend class EventDisplayer;
    public:
        BohdanAna(const BohdanAna&) = delete;
        BohdanAna& operator=(const BohdanAna&) = delete;

        marlin::Processor* newProcessor() { return new BohdanAna; }

        BohdanAna();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();

        void fillMCTrueInfo(EVENT::MCParticle* mc, EVENT::ReconstructedParticle* pfo, bool isInReconstructablePrimaryVertex, bool isInReconstructableSecondaryVertex);
        void fillTOFInfo(EVENT::Cluster* cluster, EVENT::Track* track, EVENT::CalorimeterHit* closestHit, const dd4hep::rec::Vector3D& trackPosAtCalo, const dd4hep::rec::Vector3D& trackMomAtCalo);
        void fillTrackLengthInfo(EVENT::ReconstructedParticle* pfo, EVENT::MCParticle* mc, const UTIL::LCRelationNavigator& navToSimTrackerHits);
        void fillCsvForKonrad(EVENT::Cluster* cluster, int pdg, double trueTOF, double trackLength, const dd4hep::rec::Vector3D& trackPosAtCalo, const dd4hep::rec::Vector3D& trackMomAtCalo);
        void fillTrackStates(EVENT::ReconstructedParticle* pfo, EVENT::ReconstructedParticle* refittedPfo);


        void resetVariables();
    private:
        int _nEvent{};
        float _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;

        bool _produce_csv_output{};
        std::ofstream _csv_output_file;
        int _global_pfo_number{};

        int _quarksToPythia{}; // event information!
        // True infromation of MCParticle
        int _pdg{};
        std::array<float, 3> _mcVtx{};
        std::array<float, 3> _mcMom{};
        float _d0True{};
        float _z0True{};
        float _omegaTrue{};
        float _tanLambdaTrue{};
        float _timeTrue{};
        bool _isOverlay{};
        bool _isSimulated{};
        int _imidiateParentPDG{};
        int _firstStableParentPDG{};
        bool _isBottomQuarkDecay{};
        bool _isCharmQuarkDecay{};
        bool _isHadronisationDecay{};
        bool _isInReconstructablePrimaryVertex{};
        bool _isInReconstructableSecondaryVertex{};
        bool _isReconstructed{};
        bool _hasTrack{};
        bool _hasShower{};

        //True vertex information
        int _nTracksInVertexTrue{};
        bool _isV0True{};
        bool _oppositeChargeTrue{};
        float _invMassTrue{};
        float _minEnergyTrue{};
        float _cosThetaTrue{};
        float _cosThetaIPTrue{};

        float _dEdx{};
        float _omegaIP{};
        float _omegaECAL{};
        float _tanLambdaIP{};
        float _tanLambdaECAL{};

        std::array<float, 3> _recoIpMom{};
        std::array<float, 3> _recoCaloPos{};
        std::array<float, 3> _recoCaloMom{};

        float _refittedOmegaIP{};
        float _refittedOmegaECAL{};
        float _refittedTanLambdaIP{};
        float _refittedTanLambdaECAL{};
        std::array<float, 3> _refittedRecoIpMom{};
        std::array<float, 3> _refittedRecoCaloPos{};
        std::array<float, 3> _refittedRecoCaloMom{};


        float _trackLength_IDR{};
        float _trackLength_SHA_phiLambda_IP{};
        float _trackLength_SHA_phiZed_IP{};
        float _trackLength_SHA_zedLambda_IP{};
        float _trackLength_SHA_phiLambda_ECAL{};
        float _trackLength_SHA_phiZed_ECAL{};
        float _trackLength_SHA_zedLambda_ECAL{};

        float _trackLength_IKF_phiLambda{};
        float _trackLength_IKF_phiZed{};
        float _trackLength_IKF_zedLambda{};
        float _harmonicMom_IKF_phiLambda{};
        float _harmonicMom_IKF_phiZed{};
        float _harmonicMom_IKF_zedLambda{};

        float _trackLengthToSET_IKF_phiLambda{};
        float _trackLengthToSET_IKF_phiZed{};
        float _trackLengthToSET_IKF_zedLambda{};
        float _harmonicMomToSET_IKF_phiLambda{};
        float _harmonicMomToSET_IKF_phiZed{};
        float _harmonicMomToSET_IKF_zedLambda{};

        bool _cleanTrack = true;

        //each element per time resolution from 0 to 100 ps
        int _typeClosest = -1;
        int _caloIDClosest = -1;
        int _layoutClosest = -1;
        int _layerClosest = -1;
        bool _cleanClosestHit = false;
        std::vector<float> _resolutions = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300};
        std::array<float, 15> _tofClosest{};
        std::array<float, 15> _tofAverage{};
        std::array<float, 15> _tofSETFront{};
        std::array<float, 15> _tofSETBack{};
        std::array<float, 15> _tofFit{};

        int _nHits;
        std::vector<float> _xHit;
        std::vector<float> _yHit;
        std::vector<float> _zHit;
        std::vector<float> _tHit;
        std::vector<int> _layerHit;
        std::vector<float> _energyHit;

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
        TApplication _application = TApplication("app", 0, nullptr);
};

#endif
