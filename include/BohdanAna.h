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

        void resetVariables();
    private:
        int _nEvent{};
        float _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;

        int _global_pfo_number{};
        bool _produce_csv_output = true;
        std::ofstream _csv_output_file = std::ofstream("output.csv");

        int _pdg{};
        float _dEdx{};
        float _omegaIP{};
        float _omegaECAL{};
        float _tanLambdaIP{};
        float _tanLambdaECAL{};

        std::array<float, 3> _mcMom{};
        std::array<float, 3> _recoIpMom{};
        std::array<float, 3> _recoCaloPos{};
        std::array<float, 3> _recoCaloMom{};


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
