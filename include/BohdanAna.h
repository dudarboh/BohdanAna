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
        double _bField{};
        
        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;

        int _pdg{};
        double _dEdx{};
        double _omegaIP{};
        double _omegaECAL{};
        double _tanLambdaIP{};
        double _tanLambdaECAL{};

        std::array<double, 3> _mcMom{};
        std::array<double, 3> _recoIpMom{};
        std::array<double, 3> _recoCaloPos{};
        std::array<double, 3> _recoCaloMom{};


        double _trackLength_IDR{};
        double _trackLength_SHA_phiLambda_IP{};
        double _trackLength_SHA_phiZed_IP{};
        double _trackLength_SHA_zedLambda_IP{};
        double _trackLength_SHA_phiLambda_ECAL{};
        double _trackLength_SHA_phiZed_ECAL{};
        double _trackLength_SHA_zedLambda_ECAL{};
        TrackLengthResult _trackLength_IKF_phiLambda;
        TrackLengthResult _trackLength_IKF_phiZed;
        TrackLengthResult _trackLength_IKF_zedLambda;
        bool _cleanTrack = true;

        //each element per time resolution from 0 to 100 ps
        int _typeClosest = -1;
        int _caloIDClosest = -1;
        int _layoutClosest = -1;
        int _layerClosest = -1;
        bool _cleanClosestHit = false;
        std::vector<double> _resolutions = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300};
        std::array<double, 15> _tofClosest{};
        std::array<double, 15> _tofAverage{};
        std::array<double, 15> _tofSET{};
        std::array<double, 15> _tofFit{};

        int _nHits;
        std::vector<double> _xHit;
        std::vector<double> _yHit;
        std::vector<double> _zHit;
        std::vector<double> _tHit;
        std::vector<int> _layerHit;
        std::vector<double> _energyHit;

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
        TApplication _application = TApplication("app", 0, nullptr);
};

#endif
