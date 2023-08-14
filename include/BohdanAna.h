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

        std::array<double, 3> _mcMom{};
        std::array<double, 3> _recoIpMom{};
        std::array<double, 3> _recoCaloMom{};

        double _trackLength_SHA_phiLambda_IP{};
        double _trackLength_SHA_phiZed_IP{};
        double _trackLength_SHA_zedLambda_IP{};
        double _trackLength_SHA_phiLambda_ECAL{};
        double _trackLength_SHA_phiZed_ECAL{};
        double _trackLength_SHA_zedLambda_ECAL{};
        TrackLengthResult _trackLength_IKF_phiLambda;
        TrackLengthResult _trackLength_IKF_phiZed;
        TrackLengthResult _trackLength_IKF_zedLambda;


        //each element per time resolution from 0 to 100 ps
        int _layerClosest = -1;
        std::array<double, 11> _tofClosest{};
        std::array<double, 11> _tofAverage{};
        std::array<double, 11> _tofSET{};
        std::array<double, 11> _tofFit{};

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
        TApplication _application = TApplication("app", 0, nullptr);
};

#endif
