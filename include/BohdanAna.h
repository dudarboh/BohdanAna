#ifndef BohdanAna_h
#define BohdanAna_h 1

#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TFile.h"
#include "TTree.h"
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

        TrackLengthResult _trackLength; // only latest and greatest for now, with HM momenta

        //each element per time resolution from 0 to 100 ps
        int _layerClosest = -1;
        std::array<double, 11> _tofClosest{};
        std::array<double, 11> _tofAverage{};
        std::array<double, 11> _tofSET{};
        std::array<double, 11> _tofFit{};

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;
};

#endif
